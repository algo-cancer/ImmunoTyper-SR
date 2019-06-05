import os, sys, tempfile, abc
from common import log, fasta_from_seq, SeqRecord
from Bio import SeqIO
import blasr_wrapper


#---------------------------------------------#
# Abstract base class for allele calling tool #
#---------------------------------------------#
class AlleleCaller(object):
	__metaclass__ = abc.ABCMeta

	consensus_builder = None

	@abc.abstractmethod
	def call_cluster(self, cluster, filter_function=None):
		## Generates an allele call (or list of calls) for cluster, filtered by filter_function
		##		cluster can be:
		##			1. cluster_class.Cluster instance. Then calls are added to instance attributes and nothing is returned
		##			2. List of Bio.SeqRecord.SeqRecord-like objects. Then call(s) is returned, if multiple as a list
		##			3. A Bio.Seq.Seq-like object. This is assumed to be a consensus sequence, calls are made for this sequence and returns as in 2.
		##			4. A string path to a consensus sequence. As per 3.
		## Filter function takes a call / candidate and return True if it a true call and False otherwise. Ideally has a description of the filter via __str__ implementation (like NameFunction)
		raise NotImplementedError


class MinimapCaller(AlleleCaller):

	def __init__(self, minimap=None, allele_db=None):
		if allele_db:
			self.allele_db = allele_db
		else:
			import create_allele_db
			self.allele_db = '../database/V-QUEST-reference-allele-db+no-period-references.clustalw.fasta'	# TODO implement either fetching the default database file (instead of hard coded) or make a getter that creates a temporary fasta of the allele_db loaded from create allele_db
		log.info('Using the following allele database:\n' + self.allele_db)

		if minimap:
			self.minimap = minimap
		else:
			import minimap_wrapper
			self.minimap = minimap_wrapper.minimap_wrapper()


	def call_cluster(self, consensus_seq):
		## Takes a list of string sequences and returns a consensus sequence
		
		with tempfile.NamedTemporaryFile(delete=True) as f:
			f.write(fasta_from_seq('consensus', consensus_seq))

			command = [self.minimap.src, '-cx map-pb', f.name, self.allele_db]
 
			mapping_output = self.minimap.run(command)


class BlasrCaller(AlleleCaller):

	result_filter= staticmethod(lambda x: x.qName.split('|')[1]) #function on blasr_wrapper.BlasrMapping instance to return desired value (usually resulting allele name)
	description = ['Allele database: ', 'Blasr command: ', 'Filter Function: ']
	allele_db = None
	blasr = None
	filter_function = None
	consensus_builder = None



	def __init__(self, blasr=blasr_wrapper.BlasrWrapper(), consensus_builder=None, allele_db=None, filter_function=None, result_filter=None):
		if allele_db:
			self.allele_db = allele_db
		else:
			# import create_allele_db
			self.allele_db = '../database/Complete.Human.IGHV_IMGT.Feb2018.Corey.linear.modified.fasta'	# TODO implement either fetching the default database file (instead of hard coded) or make a getter that creates a temporary fasta of the allele_db loaded from create allele_db
		log.info('Using the following allele database:\n' + self.allele_db)
		self.description.append('Allele database: '+self.allele_db)

		if blasr:
			self.blasr = blasr
		else:
			import blasr_wrapper
			self.blasr = blasr_wrapper.BlasrWrapper()
		if filter_function:
			log.info('Setting blasr result filter function')
			self.filter_function = filter_function
		else:
			self.filter_function = lambda y: y.score	# key that returns selection criteria to choose the best hit for the blasr mapping of consensus to allele database. Results will be sorted in ascending order using this key
		if consensus_builder:
			self.consensus_builder = consensus_builder
	
	def __str__(self):
		result = 'Blasr-based allele caller:\n'
		return result + '\n'.join(
			[descr + str(value) for descr, value in 
				zip(self.description, [self.allele_db, self.blasr, self.filter_function])
			if value])

	def call_cluster(self, cluster, filter_function=None, result_filter=None, temp_file_path=None):
		import tempfile

		if len(cluster) == 1:
			log.warn('Cluster {} has single read, not calling'.format(cluster.id))
			try:
				cluster.consensus_seq=None
				cluster.consensus_builder=None
				cluster.set_call(None)
				cluster.candidates = None
				cluster.candidates_method = str(self)
			except AttributeError as e:
				pass
			finally: return None


		consensus_seq = None
		consensus_seq_id = None
		f=None
		is_cluster_inst = False						# flag for filling descriptive attributes
		if hasattr(cluster, '__getitem__'):			# assumed to be list of sequences, get consensus
			try:
				if temp_file_path:
					with open(temp_file_path, 'wb') as f:
						f.write(fasta_from_seq(*zip(*[(x.id, x.seq) for x in cluster])))
				consensus_seq = self.consensus_builder.generate_consensus(temp_file_path if temp_file_path else cluster)
				if not consensus_seq:
					cluster.consensus = None
					cluster.candidates_method = str(self)
					return
				consensus_seq_id = 'cons'
				log.info('Generated consensus with:\n{}'.format(str(self.consensus_builder)))
				log.debug('Output:\n{}'.format(consensus_seq))

				try:
					cluster.consensus = consensus_seq
					cluster.consensus_method = str(self.consensus_builder)
				except AttributeError as e:
					pass
			except TypeError as e:					## No consensus builder is set
				raise ValueError('Cluster calling: list of cluster sequences provided but no consensus builder instantiated.')
		else:		
			if isinstance(cluster, basestring):		# input is path
				if os.path.exists(cluster):
					cons_path = cluster
				else:
					raise ValueError('Cluster calling input invalid. String provided but is not valid path. If trying to cast as Bio.Seq.Seq-like object')
			else:									# input is consensus seq
				consensus_seq = cluster.seq
				consensus_seq_id = cluster.id
		
		## save blasr target in all cases except path as input
		if consensus_seq:
			try:
				f = open(temp_file_path, 'wb+') if temp_file_path else tempfile.NamedTemporaryFile(delete=False)
				f.write(str(fasta_from_seq(consensus_seq_id, consensus_seq)))
				cons_path=f.name
				f.close()
			except AttributeError as e:
				raise ValueError('Cluster calling input invalid. Provide iterable of cluster sequences, path to cluster consensus or Bio.Seq.Seq-like object to call')


		## run blasr mapping of consensus_seq against allele database
		command = [self.blasr.src, '', self.allele_db, cons_path]
			
		try:
			mapping_output = self.blasr.run(*command)
		except ValueError as e:
			log.warn('Blasr returned no mapping')
			try:
				cluster.set_call(None)
				cluster.candidates = None
				cluster.candidates_method = str(self)
			except AttributeError as e:
				pass
			finally: return None

		f.close()

		## select from mapping the desired result as the call
		if not filter_function:
			filter_function = self.filter_function
		
		try:
			mapping_output = sorted(mapping_output, key=filter_function)
			cluster_call = mapping_output[0]
		except ValueError as e:
			log.error('Invalid blasr mapping value')
			log.debug('\n'.join([str(x) for x in mapping_output]))
			raise e

		if not result_filter:
			result_filter = self.result_filter
		result = result_filter(cluster_call)

		try:
			cluster.set_call([result])
			cluster.candidates = list(mapping_output)
			cluster.candidates_method = str(self)
		except AttributeError as e:
			return result