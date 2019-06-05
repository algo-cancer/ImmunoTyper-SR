from common import log, fasta_from_seq
from cluster_class import Cluster
from subprocess import Popen, PIPE, call
import abc, tempfile, os, subprocess
from Bio import SeqIO

#---------------------------------------------#
# Abstract base class for distance calulators #
#---------------------------------------------#
class DistanceCalculator(object):
	__metaclass__ = abc.ABCMeta
	matrix = None

	@abc.abstractmethod
	def generate_distances(self, reads):
		raise NotImplementedError

	##### Static Methods #####
	@staticmethod
	def filter_matrix(raw_matrix, filter_func):
		## Filters edges out based on some filtering criteria
		## Input: matrix = distance matrix
		##			filter_func = function that takes the edge distance, returns True if is passes the filter False otherwise
		## Output: returns a filtered COPY of the input matrix

		# make copy of matrix to not modify input
		raw_matrix_nodes = len(DistanceCalculator.get_nodes(raw_matrix))
		raw_matrix_edges = len(DistanceCalculator.get_edges(raw_matrix))

		result = {}

		log.info('Filtering distance matrix')
		log.debug('Before:\nSize dict: {}\tNumber nodes: {}\tNumber edges: {}'.format(len(raw_matrix), 
																					raw_matrix_nodes, 
																					raw_matrix_edges))

		for read, maps in raw_matrix.iteritems():
			for map_name, distance in maps.iteritems():
				if filter_func(distance):
					if read not in result:
						result[read] = {}
					result[read][map_name] = distance
		


		# logging
		result_nodes = len(DistanceCalculator.get_nodes(result))
		result_edges = len(DistanceCalculator.get_edges(result))

		log.debug('After:\nSize dict: {}\tNumber nodes: {}\tNumber edges: {}'.format(len(result), 
																					result_nodes, 
																					result_edges))
		log.info('''{} / {} reads removed due to all neighbouring / mapped reads 
				having distance that does not pass filter'''.format(raw_matrix_nodes-result_nodes, raw_matrix_nodes))
		log.info('{} / {} total edges removed by distance filter'.format(raw_matrix_edges-result_edges, raw_matrix_edges))

		return result

	@staticmethod
	def get_nodes(matrix):
		## Returns set of all unique nodes in provided distance matrix
		result = set(matrix.keys())
		for v in matrix.values():
			result = result.union(v.keys())
		return result
	@staticmethod
	def get_edges(matrix):
		## Returns set of all unique nodes in provided distance matrix
		result = []
		for n1, v in matrix.iteritems():
			result.extend([(n1, n2, distance) for n2, distance in v.iteritems()])
		return result

class EditDistanceCalculator(DistanceCalculator):
	filter_function = None

	def __init__(self, msa_generator):
		self.msa_generator = msa_generator

	def generate_distances(self, reads, filter_func=lambda x: True):
		## Generates distance matrix of form {read_id }	
		## reads = [Bio.SeqIO, ...] = list of ORIENTED (ie no rev-compl) reads to be clustered (ie containing genes)
			## if None uses self.reads
		## minimap = instance of MinimapWrapper object. If none uses pre-set self.minimap
		## filter_func = distance included in output if threshold(distance) = True
		import copy
		from Bio import SeqIO
		from itertools import combinations	

		if self.matrix:
			log.info('Using cached distance matrix')
			return self.matrix

		def distance(i, j):
			result = 0
			for index in range(0, len(i)):
				try:
					if i[index] != j[index]:
						result += 1
				except IndexError:
					raise IndexError('Length seq 1 ({}) != length seq 2 ({})\n{}\n{}'.format(len(i), len(j), i, j))
			return result

		if self.filter_function:
			filter_func = self.filter_function

		msa = self.msa_generator.generate_msa(reads)
		
		result = {}
		
		for seq1, seq2 in combinations(msa, 2):
			if seq1.id not in result:
				result[seq1.id] = {}
			result[seq1.id][seq2.id] = distance(seq1.seq, seq2.seq)

		if filter_func:
			result_filtered =  self.filter_matrix(copy.deepcopy(result), filter_func)
			
			result = result_filtered

		self.matrix = result
		return result



class MinimapDistanceCalculator(DistanceCalculator):

	##### Class Attributes #####
	reads = None
	minimap = None

	##### Class Methods #####
	def __init__(self, minimap=None, reads=None, filter_function=None):
		import minimap_wrapper

		self.minimap = minimap if minimap else minimap_wrapper.MinimapWrapper(params='-cx ava-pb -k14 -w3')

		self.filter_function = filter_function

	def generate_distances(self, reads=None, minimap=None, filter_func=lambda x: True):
		## Generates distance matrix of form {read_id }	
		## reads = [Bio.SeqIO, ...] = list of ORIENTED (ie no rev-compl) reads to be clustered (ie containing genes)
			## if None uses self.reads
		## minimap = instance of MinimapWrapper object. If none uses pre-set self.minimap
		## filter_func = distance included in output if threshold(distance) = True
		import copy
		from Bio import SeqIO

		if self.filter_function:
			filter_func = self.filter_function

		if self.matrix:
			log.info('Using cached distance matrix')
			result = self.matrix
			if filter_func:
				result_filtered =  self.filter_matrix(copy.deepcopy(self.matrix), filter_func)
				result = result_filtered
			return self.matrix
		
		if not minimap:
			minimap = self.minimap

		mapping = minimap.ava(reads=reads)

		result = {}
		
		mapped_reads = set() # for keeping track of mapped reads to report missing reads
		for i, line in enumerate(mapping):
					
			mapped_reads.add(line.qName)
			mapped_reads.add(line.tName)
			

			try:
				if 'NM' not in line.NM:
					raise IndexError
				NM = int(line.NM.split(':')[2])
			except IndexError as e:
				log.error('Error in Minimap output: NM field is likely missing\nmapping line:{}'.format('\t'.join(line)))
				log.debug(dir(line))
				log.debug(zip(line.header, line.attributes))
				raise ValueError()
			
			distance_value = (NM+line.qStart+(line.qLength-line.qEnd) + NM+line.tStart+(line.tLength-line.tEnd)) / float(line.qLength+line.tLength)
			
			if line.qName not in result:
				result[line.qName] = {}
			
			result[line.qName][line.tName] = distance_value
		
		# check if any reads missing from mapping
		missing = set([x.id for x in reads]).difference(mapped_reads)
		if missing:
			log.warn('{} / {} reads missing from mapping'.format(len(missing), len(list(reads))))
			log.debug('\n'.join(list(missing)))
		
		self.matrix = result
		
		if filter_func:
			result_filtered =  self.filter_matrix(copy.deepcopy(result), filter_func)
			
			result = result_filtered

		return result




#-----------------------------------------#
# Abstract base class for clustering tool #
#-----------------------------------------#
class ClusteringGenerator(object):
	__metaclass__ = abc.ABCMeta
	@abc.abstractmethod
	def cluster(self, reads):
		raise NotImplementedError


class DSF(ClusteringGenerator):


	##### Class attributes #####
	params = ''
	src = 'dsf'
	distance_calculator=None


	##### Class Methods #####
	def __init__(self, distance_calculator, params=None, dsf_src=None):
		## str(parameter) = command line arguments for dsf tool as list of strings
		## dsf_src = path to dsf executable as string
		from cluster_class import Cluster
		self.distance_calculator=distance_calculator
		self.matrix = distance_calculator.matrix
		if params:
			self.params = params
		if dsf_src:
			self.src = dsf_src
	def __str__(self):
		result = ['DSF clustering method']
		if self.distance_calculator:
			result.append('Distance calculator:\n'+str(self.distance_calculator))
		result.append('DSF command:\n' + ' '.join([self.src, self.params]))
		return '\n'.join(result)



	##### Static Methods ######
	@staticmethod
	def run(src, params, input_path, output_path, mappings=False):		  
		import sys, stat
		if not os.path.exists(src):
			log.critical('DSF executable not found at provided location: ' + src)
			raise IOError

		sys.path.append(os.path.split(src)[0])
		import dense_subgraph_finder


		command = [x for x in [src, params, '--graph', input_path, '-o', output_path] if x]

		with open('temp.sh', 'wb+') as script:
			script.write('#!/bin/bash\n'+' '.join(command))
			st = os.stat(script.name)
			os.chmod(script.name, st.st_mode | stat.S_IEXEC)
			log.debug('Script at: {}\n{}'.format(script.name, script.read()))
			process = subprocess.call('/bin/bash '+script.name, shell=True)
			process = Popen(' '.join(command), stdout=PIPE, stderr=PIPE, shell=True)
			stdout, stderr = process.communicate()

			print(stdout)
			print(stderr)


		log.info('Running Dense Subgraph Finder:\n{}'.format(' '.join(command)))

		try:
			result = DSF.parse_dsf_output(path=output_path+'/dense_subgraphs.txt', mappings=mappings)
			return result
		except IOError as e:
			raise Exception('DSF failed, see debug log')
	
	@staticmethod
	def reverse_mappings(mappings, reads=None):
		## reverses mapping dictionary from convert_adjacency_matrix. i.e. dictionary with line number as key, read id as value
		## reads = {read_id: read, ...}. If provided the values of the returned dictionary are not the read ids but the read object itself
		result = dict([(v, k) for k, v in mappings.iteritems()])
		if reads:
			result = dict([(k, reads[v]) for k, v in result.iteritems()])
		return result

	@staticmethod
	def parse_dsf_output(path, mappings=False):
		# Input: mappings = reverse of mapping as returned by convert_adjacency_matrix: form: {str(original_read_name): int(dsf_vertex_numerical_id), ...}
		#        path = path to dsf output
		# Output: Returns list of clusters, where each cluster is a list of read names (or numerical ids if mapping not provided)
		
		get_read_name = (lambda x: mappings[x]) if mappings else (lambda x: x)
		
		with open(path, 'r') as f:
			clusters = [int(x.strip()) for x in f.readlines()]
		result = [[]] * max(clusters)
		result = {}
		for read, c in enumerate(clusters):
			if c not in result:
				result[c] = []
			result[c].append(get_read_name(read))
		return result.values()
	
	def cluster(self, reads, save_input_path=None, output_dir=None, cached_output=None):
		import shutil
		def get_seq_obj(output):
			seq_mapping = dict([(x.id, x) for x in reads])
			output_seqs = map(lambda x: [seq_mapping[y] for y in x], output)
			return output_seqs


		try:	## reads are a path
			self.reads = dict([(x.id, x) for x in SeqIO.parse(reads, 'fasta')])
		except AttributeError as e: 	## reads is a list of SeqRecord-like objects
			self.reads = dict([(x.id, x) for x in reads])
		finally:
			if cached_output and self.distance_calculator.matrix:
				self.input_matrix, num_edges, mapping = self.convert_adjacency_matrix(self.distance_calculator.matrix)
				mapping = self.reverse_mappings(mapping, self.reads)
				return [Cluster(x, cluster_id=i, clustering_tool=self)
							for i, x in enumerate(sorted(self.parse_dsf_output(cached_output, mapping),
														key=lambda x: len(x), reverse=True))]
			self.distance_calculator.generate_distances(reads)
		try:
			reads.close()
		except AttributeError as e:
			pass

		self.input_matrix, num_edges, mapping = self.convert_adjacency_matrix(self.distance_calculator.matrix)
		mapping = self.reverse_mappings(mapping, self.reads)

		log.debug('Number of edges in input graph:' +str(num_edges))



		# write adjacency to file for dsf input
		def writer(matrix, num_edges):
			yield '{} {} 001\n'.format(len(matrix), num_edges)
			for neighbours in matrix:
				line = ' '.join([' '.join(map(str, (n+1, w))) for n, w in sorted(neighbours, key=lambda x: x[0])])
				yield '{}\n'.format(line)
		
		matrix_output_iterator = writer(self.input_matrix, num_edges)
		in_file = tempfile.NamedTemporaryFile(delete=True)

		try:
			if save_input_path:
				in_file = open(save_input_path, 'wb')
		except IOError as e:
			in_file = tempfile.NamedTemporaryFile(delete=True)
			log.warn('Provided DSF input matrix write path not valid, using temporary file')
		log.info('Saving dsf input file to'.format(in_file.name))
		in_file.writelines(matrix_output_iterator)
		in_file.flush()    


		# check provided output_dir is valid
		if output_dir:
			if not os.path.exists(output_dir):
				log.warn('Provied DSF output directory path not valid, using temporary directory')
				output_dir = None
			else:
				temp_dir = None
		

		# make temp output dir if no valid output dir provided
		if not output_dir:
			temp_dir = tempfile.mkdtemp()
			output_dir = temp_dir
			saved_umask = os.umask(0077)			# Ensure the file is read/write by the creator only


		# run DSF
		try:
			output = self.run(self.src, self.params, in_file.name, output_dir, mapping)
			# run dsf
		except Exception as e:	# This is just so the temp files get deleted in the case some previous unhandled exception gets raised
			raise e
		finally:
			if temp_dir:
				os.umask(saved_umask)
				shutil.rmtree(temp_dir)
			in_file.close()

		## generate instances of cluster_class.Cluster as result
		output = [Cluster(x, cluster_id=i, clustering_tool=self) for i, x in enumerate(sorted(output, key=lambda x: len(x), reverse=True))]
		return output

	def convert_adjacency_matrix(self, matrix):
		## Converts dictionary form of distance matric into a linked list representing the form used by dsf
		counter = 0
		llist = []
		mapping = {}
		num_edges = 0

		def new_vertex(name, counter):
			llist.append([])
			mapping[name] = counter
			counter += 1
			return counter

		for read, neighbours in matrix.iteritems():
			if read not in mapping:
				counter = new_vertex(read, counter)

			for n, weight in neighbours.iteritems():
				if n not in mapping:
					counter = new_vertex(n, counter)

				llist[mapping[read]].append((mapping[n], weight))
				llist[mapping[n]].append((mapping[read], weight))
				num_edges += 1

		return llist, num_edges, mapping



class MappingClustering(ClusteringGenerator):
	""" Creates clusters from read-to-allele mapping information.

	Inherits from ABC ClusteringGenerator.
	Creates high confidence clusters by using unambiguous 
	read-to-allele mappings based of the read.mapping attribute. 
	"""

	def __init__(self, expected_coverage, cluster_support_threshold=None):
		"""Inits MappingClustering class.

		Args:
			expected_coverage (int): expected read depth / coverage. 
				Used to filter out non-confident calls / clusters
			cluster_support_threshold (int): minimum number of reads 
				in a cluster assignment to be considered confident.
				If None, sets as 75% of expected_coverage

		"""
		self.expected_coverage = expected_coverage
		self.cluster_support_threshold = int(0.85*expected_coverage) if not cluster_support_threshold else cluster_support_threshold


	def cluster(self, reads, no_cov_estimation=False):
		"""Clusters reads based on their allele mappings

		Args:
			reads: iterable of reads with mapping attribute

		Returns:
			result (list): list of high confidence 
				cluster_class.Cluster objects
			ambiguous_reads: Reads that do no belong to one of the 
				above clusters
		"""

		## Filter ambigiuous reads
		ambiguous_reads = []
		clusters = {}

		for r in reads:
			sorted_mappings = sorted(r.mapping, key=self.get_mapping_errors)
			if not len(r.mapping) ==1 and self.has_ambiguous(r.mapping):
				ambiguous_reads.append(r)
			else:
				try:
					clusters[sorted_mappings[0].tName.split('|')[1]].append(r)
				except KeyError:
					clusters[sorted_mappings[0].tName.split('|')[1]] = [r]		



		# Filter reads mapped to ignored alleles
		discarded = 0
		with open('../database/ignored_alleles.txt', 'r') as f:
			ignored_alleles = [x.strip() for x in f.readlines()]
		for a in ignored_alleles:
			if a in clusters:
				discarded += len(clusters[a])
				del clusters[a]
		log.info('{} reads have unambiguous mappings to short and ignored alleles, discarding these reads'.format(discarded))

		# Calculate expected coverage
		if no_cov_estimation:
			coverage = self.expected_coverage
			log.info('Not performing coverage estimation')
		else:
			cluster_sizes = [len(v) for v in clusters.values() if len(v) > 0.5*self.expected_coverage]
			coverage = sorted(cluster_sizes)[len(cluster_sizes)/2]
			self.cluster_support_threshold = 0.85*coverage
			log.info('Calculated coverage as {}'.format(coverage))


		# Filter low-confidence clusters
		discarded = 0
		low_confidence = []
		for k, v in clusters.iteritems():
			if len(v) < self.cluster_support_threshold:
				ambiguous_reads.extend(v)
				low_confidence.append(k)
		for k in low_confidence: 
			discarded += len(clusters[k])
			del clusters[k]
		log.info('{} alleles have insufficient mapped reads to be confident, discarding these assignments'.format(discarded))
		
		# Create cluster class intances and return
		result = []
		for num, (call, c_reads) in enumerate(clusters.iteritems()):
			c = Cluster(reads=c_reads, 
									cluster_id=num,
									clustering_tool=self)
			c.set_call([call]*int(round(float(len(c_reads)) / coverage)))
			result.append(c)

		return result, ambiguous_reads, coverage


	def has_ambiguous(self, mappings):
		"""Identifies if a read's mappings are ambiguous.
				
		Read is ambiguous if the difference between the lowest and
		second-lowest error mappings is 6 or less errors. Errors
		are defined by self.get_mapping_errors.

		Args:
			mappings: iterable of read-to-allele mapping objects with 
				attribute required by self.get_mapping_errors

		Returns:
			True if ambiguous else False
		"""
		errors = sorted([self.get_mapping_errors(x) for x in mappings])
		
		if errors[1]-errors[0] <=6:
			return True
		else: 
			return False

	@staticmethod
	def get_mapping_errors(mapping):
		"""Gets mapping.tErrors

		Allows for different mapping object types through inheritance /
		override.

		Args:
			mapping: minimap_wrapper.PAFMapping-like object or with tErrors() method

		Returns:
			int: mapping.tErrors value

		Raises:
			AttributeError: mapping does not have tErrors() implemented
		"""
		return mapping.tErrors()
			
