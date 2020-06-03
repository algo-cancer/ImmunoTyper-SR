import os, sys, tempfile, abc
from subprocess import Popen, PIPE
from common import log, fasta_from_seq
from Bio import SeqIO
#-----------------------#
# Abstract base classes #
#-----------------------#

class ConsensusBuilder(object):
	""" ABC for consensus building implementations. Inherit this to ensure compatibility.
	"""
	__metaclass__ = abc.ABCMeta
	@abc.abstractmethod
	def generate_consensus(self, reads):
		"""
		Args:
			sequences: 			iterable of read-like objects
		Returns:
			string of consensus sequence
		"""
		raise NotImplementedError


class MSAGenerator(object):
	""" ABC for consensus building implementations. Inherit this to ensure compatibility.
	"""
	__metaclass__ = abc.ABCMeta
	@abc.abstractmethod
	def generate_msa(self, reads):
		"""
		Args:
			reads: 			iterable of read sequences as str
		Returns:
			msa (dict):		dict of {msa_seq_is: read-like object containing sequence, ...}. Gaps are '-'
		"""
		raise NotImplementedError


class ExternalTool(object):
	""" ABC essentially providing CLI wrapper functions when using exteranl tool for MSA / consensus generation.  

	Attr:
		params (list): 			CLI parameters for run, of form ['-p1', 'val1', '-p2', 'val2', ...]
		src (str): 				path / alias to executable
		name (str):				Name of tool (for logging)
	"""
	__metaclass__ = abc.ABCMeta

	params = []
	src = ''
	name=''

	def __init__(self, params=None, src=None):
		if params:
			self.params = params
		if src:
			self.src = src

	def __str__(self):
		return self.src + ' ' + (self.params if isinstance(self.params, basestring) else ' '.join([str(x) for x in list(self.params)]))

	def get_name(self):
		if not self.name:
			raise NotImplementedError
		else:
			return self.name

	def setup_run(self, cluster, params=None):
		"""Creates temp file with cluster reads and provides CLI-ready parameter list, using self.params 

		Args:
			cluster: 				iterable of read-like objects or path to reads
			params:					Overrides self.params
		Returns:
			run_input_path: 		path to input file
			params (str):			str of params for run
		"""

		if not params: params = self.params

		run_input_path = None
		if isinstance(cluster, str):        # input path provided
			if not os.path.exists(cluster):
				raise ValueError('Input to {} is not valid, please provide either a list of read objects or a valid path to fasta file'.format(self.get_name()))
			else:
				run_input_path = cluster
		else:       						# list of seq objects provided
			try:
				f = tempfile.NamedTemporaryFile(suffix='.fa', delete=False)
				f.write(fasta_from_seq(*zip(*[(x.id, x.seq) for x in cluster])))
				run_input_path = f.name
			except AttributeError as e:
				log.error('Input to {} is not valid, please provide either a list of read objects or a valid path to fasta file'.format(self.get_name()))
				raise e
		return run_input_path, params

	def run(self, params, input_path, parser, output_path=None):
		"""Runs the external tool 

		Args:
			parser (func):			function that takes stdout of run. If output to stdout is not possible, the path to the output needs to be provided as output_path
			output_path: 			path to output, for parser
		Returns:	output of parser
		"""

		if not os.path.exists(self.src):
			log.critical('Executable not found at provided location: ' + self.src)
			raise IOError
		command = [self.src]
		command.extend(params)
		if input_path: command.append(input_path)

		# log.debug('Running {}:\n{}'.format(self.get_name(), ' '.join(command)))

		process = Popen(command, stdout=PIPE, stderr=PIPE)
		stdout, stderr = process.communicate()

		try: 
			return parser(stdout.strip() if not output_path else output_path)
		except (IndexError, IOError, ValueError):
			log.error('{} failed'.format(self.get_name()))
			log.debug(stderr)
			log.debug(stdout)
			return None


class Spoa(ConsensusBuilder, MSAGenerator, ExternalTool):

	def parse_spoa_msa(self, data):
		return [x for x in data.split('\n') if 'Multiple sequence alignment' not in x]

	def parse_spoa_consensus(self, data):
		result = data.split('\n')[1].lower()
		if not result:
			raise ValueError
		return result


	##### Class attributes #####
	params = ''
	src = '../tools/spoa/build/bin/spoa'
	name='SPOA'

	##### Class Methods #####

	def generate_consensus(self, cluster, params=None):
		
		cluster, params = self.setup_run(cluster, params)		## TODO figure out how to close the file object

		try:
			output = self.run(['-r 0'], cluster, parser=self.parse_spoa_consensus)
		except ValueError as e:
			log.error('Consensus generation failed for cluster {}'.format(cluster.id))
		return output


	def generate_msa(self, cluster, params=None):
		cluster, params = self.setup_run(cluster, params)		## TODO figure out how to close the file object
		try:
			output = self.run(['-r 1'], cluster, parser=self.parse_spoa_msa)
		except ValueError as e:
			raise e

		# cluster_read_ids = [x.id for x in SeqIO.parse(cluster, 'fasta')]

		# result = [SeqRecord(read_id, seq) for read_id, seq in zip(cluster_read_ids, output)]
		return output

class Poa(ConsensusBuilder, MSAGenerator, ExternalTool):

	def parse_poa_pir(self, data, cluster=None):
		
		output = list(SeqIO.parse(data, 'fasta'))
		consensus = [x for x in output if 'CONSENS' in x.id]
		msa = dict([(x.id, str(x.seq)) for x in output if 'CONSENS' not in x.id])


		try:
			cluster.msa = []
			# check PIR contains same reads as cluster
			if len(set([c.id for c in cluster]).symmetric_difference(msa.keys())) > 0:
				raise IOError
			# assign MSA sequences to reads
			for read in cluster:
				read.alignment = str(msa[read.id])
				cluster.msa.append(read.alignment)
		except AttributeError as e:
			pass

		if len(consensus) > 1:
			try:
				log.warn('{} generated alignment for cluster {} returned {} consensus sequences - using only first sequence:'.format(self.name, cluster.id, len(consensus)))
			except AttributeError as e:
				log.warn('{} generated alignment returned {} consensus sequences - using only first sequence'.format(self.name, len(consensus)))
			log.debug('\n'.join([x.id for x in consensus]))
		
		if consensus:
			consensus = str(consensus[0].seq)
		
		try:
			cluster.consensus = str(consensus)
		except AttributeError:
			pass

		return consensus, msa



	##### Class attributes #####
	params = '-read_fasta {} -hb -preserve_seqorder -do_global -pir {} ../tools/poaV2/blosum80.mat'
	src = '../tools/poaV2/poa'
	name='POA'

	##### Class Methods #####

	def generate_consensus_msa(self, cluster, params=None, output_path=None):
		if not params:
			params = self.params
		
		input_path, params = self.setup_run(cluster, params)		## TODO figure out how to close the file object
		if not output_path:
			output_file = tempfile.NamedTemporaryFile(suffix='.fa', delete=False)
			output_path = output_file.name

		parser = lambda x: self.parse_poa_pir(x, cluster)
		try:
			consensus, msa = self.run(params=params.format(input_path, output_path).split(' '), 
							input_path=None,
							output_path=output_path,
							parser=parser)
		except ValueError as e:
			raise e
		return consensus, msa


	def generate_consensus(self, cluster, params=None):
		return self.generate_consensus_msa(cluster, params)[0]

	def generate_msa(self, cluster, params=None, output_path=None):
		result = self.generate_consensus_msa(cluster, params, output_path)
		return result[1]


