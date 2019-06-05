from common import log, fasta_from_seq
from itertools import izip_longest
from collections import Sequence
class MinimapWrapper(object):
	""" CLI wrapper for Minimap mapper (github.com/lh3/minimap2).

	Attr:
		params (str): 			saved minimap parameters to be used with run
		src (str): 				path / alias to minimap executable
		output_path (str):		mapping save location for caching / not using temp files
	"""
	params = '-cx map-pb'			## command arguments
	src = 'minimap' 	## path / command to access executabke
	output_path =  None


	def __init__(self, params=None, minimap_src=None, output_path=None):
		from subprocess import Popen, PIPE

		if params: self.params = params
		if minimap_src: self.src = minimap_src
		if output_path: self.output_path = output_path


	def __str__(self):
		return(self.src + self.params)

	def ava(self, reads=None, params=None):
		"""Runs minimap with query == target for all-to-all mapping

		Args:
			reads: 					iterable of read-like objects
			params:					mapping parameters. Overrides self.params
		Returns:
			result: 				result of self.run with params using self.paf_parser
		Raises:
			ValueError:				if 'ava-pb' not in params (or self.params if params not provided).
		"""
		if not params:
			params = self.params
		elif not params:
			params = 'ava-pb'
		elif 'ava' not in params:
			raise ValueError('''This MinimapWrapper has parameter configurations that do not included all-to-all mapping. 
								Please change the params attribute or manually overide the parameters when calling MinimapWrapper.ava''')

		reads_file = self.create_temp_file(fasta_from_seq(*zip(*[(x.id, x.seq) for x in reads])))

		return self.run(reads_file.name, reads_file.name, params=params)	## NOTE / TODO using tempfile.name twice MAY NOT work on some systems


	@staticmethod
	def create_temp_file(write_data=None):
		import tempfile
		result = tempfile.NamedTemporaryFile(delete=False)
		if write_data:
			result.write(write_data)
			result.flush()
		return result

	def run(self, query, target, src=None, params=None, parser=None, output_path=None):
		"""Runs minimap using subprocess.

		Args:
			query: 				iterable of read-like objects or path to fasta
			target:				iterable of read-like objects or path to fasta
			src (str)			path to minimap executable. self.src if None
			params (str):		string of minimap parameters. self.params if None
			parser (func(x)):	parser func for minimap stdout result. MinimapWrapper.paf_parser if None
			output_path (str):	cache path to save mapping result to
		
		Note:
			read-like requires 'id' and 'seq' attributes

		Returns:
			output: 			result of parser
		"""

		from subprocess import Popen, PIPE
		from os.path import exists
		
		## Check type(query), make temp file and write query seqs as needed
		if isinstance(query, basestring):
			if not exists(query):
				log.error('Provided query path is invalid, please provide a path as a string or Bio.SeqIO-like objects')
			query_path = query
			query_file = None
		else:
			# try:
			# try:
			query_file = self.create_temp_file(write_data=fasta_from_seq(*zip(*[(x.id, x.seq) for x in query])))
			# except TypeError: # is not iterable
				# query_file = self.create_temp_file(write_data=fasta_from_seq(query.id, query.seq))
			# except AttributeError as e:
			# 	log.error('Provided query input is invalid, please provide a path as a string or Bio.SeqIO-like objects')
				# raise e
			query_path = query_file.name

		## Check type(target), make temp file and write target seqs as needed
		if isinstance(target, basestring):
			if not exists(target):
				log.error('Provided target path is invalid, please provide a path as a string or Bio.SeqIO-like objects')
			target_path = target
			target_file = None
		else:
			try:
				try:
					target_file = self.create_temp_file(write_data=fasta_from_seq(*zip(*[(x.id, x.seq) for x in target])))
				except TypeError: # is not iterable
					target_file = self.create_temp_file(write_data=fasta_from_seq(target.id, target.seq))
			except AttributeError as e:
				log.error('Provided target input is invalid, please provide a path as a string or Bio.SeqIO-like objects')
				raise e
			target_path = target_file.name

		if not src:
			src = self.src
		if not params:
			params = self.params
		if not output_path:
			output_path = self.output_path
		if not parser:
			parser = MinimapWrapper.paf_parser	


		command = ' '.join([src, params, target_path, query_path])
		log.debug('Running minimap:\n{}'.format(command))
		process = Popen(command.split(), stdout=PIPE, stderr=PIPE)
		stdout, stderr = process.communicate()

		## save / cache output if needed
		if output_path:
			try:
				with open(output_path, 'wb') as f:
					f.write(stdout)
			except OSError as e:
				log.error('Provided minimap output path is not valid, output will be discarded')

		if not stdout.strip():
			log.error('Minimap returned no mapping')
			log.debug(stderr)
			log.debug(stdout)
			with open(query_path, 'r') as f:
				log.debug(f.readlines())
			with open(target_path, 'r') as f:
				log.debug(f.readlines())
			raise ValueError('Minimap returned no mapping')

		output = parser(stdout.strip())
		
		if query_file:
			query_file.close()
		if target_file:
			target_file.close()

		return output
		

	@staticmethod
	def paf_parser(output):
		output = output.strip().split('\n')
		result = []
		for i, line in enumerate(output):
			try:
				result.append(PAFMapping(line))
			except ValueError as e:
				log.error('Error in Minimap output line {}:\n{}'.format(i, line))
				raise e
		return result



class PAFMapping(Sequence):

	def __init__(self, data_line, header=['qName', 'qLength', 'qStart', 'qEnd', 'strand', 'tName', 'tLength', 'tStart', 'tEnd', 'numMatch', 'alignedLength', 'mapQV']):
		self.header = header
		self.attributes = []
		isint = set(['qLength', 'qStart', 'qEnd', 'tLength', 'tStart', 'tEnd', 'score', 'numMatch', 'alignedLength', 'mapQV'])
		try:
			data = data_line.split()
			for attribute, value in izip_longest(header, data):
				if attribute:
					if attribute in isint:
						exec('self.{}=int("{}")'.format(attribute, value))
					else:
						exec('self.{}="{}"'.format(attribute, value))
				elif value == data[-1]:
					exec('self.cigar="{}"'.format(value))
					self.header.append('Cigar')
				else:
					exec('self.{}="{}"'.format(value.split(':')[0], value))
					self.header.append(value.split(':')[0])
				self.attributes.append(value)
			self.NM = [x for x in data if 'NM' in x][0]
		except IndexError as e:
			log.error('PAF output invalid or header data fields invalid')
			log.debug('Line:\n' + str(data))
			log.debug('Header data fields: \n' + str(header))
			raise e


	def get_NM(self):
		return int(self.NM.split(':')[2])


	## Get error functions: if norm=True it is normalized relative to the sequence length
	def tErrors(self, norm=False):
		e = self.get_NM() + self.tStart + (self.tLength-self.tEnd)
		if norm:
			e = e / float(self.tLength)
		return e

	def qErrors(self, norm=False):	
		e = self.get_NM() + self.qStart + (self.qLength-self.qEnd)
		if norm:
			e = e / float(self.qLength)
		return e

	def total_errors(self, norm=False):
		e = self.tErrors() + self.qErrors()
		if norm:
			e = e / float(self.tLength + self.qLength)
		return e

	def set_quality(self, val):
		self.quality = val
		self.attributes.append(self.quality)
		self.header.append('Quality')

	def __str__(self):
		return " ".join([str(x) for x in self.attributes])

	def __getitem__(self, i):
		return self.attributes[i]
	def __len__(self):
		return self.alignedLength
		
	def get_header(self):
		return " ".join(self.header)