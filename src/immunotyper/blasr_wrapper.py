from common import log, fasta_from_seq
import os

class BlasrWrapper(object):

	##### Class Attributes #####
	params = ''			## command arguments
	src = 'blasr' 	## path / command to access executabke
	output_path =  None

	##### Class Methods ######
	def __init__(self, params=None, blasr_src=None, output_path=None):
		from subprocess import Popen, PIPE

		if params: self.params = params
		if blasr_src: self.src = blasr_src
		if output_path: self.output_path = output_path


	def __str__(self):
		return(' '.join([self.src, self.params]))


	@staticmethod
	def run(src, params, query, target, output_path=None, parser=None):
		## Just runs the command using subprocess
		## parser must be a function that handles the raw blasr stdout string
		from os.path import exists
		def create_temp_file(write_data=None):
			import tempfile
			result = tempfile.NamedTemporaryFile(delete=False)
			if write_data:
				result.write(write_data)
			return result

		from subprocess import Popen, PIPE

		if not parser:
			parser = BlasrWrapper.parse_psl

		if isinstance(query, basestring):
			if not exists(query):
				log.error('Provided blasr query path is invalid, please provide a path as a string or Bio.SeqIO-like objects')
			query_file = None
		else:
			try:
				query_file = create_temp_file(fasta_from_seq(*zip(*[(x.id, x.seq) for x in query])))
				query = query_file.name
			except AttributeError as e:
				log.error('Provided blasr query input is invalid, please provide a path as a string or Bio.SeqIO-like objects')
				raise e

		if isinstance(target, basestring):
			if not exists(target):
				log.error('Provided blasr target path is invalid, please provide a path as a string or Bio.SeqIO-like objects')
			target_file = None
		else:
			try:
				target_file = create_temp_file(fasta_from_seq(*zip(*[(x.id, x.seq) for x in target])))
				target = target_file.name
			except AttributeError as e:
				log.error('Provided blasr target input is invalid, please provide a path as a string or Bio.SeqIO-like objects')
				raise e

		if params:
			command = [src, '-m', str(5), '--header'] + params + [query, target]
		else:
			command = [src, '-m', str(5), '--header', query, target]

		if query_file:
			query_file.close()
		if target_file:
			target_file.close()

		log.debug('Running blasr:\n{}'.format(' '.join(map(str, command))))

		process = Popen([str(x).strip() for x in command], stdout=PIPE, stderr=PIPE)
		stdout, stderr = process.communicate()

		if 'is empty' in stderr:
			log.error(stderr)
			with open(query, 'r') as f:
				log.debug('{}\n{}'.format(query, f.read()))
			with open(target, 'r') as f:
				log.debug('{}\n{}'.format(target, f.read()))
			raise ValueError(stdout)

		try:
			if not stdout.split('\n'):
				raise ValueError
			output = parser(stdout.strip())
			if not output:
				raise ValueError
		except ValueError as e:
			log.warn('Blasr return no mapping')
			log.debug(stderr)
			log.debug(stdout)
			raise e

		if output_path:
			try:
				with open(output_path, 'wb') as f:
					f.write(stdout)
			except OSError as e:
				log.error('Provided blasr output path is not valid, output will be discarded')

		return output			

	@staticmethod
	def parse_psl(mapping_output):
		result = []
		header = None
		for line in mapping_output.split('\n'):
			if 'qName' in line:
				header = line.split()
				continue
			if 'blasr' in  line:
				continue
			result.append(BlasrMapping(line) if not header else BlasrMapping(line, header))
		return result


class BlasrMapping():
	def __init__(self, data_line, header=['qName', 'qLength', 'qStart', 'qEnd', 'qStrand', 'tName', 'tLength', 'tStart', 'tEnd', 'tStrand', 'score', 'numMatch', 'numMismatch', 'numIns', 'numDel', 'mapQV', 'qAlignedSeq', 'matchPattern', 'tAlignedSeq']):
		if isinstance(header, str):
			header = header.split()
		self.header = header
		self.attributes = []
		isint = set(['qLength', 'qStart', 'qEnd', 'tLength', 'tStart', 'tEnd', 'score', 'numMatch', 'numMismatch', 'numIns', 'numDel', 'mapQV'])

		try:
			for attribute, value in zip(header, data_line.split()):
				if attribute in isint:
					exec('self.{}=int("{}")'.format(attribute, value))
				else:
					exec('self.{}="{}"'.format(attribute, value))
				self.attributes.append(value)
		except IndexError as e:
			log.error('Blasr output invalid or header data fields invalid')
			log.debug('Line:\n' + data_line)
			log.debug('Header data fields: \n' + ' '.join(header))
			raise e

	def __str__(self):
		return " ".join(self.attributes[:-4])
		
	def get_header(self):
		return " ".join(self.header[:-4])