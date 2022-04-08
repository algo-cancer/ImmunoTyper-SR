from .common import log, fasta_from_seq
from collections import Sequence
import pysam

class BwaWrapper(object):
    """ CLI wrapper for bwa mapper.

    Attr:
        params (str): 			saved bwa parameters to be used with run
        src (str): 				path / alias to bwa executable
        output_path (str):		mapping save location for caching / not using temp files
    """
    params = ''			## command arguments
    src = 'bwa' 	## path / command to access executabke
    output_path =  None


    def __init__(self, params=None, bwa_src='bwa', output_path=None):
        from subprocess import Popen, PIPE

        self.params = params
        self.src = bwa_src
        if output_path: self.output_path = output_path


    def __str__(self):
        return(self.src + self.params)

    @staticmethod
    def create_temp_file(write_data=None, delete=True):
        import tempfile
        result = tempfile.NamedTemporaryFile(delete=delete, mode="w")
        if write_data:
            result.write(write_data)
            result.flush()
        return result

    def map(self, query, target, src=None, params=None, parser=None, output_path=None):
        """Runs bwa using subprocess.

        Args:
            query: 				iterable of read-like objects or path to fasta
            target:				path to fasta. Must be in same directory as index.
            src (str)			path to bwa executable. self.src if None
            params (str):		string of bwa parameters. self.params if None
            parser (func(x)):	parser func for bwa stdout result. bwaWrapper.paf_parser if None
            output_path (str):	cache path to save mapping result to
        
        Note:
            read-like requires 'id' and 'seq' attributes

        Returns:
            output: 			result of parser
        """

        from subprocess import Popen, PIPE
        from os.path import exists
        
        ## Check type(query), make temp file and write query seqs as needed
        query_file = None
        if isinstance(query, str):
            if not exists(query):
                log.error('Provided query path is invalid, please provide a path as a string or Bio.SeqIO-like objects')
            query_path = query
        else:
            query_file = self.create_temp_file(write_data=fasta_from_seq(*zip(*[(x.id, x.seq) for x in query])))
            query_path = query_file.name

        ## Check type(target), make temp file and write target seqs as needed
        if not (isinstance(target, str) and exists(target)):
            log.error('Provided target path is invalid, please provide a path as a string')
        target_path = target

        if not src:
            src = self.src
        if not params:
            params = self.params
        if not output_path:
            output_path = self.output_path
        if not parser:
            parser = self.sam_parser	

        # make output file if needed
        if not output_path:
            output_file = self.create_temp_file()
            output_path = output_file.name

        # run commmand
        command = ' '.join([src, 'mem', params, target_path, query_path, '>', output_path])
        log.debug('Running bwa mem:\n{}'.format(command))
        print('Running bwa mem:\n{}'.format(command))
        # process = Popen(command.split(), stdout=PIPE, stderr=PIPE, shell=True)
        # stdout, stderr = process.communicate()
        
        # log.debug(stdout)
        # log.debug(stderr)
        # print(stdout)
        # print(stderr)

        import os
        os.system(command)


        if query_file:
            query_file.close()

        return parser(output_path)
        

    @staticmethod
    def sam_parser(output_path):
        file = pysam.AlignmentFile(output_path, 'rb')
        return file