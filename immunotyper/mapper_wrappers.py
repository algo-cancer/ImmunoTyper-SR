from os.path import exists
import os, abc, pysam
from src.common import log, fasta_from_seq

class MappingWrapper():
    """ Abstract classe for mapping tool wrappers.

    Attr:
        params (str): 			saved bwa parameters to be used with run
        src (str): 				path / alias to bwa executable
        output_path (str):		mapping save location for caching / not using temp files
    """
    __metaclass__ = abc.ABCMeta

    params = ''			## command arguments
    src = '' 	## path / command to access executabke
    output_path =  None

    def __init__(self, src=None, params=None, output_path=None):
        self.params = params
        if src: self.src = src
        if output_path: self.output_path = output_path

    # def __str__(self):
    #     return(self.src + self.params)

    @staticmethod
    def create_temp_file(write_data=None):
        import tempfile
        result = tempfile.NamedTemporaryFile(delete=False, mode="w")
        if write_data:
            result.seek(0)
            result.write(write_data)
            result.truncate()
            result.flush()
        return result
    
    @staticmethod
    def get_primary_mapping(mapping_list):
        '''Takes a list of pysam.AlignedSegment objects, returns the one that is primary mapping
        Raises ValueError if it does not have,or has multiple, primary mappings'''
        primary_mapping = [x for x in mapping_list if not x.is_secondary and not x.is_supplementary]
        if len(primary_mapping) > 1:
            raise ValueError('More than 1 primary mapping:\n{}'.format(str([str(x) for x in mapping_list])))
        if len(primary_mapping) == 0:
            raise ValueError('No primary mappings:\n{}'.format(str([str(x) for x in mapping_list])))
        return primary_mapping[0]

    def map(self, query, target, src=None, params=None, parser=None, output_path=None, *args, **kwargs):
        """Runs src using os.system. Depends on tool-specific implementation of self.build_command()

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
        if not isinstance(target, str):
            log.error('Provided target path is invalid, please provide a path as a string')
            raise ValueError('Provided target path is invalid, please provide a path as a string')
        if not  exists(target):
            log.error('Provided target file does not exist: %s ' % target)
            raise ValueError('Provided target file does not exist: %s ' % target)

        target_path = target

        if not src:
            src = self.src
        if not params:
            params = self.params if self.params else ''
        if not output_path:
            output_path = self.output_path
        if not parser:
            parser = self.sam_parser	

        # make output file if needed
        if not output_path:
            output_file = self.create_temp_file()
            output_path = output_file.name
        else:
            output_file = None

        # run commmand
        command = self.build_command(src, params, query_path, target_path, output_path, *args, **kwargs)

        log.debug('Running {}:\n{}'.format(self.src, command))
        print('Running {}:\n{}'.format(self.src, command))
        os.system(command)

        if query_file:
            query_file.close()

        result = parser(output_path)
        if output_file:
            output_file.close()
        return result
        
    @abc.abstractmethod
    def build_command(self, src, params, query_path, target_path, output_path, *args, **kwargs):
        raise NotImplementedError

    @staticmethod
    def sam_parser(output_path):
        file = pysam.AlignmentFile(output_path, 'rb')
        return file
    
class BwaWrapper(MappingWrapper):
    src = 'bwa'
    def build_command(self, src, params, query_path, target_path, output_path):
        return ' '.join([src, 'mem', params, target_path, query_path, '>', output_path])

class BowtieWrapper(MappingWrapper):
    src = 'bowtie2'
    def build_command(self, src, params, query_path, target_path, output_path):
        return ' '.join([src, params, '-x', target_path, '-U', query_path, '-S', output_path])

class MrsFast(MappingWrapper):
    src = 'mrsfast'
    def build_command(self, src, params, query_path, target_path, output_path):
        return ' '.join([src, '--search', target_path, '--seq', query_path, '-o', output_path, params])
    
class BwaWrapperBiowulf(BwaWrapper):
    def build_command(self, src, params, query_path, target_path, output_path):
        module_laod = 'module load bwa && '
        return module_laod + super().build_command(src, params, query_path, target_path, output_path)

class BowtieWrapperBiowulf(BowtieWrapper):
    def build_command(self, src, params, query_path, target_path, output_path):
        module_laod = 'module load bowtie/2 && '
        return module_laod + super().build_command(src, params, query_path, target_path, output_path)
