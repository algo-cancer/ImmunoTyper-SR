import abc, os, Bio, sys
from .mapper_wrappers import BwaWrapper, BowtieWrapper
from abc import abstractmethod
from os.path import splitext, exists, join
from .common import SeqRecord, log, Read

class Filter(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def filter_reads(self, reads, *agrs, **kwargs):
        ''' Takes list/generator of SeqRecord-like objects, returns [positive reads], [negative reads]'''
        raise NotImplementedError

class MappingFilter(Filter):
    __metaclass__ = abc.ABCMeta
    mapper = None           # instance of MappingWrapper class
    mapping_params = None
    reference_path = None
    write_cache_path = None
    load_cache_path = None

    def __init__(self, mapper, mapping_params=None, reference_path=None, write_cache_path=None, load_cache_path=None):
        self.mapper = mapper
        self.mapping_params = mapping_params
        self.reference_path = reference_path
        self.write_cache_path = write_cache_path
        self.load_cache_path = load_cache_path

    def filter_reads(self, reads, *agrs, **kwargs):
        if isinstance(reads, str):
            if not exists(reads):
                log.error('Provided query path is invalid, please provide a path as a string or Bio.SeqIO-like objects')
                raise ValueError('Provided query path is invalid, please provide a path as a string or Bio.SeqIO-like objects')
            input_reads = reads
            self.reads = {}
        else:
            self.reads = dict([(r.id, r) for r in reads])
            input_reads = self.reads.values()
        if self.reads:
            log.info(f'Filtering {len(self.reads)} reads')
        self.positive, self.negative = self.filter(self.map(input_reads))
        log.info(f'Of {len(self.reads)}, {len(self.positive)} passed and {len(self.negative)} failed filter')

        return self.positive, self.negative

    
    def map(self, reads):
        '''Args
        reads       Can be str path to fasta/q, or list of Bio.SeqRecord-like objects'''

        mapping_results = None
        if self.load_cache_path:
            try:
                log.info('Trying to load mapping results cache from {}'.format(self.load_cache_path))
                mapping_results = self.mapper.sam_parser(self.load_cache_path)
            except (ValueError, IOError):
                log.info('Loading failed!')

        if not mapping_results:
            if self.write_cache_path:
                log.info('Saving mapping results cache to {}'.format(os.path.join(self.write_cache_path)))
            mapping_results = self.mapper.map(reads, self.reference_path, 
                                                    params=self.mapping_params,
                                                    output_path=self.write_cache_path)
        return self.parse_mapping(reads, mapping_results)

    def parse_mapping(self, reads, mapping_results):
        self.mapping_results = mapping_results
                
        for mapping in self.mapping_results:
            if not mapping.is_unmapped:
                try:
                    self.reads[mapping.query_name].add_mapping(mapping)
                except KeyError:
                    self.reads[mapping.query_name] = SeqRecord(mapping.query_name, mapping.query_sequence)
                    self.reads[mapping.query_name].mappings = []
                    self.reads[mapping.query_name].add_mapping(mapping)
                if not mapping.is_secondary and not mapping.is_supplementary:
                    self.reads[mapping.query_name].primary_mapping = mapping
            else:
                try:
                    self.reads[mapping.query_name].primary_mapping = None
                except KeyError as e:
                    self.reads[mapping.query_name] = SeqRecord(mapping.query_name, mapping.query_sequence)
                    self.reads[mapping.query_name].mappings = []
                    self.reads[mapping.query_name].primary_mapping = None

        try:
            no_mapping = [k for k, v in reads.items() if not v.mappings]
            if no_mapping:
                log.info('%i reads did not appear in mapping results or were unmapped' % len(no_mapping))
                log.debug('\n'.join(no_mapping))
                for r in no_mapping:
                    reads[r].primary_mapping = None
        except AttributeError:
            pass    

        return list(self.reads.values())


    def filter(self, reads, key=lambda r: r.mappings):
        '''Filters reads using self.read_is_positive
        Args
            reads               List of Bio.SeqRecord-like objects
            key                 Function that returns the mappings when given a read
        Returns
            positive            List of reads that pass filter
            negative            List of reads that do not pass filter
        '''
        positive = []
        negative = []

        for read in reads:
            try:
                if self.read_is_positive(read, key(read)):
                    positive.append(read)
                else:
                    negative.append(read)
            except Exception as e:
                log.error('\nError with read %s' % read.id)
                raise
                
        return positive, negative
    
    @abstractmethod
    def read_is_positive(self, read, mappings):
        '''Takes a read and a list of pysam.AlignedSegment objects representing all mappings for the given read, 
        Returns True if read passes criteria determining it is an IGHV read,
        False otherwise'''
        raise NotImplementedError

    @staticmethod
    def get_primary_mapping(mapping_list):
        '''Takes a list of pysam.AlignedSegment objects, returns the one that is primary mapping
        Raises ValueError if it does not have,or has multiple, primary mappings'''
        primary_mapping = [x for x in mapping_list if not x.is_secondary and not x.is_supplementary]
        if len(primary_mapping) > 1:
            s = "\n".join([str(x)+f'\t Reference = {x.reference_name}' for x in primary_mapping])
            log.error(f'More than 1 primary mapping:\n{s}')
            raise ValueError('More than 1 primary mapping')
        if len(mapping_list) and not len(primary_mapping):
            raise ValueError('Provided mappings have no primary mapping'.format(str([str(x) for x in mapping_list])))
        return None if not primary_mapping else primary_mapping[0]
    

class FlankingFilter(MappingFilter):
    '''Remove all mappings with clipped alignments that have < minimum_coding_bases or clipped alignment is not at an end of reference.
    Flags read as False for filter() if all mappings removed.
    Sets read.raw_mappings '''
    minimum_coding_bases=50
    read_end_tolerance = 3         # when considering if an alignment is at the end of a reference (for clipped mappings), True if within this number of bases from the start/end

    def __init__(self, mapper, mapping_params=None, reference_path=None, write_cache_path=None, load_cache_path=None, minimum_coding_bases=50):
        self.minimum_coding_bases = minimum_coding_bases
        log.info('Using {} as minimum number of coding base threshold for mapping filter'.format(self.minimum_coding_bases))
        log.info('Allowing clipped mappings to start within {} bases of a reference start or end'.format(self.read_end_tolerance))
        super(FlankingFilter, self).__init__(mapper, mapping_params, reference_path, write_cache_path, load_cache_path)

    def read_is_positive(self, read, mappings):
        '''Returns true if primary mapping matches self.mapping_is_positive'''
        read.raw_mappings = read.mappings
        read.mappings = []
        for m in mappings:
            if self.mapping_is_positive(read, m):
                read.add_mapping(m)
        return True if read.mappings else False

    def mapping_is_positive(self, read, m):
        '''Returns true if mapping m is fully aligned or has an alignment of AT LEAST self.minimum_coding_bases at the beginning or end of the reference
        Also sets read.end_mappings (set) as the mappings that are clipped at beginning/end'''
        if not m or m.is_unmapped:
            return False
        if m.query_alignment_length == len(read):   # mapping is not clipped
            return True
        if m.query_alignment_length > self.minimum_coding_bases:     # mapping has sufficient aligned bases
            if (m.reference_start <= self.maximum_start) or (m.reference_end >= self.minimum_end(m)):    # alignment is at end of reference
                try:
                    read.end_mappings.add(m)
                except AttributeError:
                    read.end_mappings = set([m])
                return True
            else:
                return False

    def get_reference_length(self, mapping):
        return self.mapping_results.lengths[mapping.reference_id]

    @property
    def maximum_start(self):
        return self.read_end_tolerance
    
    def minimum_end(self, mapping):
        return self.get_reference_length(mapping) - self.read_end_tolerance

class BwaFlankingFilter(FlankingFilter):
    def __init__(self, mapping_params=None, reference_path=None, write_cache_path=None, load_cache_path=None, minimum_coding_bases=50):
        super(BwaFlankingFilter, self).__init__(BwaWrapper(params='-a', output_path=write_cache_path), mapping_params, reference_path, write_cache_path, load_cache_path, minimum_coding_bases)

class BowtieFlankingFilter(FlankingFilter):
    def __init__(self, mapping_params=None, reference_path=None, write_cache_path=None, load_cache_path=None, minimum_coding_bases=50):
        super(BowtieFlankingFilter, self).__init__(BowtieWrapper(params='-a --very-sensitive-local -f', output_path=write_cache_path), mapping_params, reference_path, write_cache_path, load_cache_path, minimum_coding_bases)


class TtnMappingFilter(FlankingFilter):
    '''Modification of FlankingFilter to filter TTN reads. Changes functionality to take a subset of a WGS->GRCh mapping
    instead of performing the mapping'''

    def __init__(self, region_chrom, region_start, region_end, mapper=None, mapping_params=None, reference_path=None, write_cache_path=None, load_cache_path=None, minimum_coding_bases=50):
        super().__init__(mapper=None, mapping_params=None, reference_path=None, write_cache_path=None, load_cache_path=None, minimum_coding_bases=50)
        self.region_chrom, self.region_start, self.region_end = region_chrom, region_start, region_end
    
    @property
    def maximum_start(self):
        return self.read_end_tolerance + self.region_start
    
    def minimum_end(self, mapping):
        return self.region_end - self.read_end_tolerance
