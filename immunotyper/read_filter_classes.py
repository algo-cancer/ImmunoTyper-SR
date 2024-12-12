import abc, os, Bio, sys, tempfile
from collections import defaultdict
from .mapper_wrappers import BwaWrapper, BowtieWrapper
from abc import abstractmethod
from os.path import splitext, exists, join
from .common import SeqRecord, log, Read, fasta_from_seq, ReadFlanking

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
    
    def make_read(self, *args, **kwargs):
        return Read(*args, **kwargs)

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
            log.debug(f'Filtering {len(self.reads)} reads')
        self.positive, self.negative = self.filter(self.map(input_reads))
        log.debug(f'Of {len(self.reads)}, {len(self.positive)} passed and {len(self.negative)} failed filter')

        return self.positive, self.negative

    
    def map(self, reads):
        '''Args
        reads       Can be str path to fasta/q, or list of Bio.SeqRecord-like objects'''

        mapping_results = None
        if self.load_cache_path:
            try:
                log.debug('Trying to load mapping results cache from {}'.format(self.load_cache_path))
                mapping_results = self.mapper.sam_parser(self.load_cache_path)
            except (ValueError, IOError):
                log.debug('Loading failed!')

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
                    self.reads[mapping.query_name] = self.make_read(mapping.query_name, mapping.query_sequence)
                    self.reads[mapping.query_name].mappings = []
                    self.reads[mapping.query_name].add_mapping(mapping)
                if not mapping.is_secondary and not mapping.is_supplementary:
                    self.reads[mapping.query_name].primary_mapping = mapping
            else:
                try:
                    self.reads[mapping.query_name].primary_mapping = None
                except KeyError as e:
                    self.reads[mapping.query_name] = self.make_read(mapping.query_name, mapping.query_sequence)
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
        positive = set()
        negative = set()

        for read in reads:
            try:
                if self.read_is_positive(read, key(read)):
                    positive.add(read)
                else:
                    negative.add(read)
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
    

    def __init__(self, mapper, mapping_params=None, reference_path=None, write_cache_path=None, load_cache_path=None, minimum_coding_bases=50, secondary_filter=True):
        self.perform_secondary_filter = secondary_filter
        self.minimum_coding_bases = minimum_coding_bases
        log.debug('Using {} as minimum number of coding base threshold for mapping filter'.format(self.minimum_coding_bases))
        log.debug('Allowing clipped mappings to start within {} bases of a reference start or end'.format(self.read_end_tolerance))
        self.secondary_filter_reads = dict()
        self.secondary_filter_mappings = defaultdict(dict)
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
            if m.reference_start == 0 or m.reference_end == self.get_reference_length(m):    # alignment is at end of reference
                return True
            else:
                self.secondary_filter_reads[read.id] = read
                self.secondary_filter_mappings[m.query_name][m.reference_name] = m
        return False
    
    def filter(self, reads, key=lambda r: r.mappings):

        positive, negative = super().filter(reads, key)
        
        if self.secondary_filter_reads and self.perform_secondary_filter:
            secondary_positive = self.secondary_clipping_filter(self.secondary_filter_reads)
            log.debug(f'Of {len(self.secondary_filter_reads)}, {len(secondary_positive)} passed secondary filter')
            positive.update(secondary_positive)
            negative = negative - secondary_positive
        
        return list(positive), list(negative)
    
    secondary_params = '-a -L 1000'
    def secondary_clipping_filter(self, reads_dict):
        '''Filters on reads_dict of reads with mapping with a suspicious soft clipping, 
        and remaps with query and target flipped and with a larger soft clipping tolerance to 
        force extension to the end of th reference'''

        if not reads_dict:
            return set()
        
        log.info(f'Performing secondary clipping filter on {len(reads_dict)} reads')

        tempdir =  tempfile.TemporaryDirectory()

        reads_fasta = open(os.path.join(tempdir.name, 'reads.fasta'), 'w')
        reads_fasta.write(fasta_from_seq(*zip(*[(x.id, x.seq) for x in reads_dict.values()])))
        self.mapper.index_reference(reads_fasta.name)
        self.secondary_mappings = list(self.mapper.map(self.reference_path, reads_fasta.name, params=self.secondary_params, output_path=os.path.join(tempdir.name, 'mappings.sam')))
        
        positive = set()
        for m in self.secondary_mappings:
            if m.is_unmapped:
                continue
            swapped_m = SwappedAlignmentSegment(m)
            r = reads_dict[m.reference_name]
            try:
                r.secondary_mappings.append(swapped_m)
            except AttributeError:
                r.secondary_mappings = [swapped_m]
            if self.mapping_is_positive(r, swapped_m) and (not m.reference_name in self.secondary_filter_mappings[r.id] or m.query_alignment_length > self.secondary_filter_mappings[r.id][m.reference_name].query_alignment_length):
                if r.primary_mapping:
                    swapped_m.is_secondary = True
                r.add_mapping(swapped_m)
                positive.add(r)
        
        tempdir.cleanup()

        return positive


    def get_reference_length(self, mapping):
        if isinstance(mapping, SwappedAlignmentSegment):
            return mapping.mapping.infer_read_length()
        return self.mapping_results.lengths[mapping.reference_id]

    @property
    def maximum_start(self):
        return self.read_end_tolerance
    
    def minimum_end(self, mapping):
        return self.get_reference_length(mapping)
    
class BwaFlankingFilter(FlankingFilter):
    def __init__(self, mapping_params=None, reference_path=None, write_cache_path=None, load_cache_path=None, minimum_coding_bases=50):
        super(BwaFlankingFilter, self).__init__(BwaWrapper(params='-a', output_path=write_cache_path), mapping_params, reference_path, write_cache_path, load_cache_path, minimum_coding_bases)

class BowtieFlankingFilter(FlankingFilter):
    def __init__(self, mapping_params=None, reference_path=None, write_cache_path=None, load_cache_path=None, minimum_coding_bases=50):
        super(BowtieFlankingFilter, self).__init__(BowtieWrapper(params='-a --end-to-end --very-sensitive  --n-ceil C,100,0 --np 0 --ignore-quals --mp 2,2 --score-min C,-50,0 -L 10', output_path=write_cache_path), mapping_params, reference_path, write_cache_path, load_cache_path, minimum_coding_bases)


class FlankingDatabaseFilter(FlankingFilter):
    def __init__(self, allele_db, *args, **kwargs):
        self.allele_db = allele_db
        super().__init__(*args, **kwargs)
    
    def mapping_is_positive(self, read, m):
        '''Returns true if mapping m is fully aligned or has an alignment of AT LEAST self.minimum_coding_bases at the beginning or end of the reference
        Also sets read.end_mappings (set) as the mappings that are clipped at beginning/end'''
        if not m or m.is_unmapped:
            return False

        allele_id = read.reference_id_parser(m.reference_name)
        if not self.allele_db[allele_id].has_flanking:
            return super().mapping_is_positive(read, m)
        
        if (m.reference_end <= self.allele_db[allele_id].coding_end) and (m.reference_start >= self.allele_db[allele_id].coding_start) and m.query_alignment_length == len(read):   # mapping is not clipped and wholely in coding sequence
            return True

        if m.reference_start <= self.allele_db[allele_id].coding_start:      # mapping is upstream
            if m.reference_end < self.allele_db[allele_id].coding_start:       # no coding bases
                return False
            if len(self.allele_db[allele_id].upstream_flanking) == 0:   # no flanking sequence
                return super().mapping_is_positive(read, m)
            if m.query_alignment_length == len(read):       # unclipped
                return True               
            if m.reference_start <= self.maximum_start and m.query_alignment_length >= self.minimum_coding_bases:     # clipped at start of sequence
                return True
        
        if m.reference_end >= self.allele_db[allele_id].coding_end:              # mapping is downstream
            if m.reference_start > self.allele_db[allele_id].coding_end:       # no coding bases
                return False
            if len(self.allele_db[allele_id].downstream_flanking) == 0:   # no flanking sequence
                return super().mapping_is_positive(read, m)
            if m.query_alignment_length == len(read):       # unclipped
                return True               
            if m.reference_end >= self.minimum_end(m) and m.query_alignment_length >= self.minimum_coding_bases:     # clipped at end of sequence
                return True
            
        return False

   
    def make_read(self, *args, **kwargs):
        return ReadFlanking(self.allele_db, *args, **kwargs)

class BwaFlankingDatabaseFilter(FlankingDatabaseFilter):
    def __init__(self, allele_db, mapping_params=None, reference_path=None, write_cache_path=None, load_cache_path=None, minimum_coding_bases=50):
        super().__init__(allele_db, BwaWrapper(params='-a', output_path=write_cache_path), mapping_params, reference_path, write_cache_path, load_cache_path, minimum_coding_bases)


class BowtieFlankingDatabaseFilter(FlankingDatabaseFilter):
    def __init__(self, allele_db, mapping_params=None, reference_path=None, write_cache_path=None, load_cache_path=None, minimum_coding_bases=50):
        super().__init__(allele_db, BowtieWrapper(params='-a --end-to-end --very-sensitive -f  --n-ceil C,100,0 --np 0 --ignore-quals --mp 2,2 --score-min C,-50,0 -L 10', output_path=write_cache_path), mapping_params, reference_path, write_cache_path, load_cache_path, minimum_coding_bases)


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


class SwappedAlignmentSegment(object):
    '''Wrapper for pysam.AlignedSegment that swaps the query and reference sequences and start and end values'''

    def __init__(self, mapping):
        self.mapping = mapping
        self.get_tag = mapping.get_tag
        self.reference_sequence = mapping.query_sequence
        self.query_name = mapping.reference_name
        self.reference_name = mapping.query_name
        self.query_alignment_start = mapping.reference_start
        self.query_alignment_end = mapping.reference_end
        self.reference_start = mapping.query_alignment_start
        self.reference_end = mapping.query_alignment_end
        self.query_alignment_length = mapping.reference_length
        self.reference_length = mapping.query_length
        self.is_reverse = mapping.is_reverse
        self.is_secondary = mapping.is_secondary
        self.is_supplementary = mapping.is_supplementary
        self.is_unmapped = mapping.is_unmapped
        self.mapping_quality = mapping.mapping_quality
        self.cigarstring = mapping.cigarstring
        self.cigartuples = mapping.cigartuples
        self.tags = mapping.tags
        self.template_length = mapping.template_length
        self.query_qualities = mapping.query_qualities
        self.query_alignment_qualities = mapping.query_qualities
        self.query_sequence_length = mapping.reference_length
        self.query_qualities = mapping.query_qualities
        self.query_length = mapping.reference_length
        self.query_name = mapping.reference_name
        self.reference_end = mapping.query_alignment_end
        self.reference_length = mapping.query_length
        self.reference_name = mapping.query_name
        self.reference_start = mapping.query_alignment_start
