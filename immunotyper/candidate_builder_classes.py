from __future__ import annotations
from .common import log, create_temp_file, suppress_stdout, Read, fasta_from_seq
import abc, pysam, os
from abc import abstractmethod
from .mapper_wrappers import BwaWrapper, BowtieWrapper
from collections import OrderedDict, defaultdict
from tempfile import TemporaryDirectory

from typing import TYPE_CHECKING, List
if TYPE_CHECKING:
    from .models import ModelMeta

class Landmark(object):
    def __init__(self, allele_db_landmark):
        self.pos, self.exp_cov, self.stdevs = allele_db_landmark
        self.reads = []

class Candidate(object):
    def __init__(self, name, allele_db):
        self.id = name
        self.allele_db = allele_db

        self.reads = []
        self.reads_dict = {}

        self.is_ignored = allele_db[self.id].is_ignored
        self.is_functional = allele_db[self.id].is_functional
    
    def add_read(self, read_candidate):
        self.reads.append(read_candidate)
        self.reads_dict[read_candidate.id] = read_candidate
        read_candidate.candidate = self

    @staticmethod
    def reference_id_parser(raw_id):
        return raw_id.split('|')[1]
    
    @abstractmethod
    def num_copies(self) -> int:
        """Returns number of called copies

        Returns:
            int: Number of called copies
        """
        return NotImplementedError

    def __eq__(self, __o: object) -> bool:
        return (self.id == __o.id)
    
    def __hash__(self) -> int:
        return hash(self.id)

class LandmarksNotSetError(Exception):
    """Raise instead of key error when provided allele is not in database"""
    pass


class LandmarkCandidate(Candidate):
    """Candidate for ILP model
    """

    def __init__(self, name, allele_db):
        super().__init__(name, allele_db)
        try:
            self.landmark_groups = [[Landmark(landmark) for landmark in group] for group in allele_db[self.id].landmark_groups] #[read]
            self.landmarks = [l for g in self.landmark_groups for l in g]
        except KeyError:
            raise KeyError('Follow mapping has reference name of {} not in allele database\n{}'.format(self.id, str(self.mapping)))
        except (AttributeError, TypeError):
            raise LandmarksNotSetError(f"Landmarks are not set for {self.id}")

    @property
    def num_copies(self) -> int:
        return int(self.model.get_value(self.copy_multipler))
    
    @property
    def model(self):
        ''''''
        try:
            return self._model
        except AttributeError as e:
            raise NotImplementedError(f"Model attribute in Candidate class needs to be set by model instance: {self.id}  {str(e)}")
    
    @model.setter
    def model(self, val: ModelMeta):
        self._model = val

    def add_read(self, read_candidate):
        super().add_read(read_candidate)
        for landmark in self.landmarks:
            if read_candidate.covers_position(landmark.pos, self.id):
                landmark.reads.append(read_candidate)

# Classes for EM Model


class NovelVariant(object):
    """Small class for each novel variant
    """
    def __init__(self, candidate: Candidate, pos: int, value: str, ref_value: str) -> None:
        #TODO doctring

        self.candidate = candidate
        self.pos = pos
        self.value = value.upper()
        self.ref_value = ref_value.upper()
        self.is_wildtype = True if self.ref_value == self.value else False

        self.make_depth_expectations()
        
        self.reads = []
        self.read_positions = {} # {read_id (str): position in read.seq of variant (int)}
        
    @property
    def read_support(self):
        return len(self.reads)

    @property
    def num_copies(self) -> int:
        return int(self.candidate.model.get_value(self.copy_multipler))

    def make_depth_expectations(self):
        """Makes and sets self.exp_cov and self.stdevs using self.candidate.allele_db - see Landmark object
        """
        _, self.exp_cov, self.stdevs = self.candidate.allele_db.calculate_landmark_expection(self.candidate.allele_db[self.candidate.id], self.pos)


    def add_read(self, read: Read, variant_position: int):
        self.reads.append(read)
        self.read_positions[read.id] = variant_position
        
        # add variant to read object for reverse lookup
        read.add_variant(self)
    
    def __str__(self) -> str:
        return f"{self.candidate.id} variant at pos {self.pos}, ref base {self.ref_value}, alt base {self.value}"
    
    def __eq__(self, __o: object) -> bool:
        return (self.candidate == __o.candidate) and (self.pos == __o.pos) and (self.value == __o.value) and (self.ref_value == __o.ref_value)
    
    def __hash__(self) -> str:
        return hash(str(self))

            

class FlankingCandidate(LandmarkCandidate):
    """Uses mapper, pysam.pileup to determine alignment errors isolated to coding region
    saving them as attributes
    """

    def __init__(self, name, allele_db, mapper=BowtieWrapper(params='--end-to-end --very-sensitive --n-ceil C,100,0 --np 0 --ignore-quals --mp 2,2 --score-min C,-50,0 -L 10',
                                                            output_sorted_bam=True)):
        super().__init__(name, allele_db)
        
        # set novel-variant parameters
        self.mapper = mapper
        try:
            self.seq = allele_db[self.id].seq_with_flanking
        except AttributeError:
            self.seq = allele_db[self.id].seq
        self.allele = allele_db[self.id]
    
            
    def calculate_coding_distances(self, cache_dir=None) -> None:
        
        if not self.allele.has_flanking:
            raise ValueError(f"Allele {self.id} does not have flanking sequence")

        if not self.reads: raise TypeError(f"No reads assigned to {self.id}")
        # build mapping output tempfiles and map with mapper
        if not cache_dir:
            allele_ref_seq_dir_obj = TemporaryDirectory(delete=True)
            allele_ref_seq_dir = allele_ref_seq_dir_obj.name
        else:
            allele_ref_seq_dir = cache_dir
        fasta_path = os.path.join(allele_ref_seq_dir, self.id.replace('/', '#').replace('(', '-').replace(')', '-')+'.fa')
        bam_safe_allele_id = self.id.replace('(', '-').replace(')', '-')
        bam_file_path = fasta_path.replace('.fa', '.bam')
        bam_file_index = bam_file_path.replace('.bam', '.bai')
        
        try: # load cache
            mapping_file = pysam.AlignmentFile(bam_file_path)
        except (FileNotFoundError, ValueError): # make mapping
            
            with suppress_stdout():
                # write fasta with self seq as only entry, bwa index
                with open(fasta_path, 'w') as f:
                    f.write(fasta_from_seq(bam_safe_allele_id, self.allele.seq_with_flanking))
                self.mapper.index_reference(fasta_path)
                try:
                    # map to self reads
                    mapping_file = self.mapper.map(self.reads, fasta_path,
                                            output_path=bam_file_path)
                except ValueError: # no reads map
                    log.error(f"Invalid BAM mapping for allele candidate {self.id}")
                    return
            
            pysam.index(bam_file_path, bam_file_index)

        # Iterate through pileup of bam
        for pileupcolumn in mapping_file.pileup(): # position
            if pileupcolumn.reference_pos < self.allele.coding_start or pileupcolumn.reference_pos > self.allele.coding_end:  # pileupcolumn reference position is not coding position
                continue
            
            # sanity check
            if self.reads[0].reference_id_parser(pileupcolumn.reference_name) != bam_safe_allele_id:
                raise ValueError(f"Novel variant mapping for {self.id} has pileup to {pileupcolumn.reference_name}")

            for pileupread in pileupcolumn.pileups: # reads piled on position
                if pileupread.is_del or self.seq[pileupcolumn.reference_pos] != pileupread.alignment.query_sequence[pileupread.query_position]: # ref_base_value != read_base_value
                    try:
                        self.reads_dict[pileupread.alignment.query_name].coding_distances[self.id] += 1
                    except KeyError:
                        self.reads_dict[pileupread.alignment.query_name].coding_distances[self.id] = 1
                    except AttributeError:
                        self.reads_dict[pileupread.alignment.query_name].coding_distances = {self.id: 1}
                    
        if not cache_dir:
            allele_ref_seq_dir.cleanup() # remove all temp files
        

class CandidateBuilder(object):
    __metaclass__ = abc.ABCMeta
    candidates = {}         # index of candidates
    no_candidates = []      # reads that do not have any candidates

    def __init__(self, read_length, allele_db, candidate_class=LandmarkCandidate):
        self.read_length = read_length
        self.allele_db = allele_db
        self.candidate_class = candidate_class
    
    def make_candidates(self, reads):
        log.debug('Adding candidates to reads...')
        self.input_reads = reads
        self.candidates = {}
        for read in reads:
            read.candidates = []

            read_assignments = self.make_read_assignments(read)

            for allele in read_assignments:
                try:
                    candidate = self.candidates[allele]
                except KeyError:
                    try:
                        candidate = self.candidate_class(allele, self.allele_db)
                        self.candidates[allele] = candidate
                    except AttributeError as e:
                        log.warn(f'Unable to make candidate {allele} (length={len(self.allele_db[allele])})\n{str(e)}')
                        raise
                    except LandmarksNotSetError:
                        log.warn(f"Landmarks not set for {allele}, skipping...")
                        continue
                candidate.add_read(read)
                
                try:
                    read.candidates.append(candidate)
                except AttributeError:
                    read.candidates = [candidate]

        log.debug(f"Made {len(self.candidates)} candidates with {len([rc for c in self.candidates.values() for rc in c.reads])} read-allele pairs")

        return list(self.candidates.values())

    @abstractmethod
    def make_read_assignments(self, read):
        '''Takes a read object, return list of ReadAssignment objects for that read'''
        raise NotImplementedError


class BwaMappingsCandidateBuilder(CandidateBuilder):

    def make_read_assignments(self, read):
        try:
            return [read.reference_id_parser(m.reference_name) for m in read.mappings]
        except RecursionError:
            print(read)
            raise


class FlankingCandidateBuilder(BwaMappingsCandidateBuilder):
    def __init__(self, read_length, allele_db, candidate_class=FlankingCandidate):
        #TODO docstring
        super().__init__(read_length, allele_db, candidate_class)

    def make_candidates(self, reads, cache_dir: str = None):
        candidates = super().make_candidates(reads)
        print(f"Making coding distances for {len(candidates)} candidates...")
        
        for c in candidates:
            # try:
            c.calculate_coding_distances(cache_dir)
            # except AttributeError:
            #     pass

        return candidates




class LandmarkFeasibleCandidateFilter(object):
    def __init__(self, allele_db, num_landmarks, num_landmark_groups, expected_coverage, max_coverage_variation, read_length, minimum_coding_bases):
        self.num_landmarks = num_landmarks
        self.num_landmark_groups = num_landmark_groups
        self.read_length = read_length
        self.minimum_coding_bases = minimum_coding_bases
        self.expected_coverage = expected_coverage
        self.max_coverage_variation = max_coverage_variation
        self.allele_db = allele_db

        if self.allele_db.num_landmarks != num_landmarks:
            self.allele_db.make_landmarks(num_landmarks, read_length, expected_coverage, minimum_coding_bases, num_landmark_groups)
    
    def filter_feasible(self, candidates):
        '''Filters out Candidate instances that have insufficient depth coverage at landmark positions to be a feasible allele in the ILP'''
        log.info(f'Filtering candidate alleles with infeasible landmark coverage support...')

        self.positive = []
        self.negative = []
        for candidate in candidates:
            is_pos = True
            for group in candidate.landmark_groups:
                total_error = 0
                for landmark in group:
                    # total_error += max(0, landmark.exp_cov-len(landmark.reads))/(float(landmark.exp_cov)/len(candidate.landmark_groups))
                    total_error += max(0, landmark.exp_cov-len(landmark.reads))/(float(landmark.exp_cov))
                if not total_error < self.max_coverage_variation*round(self.num_landmarks/self.num_landmark_groups):
                    is_pos = False
                    break
            if is_pos:  
                self.positive.append(candidate)
                candidate.is_feasible = True
            else:
                self.negative.append(candidate)
                candidate.is_feasible = False
        log.info(f"{len(self.positive)} ({sum(self.allele_db[p.id].is_functional for p in self.positive)} functional) candidates after filter ({len(self.negative)} / {len(self.candidates)} removed)")
        return self.positive

class MappingCandidateBuilderFeasibleLandmark(BwaMappingsCandidateBuilder, LandmarkFeasibleCandidateFilter):
    def __init__(self, read_length, allele_db, num_landmarks, num_landmark_groups, expected_coverage, max_coverage_variation, minimum_coding_bases):
        BwaMappingsCandidateBuilder.__init__(self, read_length, allele_db)
        LandmarkFeasibleCandidateFilter.__init__(self, allele_db, num_landmarks, num_landmark_groups, expected_coverage, max_coverage_variation, read_length, minimum_coding_bases)
    
    def make_candidates(self, reads):
        super().make_candidates(reads)
        self.filtered_candidates = self.filter_feasible(self.candidates.values())

        positive_reads = []
        negative_reads = []
        for read in reads:
            feasible_assignments = [assign_cand for assign_cand in read.candidates if assign_cand.is_feasible]
            if feasible_assignments:
                positive_reads.append(read)
                read.candidates = feasible_assignments
            else:
                negative_reads.append(read)
        
        log.info(f"Removed {len(negative_reads)} / {len(reads)} reads that only had infeasible allele candidates")
        self.positive_reads = positive_reads
        self.negative_reads = negative_reads
        return positive_reads, self.filtered_candidates
