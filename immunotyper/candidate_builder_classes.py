from __future__ import annotations
from .common import log, create_temp_file, suppress_stdout
import abc, pysam, os
from abc import abstractmethod
from .mapper_wrappers import BwaWrapper
from collections import OrderedDict

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from .models import ModelMeta

class Landmark(object):
    def __init__(self, allele_db_landmark):
        self.pos, self.exp_cov, self.stdevs = allele_db_landmark
        self.reads = []

class Candidate(object):
    def __init__(self, name, allele_db):
        self.id = name

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
            raise AttributeError(f"Landmarks are not set for {self.id}")

    @property
    def num_copies(self) -> int:
        return int(self.model.get_value(self.copy_multipler))
    
    @property
    def model(self):
        ''''''
        try:
            return self._model
        except AttributeError:
            raise NotImplementedError("Model attribute in Candidate class needs to be set by model instance")
    
    @model.setter
    def model(self, val: ModelMeta):
        self._model = val

    def add_read(self, read_candidate):
        super().add_read(read_candidate)
        for landmark in self.landmarks:
            if read_candidate.covers_position(landmark.pos, self.id):
                landmark.reads.append(read_candidate)



class CandidateBuilder(object):
    __metaclass__ = abc.ABCMeta
    candidates = {}         # index of candidates
    no_candidates = []      # reads that do not have any candidates

    def __init__(self, read_length, allele_db, candidate_class=LandmarkCandidate):
        self.read_length = read_length
        self.allele_db = allele_db
        self.candidate_class = candidate_class
    
    def make_candidates(self, reads):
        log.info('Adding candidates to reads...')
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
                        log.info(f'Unable to make candidate {allele} (length={len(self.allele_db[allele])})\n{str(e)}')
                        continue
                candidate.add_read(read)
                
                try:
                    read.candidates.append(candidate)
                except AttributeError:
                    read.candidates = [candidate]

        log.info(f"Made {len(self.candidates)} candidates with {len([rc for c in self.candidates.values() for rc in c.reads])} read-allele pairs")

        return list(self.candidates.values())

    @abstractmethod
    def make_read_assignments(self, read):
        '''Takes a read object, return list of ReadAssignment objects for that read'''
        raise NotImplementedError

class BwaMappingsCandidateBuilder(CandidateBuilder):

    def make_read_assignments(self, read):
        return [read.reference_id_parser(m.reference_name) for m in read.mappings]




# Classes for EM Model

class NovelVariant(object):
    """Small class for each novel variant
    """
    def __init__(self, candidate: Candidate, pos: int, value: str, ref_value: str) -> None:
        #TODO doctring

        self.candidate = candidate
        self.pos = pos
        self.value = value
        self.ref_value = ref_value
        
        self.reads = []
        self.read_positions = {} # {read_id (str): position in read.seq of variant (int)}
        
    @property
    def read_support(self):
        return len(self.reads)

    def add_read(self, read, variant_position):
        self.reads.append(read)
        self.read_positions[read.id] = variant_position
        
        # add variant to read object for reverse lookup
        try:
            read.variants.append(self)
        except AttributeError:
            read.variants = [self]
            

class NovelVariantCandidate(Candidate):
    """Uses mapper, pysam.pileup to call putative novel variants according to parameters set in __init__
    saving them as attributes
    """

    def __init__(self, name, allele_db, mapper=BwaWrapper(output_sorted_bam=True)):
        super().__init__(name, allele_db)
        
        # set novel-variant parameters
        self.mapper = mapper
        self.seq = allele_db[self.id].seq
    
    

            
    def call_variants(self, allele_reference_dir) -> None:
        #TODO docstring

        if not self.reads: raise TypeError(f"No reads assigned to {self.id}")
        self.variants_dict = OrderedDict()   # {pos (int): [NovelVariant instances]}
        self.variants = []

        # build mapping output tempfiles and map with mapper
        bam_file = create_temp_file(delete=True)
        bam_file_index = create_temp_file(delete=True, suffix='.bai', prefix=bam_file.name)
        with suppress_stdout():
            try:
                mapping_file = self.mapper.map(self.reads, os.path.join(allele_reference_dir, self.id.replace('/', '#')+'.fa'),
                                        output_path=bam_file.name)
            except ValueError: # no reads map
                log.error(f"Invalid BAM mapping for allele candidate {self.id}")
                return
        pysam.index(bam_file.name, bam_file_index.name)


        # Iterate through pileup of bam
        for pileupcolumn in mapping_file.pileup():
            
            # sanity check
            if self.reads[0].reference_id_parser(pileupcolumn.reference_name) != self.id:
                raise ValueError(f"Novel variant mapping for {self.id} has pileup to {pileupcolumn.reference_name}")

            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    ref_base = self.seq[pileupcolumn.reference_pos]
                    read_base = pileupread.alignment.query_sequence[pileupread.query_position]
                    if read_base != ref_base:
                        try:
                            self.variants_dict[pileupcolumn.reference_pos][read_base].add_read(self.reads_dict[pileupread.alignment.query_name], pileupread.query_position)
                        except KeyError:
                            v = NovelVariant(self, pileupcolumn.reference_pos, read_base, ref_base)
                            v.add_read(self.reads_dict[pileupread.alignment.query_name], pileupread.query_position)
                            self.variants.append(v)
                            try:
                                self.variants_dict[pileupcolumn.reference_pos][read_base] = v
                            except KeyError:
                                self.variants_dict[pileupcolumn.reference_pos] = {read_base: v}


class EMCandidateBuilder(BwaMappingsCandidateBuilder):
    def __init__(self, read_length, allele_db, allele_reference_dir, candidate_class=NovelVariantCandidate):
        #TODO docstring
        super().__init__(read_length, allele_db, candidate_class)
        self.allele_reference_dir = allele_reference_dir

    def make_candidates(self, reads):
        candidates = super().make_candidates(reads)
        
        log.info("Calling putative novel variants for candidates...")
        for c in candidates:
            c.call_variants(self.allele_reference_dir)

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
