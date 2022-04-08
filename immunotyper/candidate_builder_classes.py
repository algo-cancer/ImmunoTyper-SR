from .common import log
import abc
from abc import abstractmethod


class CandidateBuilder(object):
    __metaclass__ = abc.ABCMeta
    candidates = {}         # index of candidates
    no_candidates = []      # reads that do not have any candidates

    def __init__(self, read_length, allele_db):
        self.read_length = read_length
        self.allele_db = allele_db
    
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
                        candidate = Candidate(allele, self.allele_db)
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

class Candidate():
    def __init__(self, name, allele_db):
        self.id = name

        self.reads = []

        self.is_ignored = allele_db[self.id].is_ignored
        try:
            self.landmark_groups = [[CandidateLandmark(landmark) for landmark in group] for group in allele_db[self.id].landmark_groups] #[read]
            self.landmarks = [l for g in self.landmark_groups for l in g]
        except KeyError:
            raise KeyError('Follow mapping has reference name of {} not in allele database\n{}'.format(self.id, str(self.mapping)))
        except (AttributeError, TypeError):
            raise AttributeError(f"Landmarks are not set for {self.id}")
    
    def add_read(self, read_candidate):
        self.reads.append(read_candidate)
        for landmark in self.landmarks:
            if read_candidate.covers_position(landmark.pos, self.id):
                landmark.reads.append(read_candidate)
        read_candidate.candidate = self

    @staticmethod
    def reference_id_parser(raw_id):
        return raw_id.split('|')[1]

class CandidateLandmark(object):
    def __init__(self, allele_db_landmark):
        self.pos, self.exp_cov, self.stdevs = allele_db_landmark
        self.reads = []
