from __future__ import annotations
from abc import ABCMeta, abstractmethod
from typing import List
from .solvers import SolverMeta, GurobiSolver
from .common import Read, log

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from .candidate_builder_classes import Candidate


class ModelMeta(metaclass=ABCMeta):

    @property
    def candidates(self) -> List[Candidate]:
        """Getter for candidates set by self.build

        Returns:
            list: list of self.candidate_class instances
        """

        return self._candidates


    def build(self, filtered_reads: List[Read], candidates: List[Candidate]) -> None:
        """Build model

        Args:
            filtered_reads ([common.Read]): List of common.Read objects. Include candidate attributes.
            candidates ([candidate_builder_class.Candidate]): List of candidate_builder_class.Candidate objects
        """
        for c in candidates:
            c.model = self
        self._candidates = candidates
        self.reads = filtered_reads

        return self._build(filtered_reads, candidates)
    
    @abstractmethod
    def _build(self, filtered_reads: List[Read], candidates: List[Candidate]) -> None:
        """Private function to implement in child classes to build actual model
        """
        raise NotImplementedError


    def get_allele_calls(self, no_cnv=False, functional_only=False) -> list:
        """Returns list of called alleles, where the number of copies equals

        Args:
            no_cnv (bool, optional): If True, alleles are unique in output, else number of instances of allele = number of copies. Defaults to False.
            functional_only (bool, optional): If true, only return functional alleless

        Returns:
            list: List of allele id strings
        """
        called = []
        for c in self.candidates:
            if functional_only and not c.is_functional:
                continue
            for copy in range(c.num_copies):
                called.append(c.id)
        called = sorted(called)

        return called if not no_cnv else list(set(called))

    def write_allele_calls(self, output_path: str, functional_only=False) -> None:
        """Writes allele calls to file, with identical alleles tab-separated on same line

        Args:
            output_path (str): Path to output file
            functional_only (bool, optional): If true, only output functional alleles. Defaults to False.
        """
        called_alleles = self.get_allele_calls(functional_only=functional_only)
        
        with open(output_path, 'w') as f:
            for allele in sorted(called_alleles):
                line_alleles = [allele]
                if hasattr(self, 'allele_db'):
                    identical_alleles = self.allele_db.get_identical_alleles(allele)
                    if identical_alleles:
                        line_alleles.extend(identical_alleles)
                f.write('\t'.join(line_alleles) + '\n')

class ShortReadModelTotalErrorDiscardObj(ModelMeta):

    def __init__(self, solver: SolverMeta = GurobiSolver, 
                 stdev_coefficient=1.5, 
                 min_landmark_depth=.3, 
                 num_landmarks=6, 
                 num_landmark_groups=6, 
                 maxcopy=4, 
                 sequencing_error_rate=0.02, 
                 discard_penalty_multiplier=2):
        
        super().__init__() # init ModelMeta

        # set solver as another parent class to access functions directly
        self.__class__ = type(self.__class__.__name__,
                              (solver, ModelMeta),
                              dict(self.__class__.__dict__))
        
        super(self.__class__, self).__init__('ShortReadModelTotalErrorDiscardObj') #init solver

        # set attributes
        self.maxcopy=maxcopy
        self.stdev_coefficient = stdev_coefficient
        self.min_landmark_depth = min_landmark_depth
        self.num_landmark_groups = num_landmark_groups
        self.sequencing_error_rate = sequencing_error_rate
        self.discard_penalty_multiplier = discard_penalty_multiplier
        self.num_landmarks = num_landmarks


    def _build(self, reads, candidates):
        objective = self.LinExpr()		    # objective function
        objective_array = []                # for holding variables in the objective function

        assignment = {}				    # {read.id: {c.id: D[cand][read] variable}}
        discard_error = self.LinExpr()
        discard_error_array = []
        already_seen_reads = set()
        self.depth_errors = {}
        self._candidates = candidates

        log.debug(f'Building model for {len(reads)} reads and {len(candidates)} candidates')

        num_candidates_processed = 0
        for allele_candidate in candidates:
            allele_candidate.copies = [self.add_var(type='b', name=f'copies|{x}|{allele_candidate.id}') for x in range(self.maxcopy+1)]
            allele_candidate.copy_multipler = self.quicksum([i*copy_var for i, copy_var in enumerate(allele_candidate.copies)])
            
            #constraint (1) - One-hot encoding of all possible copy numbers for each allele
            self.add_constr(self.quicksum(allele_candidate.copies) == 1, name=f"constraint_1-{allele_candidate.id}")

            allele_candidate.assignment = {}

            for read in allele_candidate.reads: # all reads that have allele_candidate as a candidate
                if read.id not in already_seen_reads:
                    already_seen_reads.add(read.id)
                    
                    # make discard variables
                    read.discard_var = self.add_var(type='b', name='discard|{}'.format(read.id))
                    read.discard_penalty = round(read.primary_mapping.query_alignment_length*self.sequencing_error_rate*self.discard_penalty_multiplier)
                    discard_error = discard_error + read.discard_penalty*read.discard_var
                    discard_error_array.append(read.discard_penalty*read.discard_var)

                    # reset attribute values
                    assignment[read.id] = {}
                    read.assignment = {}
                    read.assignment_expr = self.LinExpr()

                # add read-candidate assignment variables (d_var)
                d_var = self.add_var(type='b', name='D|{}|{}'.format(read.id, allele_candidate.id))
                try:
                    assignment[read.id][allele_candidate.id] = d_var
                except KeyError:
                    assignment[read.id] = {allele_candidate.id: d_var}
                allele_candidate.assignment[read.id] = d_var
                read.assignment[allele_candidate.id] = d_var
                objective  = objective + read.get_distance(allele_candidate)*d_var
                objective_array.append(read.get_distance(allele_candidate)*d_var)
            
            # constraint (3): enforce assigned reads -> called allele
            self.add_constr(self.quicksum(allele_candidate.copies[1:])*len(reads) >= self.quicksum(allele_candidate.assignment.values()), name=f"constr_4-{allele_candidate.id}")

            # calculate depth error for allele_candidate
            allele_candidate.read_coverage = {}
            for i, group in enumerate(allele_candidate.landmark_groups):
                
                depth_error = self.LinExpr()
                depth_error_array = []          # for linear expression suming expected/assigned depth differences for landmarks in group  
                
                lambdas = []
                
                for j, landmark in enumerate(group):

                    depth_diff = self.LinExpr()
                    depth_diff_array = []               # for linear expression summing expected/assigned depth differences for this landmark

                    # add assigned read coverage for landmark to depth_diff_array
                    reads_covering_pos = [r.assignment[allele_candidate.id] for r in landmark.reads]		# get d_var for reads covering pos		
                    allele_candidate.read_coverage[landmark.pos] = reads_covering_pos                       # save as attr for post analysis                    
                    for r in reads_covering_pos:
                        depth_diff = depth_diff + r
                        depth_diff_array.append(r)
                    
                    # add expected depth to depth_diff_array
                    depth_diff = depth_diff + -landmark.exp_cov*allele_candidate.copy_multipler
                    depth_diff_array.append(-landmark.exp_cov*allele_candidate.copy_multipler)

                    # make pseudo-abs value of depth_diff_array
                    diff_abs = self.add_var(name=f'diff_abs_{allele_candidate.id}_group{i}_landmark{j}', lb=0, type='i')
                    self.add_constr(diff_abs + self.quicksum(depth_diff_array) >= 0, name=f'abs_val_1-{allele_candidate.id}_group{i}_landmark{j}')
                    self.add_constr(diff_abs - self.quicksum(depth_diff_array) >= 0, name=f'abs_val_2-{allele_candidate.id}_group{i}_landmark{j}')

                    depth_error = depth_error + diff_abs
                    depth_error_array.append(diff_abs)          # add this groups depth error to allele sum

                    # constraint 5
                    self.add_constr(self.quicksum(reads_covering_pos) >= self.min_landmark_depth*landmark.exp_cov*allele_candidate.copy_multipler, name=f"constr_5-{allele_candidate.id}_group{i}_landmark{j}")
                    
                    lambdas.append([self.stdev_coefficient*x for x in landmark.stdevs])
                

                
                # constrain depth error sum for this group by sum of standard deviation difference allowance across all landmarks in group
                lambdas = zip(*lambdas) # group landmark stdev allowances by copy number
                self.add_constr(self.quicksum(depth_error_array) <= (self.quicksum([c*sum(lam) for c, lam in zip(allele_candidate.copies, lambdas)])), name=f'depth_error-group_{i}-{allele_candidate.id}')

            allele_candidate.depth_error = depth_error
            self.depth_errors[allele_candidate.id] = depth_error
        
        
        for read in reads:
            # add constraint 2: prevent discarded reads from being assigned and enforce assignment to at most 1 copy
            try:
                self.add_constr(self.quicksum(read.assignment.values()) == 1-read.discard_var, name='Constr_5-{}'.format(read.id))		# Constraint 5
            except AttributeError:
                # ignore reads with no candidates
                continue

        objective_array.extend(discard_error_array)
        self.objective = self.quicksum(objective_array)
        self.discard_penalty = self.quicksum(discard_error_array)




class ShortReadModelTotalErrorDiscardObjAlleleClusters(ShortReadModelTotalErrorDiscardObj):
    
    def __init__(self, clusters, *args, **kwargs):
        self.clusters = clusters

        super().__init__(*args, **kwargs)

    def _build(self, reads, candidates):
        self.reads = reads
        objective = self.LinExpr()		    # objective function
        objective_array = []                # for holding variables in the objective function

        assignment = {}				    # {read.id: {c.id: D[cand][read] variable}}
        discard_error = self.LinExpr()
        discard_error_array = []
        already_seen_reads = set()
        self.depth_errors = {}
        self._candidates = candidates

        log.debug(f'Building model for {len(reads)} reads and {len(candidates)} candidates')

        num_candidates_processed = 0
        for allele_candidate in candidates:
            allele_candidate.copies = [self.add_var(type='b', name=f'copies|{x}|{allele_candidate.id}') for x in range(self.maxcopy+1)]
            allele_candidate.copy_multipler = self.quicksum([i*copy_var for i, copy_var in enumerate(allele_candidate.copies)])
            
            #constraint (1) - One-hot encoding of all possible copy numbers for each allele
            self.add_constr(self.quicksum(allele_candidate.copies) == 1, name=f"constraint_1-{allele_candidate.id}")

            allele_candidate.assignment = {}

            for read in allele_candidate.reads: # all reads that have allele_candidate as a candidate
                if read.id not in already_seen_reads:
                    already_seen_reads.add(read.id)
                    
                    # make discard variables
                    read.discard_var = self.add_var(type='b', name='discard|{}'.format(read.id))
                    read.discard_penalty = round(read.primary_mapping.query_alignment_length*self.sequencing_error_rate*self.discard_penalty_multiplier)
                    discard_error = discard_error + read.discard_penalty*read.discard_var
                    discard_error_array.append(read.discard_penalty*read.discard_var)

                    # reset attribute values
                    assignment[read.id] = {}
                    read.assignment = {}
                    read.assignment_expr = self.LinExpr()

                # add read-candidate assignment variables (d_var)
                d_var = self.add_var(type='b', name='D|{}|{}'.format(read.id, allele_candidate.id))
                try:
                    assignment[read.id][allele_candidate.id] = d_var
                except KeyError:
                    assignment[read.id] = {allele_candidate.id: d_var}
                allele_candidate.assignment[read.id] = d_var
                read.assignment[allele_candidate.id] = d_var
                objective  = objective + read.get_distance(allele_candidate)*d_var
                objective_array.append(read.get_distance(allele_candidate)*d_var)
            
            # constraint (3): enforce assigned reads -> called allele
            self.add_constr(self.quicksum(allele_candidate.copies[1:])*len(reads) >= self.quicksum(allele_candidate.assignment.values()), name=f"constr_4-{allele_candidate.id}")

            # calculate depth error for allele_candidate
            allele_candidate.read_coverage = {}
            for i, group in enumerate(allele_candidate.landmark_groups):
                
                depth_error = self.LinExpr()
                depth_error_array = []          # for linear expression suming expected/assigned depth differences for landmarks in group  
                
                lambdas = []
                
                for j, landmark in enumerate(group):

                    depth_diff = self.LinExpr()
                    depth_diff_array = []               # for linear expression summing expected/assigned depth differences for this landmark

                    # add assigned read coverage for landmark to depth_diff_array
                    reads_covering_pos = [r.assignment[allele_candidate.id] for r in landmark.reads]		# get d_var for reads covering pos		
                    allele_candidate.read_coverage[landmark.pos] = reads_covering_pos                       # save as attr for post analysis                    
                    for r in reads_covering_pos:
                        depth_diff = depth_diff + r
                        depth_diff_array.append(r)
                    
                    # add expected depth to depth_diff_array
                    depth_diff = depth_diff + -landmark.exp_cov*allele_candidate.copy_multipler
                    depth_diff_array.append(-landmark.exp_cov*allele_candidate.copy_multipler)

                    # make pseudo-abs value of depth_diff_array
                    diff_abs = self.add_var(name=f'diff_abs_{allele_candidate.id}_group{i}_landmark{j}', lb=0, type='i')
                    self.add_constr(diff_abs + self.quicksum(depth_diff_array) >= 0, name=f'abs_val_1-{allele_candidate.id}_group{i}_landmark{j}')
                    self.add_constr(diff_abs - self.quicksum(depth_diff_array) >= 0, name=f'abs_val_2-{allele_candidate.id}_group{i}_landmark{j}')

                    depth_error = depth_error + diff_abs
                    depth_error_array.append(diff_abs)          # add this groups depth error to allele sum

                    # constraint 5
                    self.add_constr(self.quicksum(reads_covering_pos) >= self.min_landmark_depth*landmark.exp_cov*allele_candidate.copy_multipler, name=f"constr_5-{allele_candidate.id}_group{i}_landmark{j}")
                    
                    lambdas.append([self.stdev_coefficient*x for x in landmark.stdevs])
                

                
                # constrain depth error sum for this group by sum of standard deviation difference allowance across all landmarks in group
                lambdas = zip(*lambdas) # group landmark stdev allowances by copy number
                self.add_constr(self.quicksum(depth_error_array) <= (self.quicksum([c*sum(lam) for c, lam in zip(allele_candidate.copies, lambdas)])), name=f'depth_error-group_{i}-{allele_candidate.id}')

            allele_candidate.depth_error = depth_error
            self.depth_errors[allele_candidate.id] = depth_error
        
        
        for read in reads:
            # add constraint 2: prevent discarded reads from being assigned and enforce assignment to at most 1 copy
            try:
                self.add_constr(self.quicksum(read.assignment.values()) == 1-read.discard_var, name='Constr_5-{}'.format(read.id))		# Constraint 5
            except AttributeError:
                # ignore reads with no candidates
                continue

        objective_array.extend(discard_error_array)
        self.objective = self.quicksum(objective_array)
        self.discard_penalty = self.quicksum(discard_error_array)

        candidates_dict = dict((c.id, c) for c in self._candidates)

        for cluster in self.clusters:
            cluster.cluster_total_copies = self.quicksum([candidates_dict[a.id].copy_multipler for a in cluster.alleles if a.id in candidates_dict])
            self.add_constr(cluster.cluster_total_copies <= cluster.cnv_call+1, name=f"cluser-{cluster.id}-constr")
            self.add_constr(cluster.cluster_total_copies >= cluster.cnv_call-1, name=f"cluser-{cluster.id}-constr")
