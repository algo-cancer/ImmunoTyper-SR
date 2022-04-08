# 786
# LICENCE
# Copyright (c) 2017, Simon Fraser University, Massachusetts Institute of Technology,
#                     Indiana University Bloomington.
# All rights reserved.
# Redistribution and use in both source and binary forms, without modification,
# is permitted provided that the following conditions are met:
# - This software shall NOT BE USED in any commercial environment,
#   unless explicitely allowed by the authors in the writing.
# - Redistributions in both source and binary form must reproduce the above
#   copyright notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# - Neither the name of the Simon Fraser University, Indiana University Bloomington,
#   Massachusetts Institute of Technology, nor the names of its contributors may be used
#   to endorse or promote products derived from this software without specific prior written permission.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
# TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
# OF SUCH DAMAGE.
#
# Additional restrictions, as specified at http://compbio.cs.sfu.ca/nwp-content/projects/aldy, do apply.

import importlib
from collections import Counter
import gurobipy
from gurobipy import LinExpr, GRB, QuadExpr
from .common import log, SeqRecord, get_columns

def model(name, solver):
    def test_gurobi(name):
        try:
            _model = Gurobi(name)
            log.warn('Using Gurobi')
        except ImportError:
            log.warn('Gurobi not found. Please install Gurobi and gurobipy Python package.')
            _model = None
        return _model

    # CPLEX does not support our quadratic models
    # def test_cplex(name):
    #   try:
    #       _model = CPLEX(name)
    #       log.warn('Using CPLEX')
    #   except:
    #       log.warn('CPLEX not found. Please install CPLEX and the following Python packages: cplex and docplex.')
    #       _model = None
    #   return _model

    def test_scip(name):
        try:
            _model = SCIP(name)
            log.warn('Using SCIP')
        except ImportError:
            log.warn('SCIP not found. Please install SCIP and pyscipopt Python package.')
            _model = None
        return _model

    if solver == 'any':
        _model = test_gurobi(name)
        # if _model is None:
        #   _model = test_cplex(name)
        if _model is None:
            _model = test_scip(name)
        if _model is None:
            raise Exception('No IP solver found. Aldy cannot solve any problems without matching IP solver. Please try installing Gurobi or SCIP.')
        return _model
    else:
        fname = 'test_' + solver
        if fname in locals():
            return locals()[fname](name)
        raise Exception('IP solver {} is not supported'.format(solver))


def get_all_solutions(c, var, opt, candidates, iteration=0, mem=[]):
    if candidates in mem:
        return []
    solutions = []
    for a in candidates:
        c.changeUb(var[a], 0)
        try:
            status, obj = c.solve()
            if status == 'optimal' and abs(obj - opt) < 1e-6:
                new_solution = set(a for a, y in var.items() if round(c.getValue(y)) > 0)
                new_candidates = (set(candidates) - set([a])) | new_solution
                solutions.append(frozenset(new_solution))
                solutions += get_all_solutions(c, var, opt, new_candidates, iteration + 1, mem)
        except NoSolutionsError:
            pass
        c.changeUb(var[a], 1)
        mem.append(candidates)
    return solutions


class NoSolutionsError(Exception):
    pass



class Gurobi():
    def __init__(self, name):
        self.gurobipy = importlib.import_module('gurobipy')

        self.INF = self.gurobipy.GRB.INFINITY
        self.GUROBI_STATUS = {
            getattr(self.gurobipy.GRB.status, v): v
            for v in dir(self.gurobipy.GRB.status)
            if v[:2] != '__'
        }
        self.model = self.gurobipy.Model(name)
        self.model.reset()

    def addConstr(self, *args, **kwargs):
        return self.model.addConstr(*args, **kwargs)

    def addGenConstrAbs(self, *args, **kwargs):
        return self.model.addGenConstrAbs(*args, **kwargs)

    def addVar(self, *args, **kwargs):
        if 'vtype' in kwargs and kwargs['vtype'] == 'B':
            kwargs['vtype'] = self.gurobipy.GRB.BINARY
        v = self.model.addVar(*args, **kwargs)
        self.model.update()
        return v

    def quicksum(self, expr):
        return self.gurobipy.quicksum(expr)

    def update(self):
        self.model.update()

    @staticmethod
    def varName(var):
        return var.varName

    def abssum(self, variables, coeffs=None):
        total = 0
        for i, v in enumerate(variables):
            coeff = 1 if coeffs is None else coeffs[i]
            absvar = self.addVar()
            total += absvar * coeff
            self.addConstr(absvar + v >= 0)
            self.addConstr(absvar - v >= 0)
        return total

    def solve(self, objective=None, method='min', threads=6, timeLimit=4000, log_path=None, MIPGap=None): # ret obj val
        if objective is not None:
            self.objective = objective
        self.model.setObjective(
            self.objective,
            self.gurobipy.GRB.MINIMIZE if method == 'min' else self.gurobipy.GRB.MAXIMIZE
        )
        self.model.update()

        if timeLimit: self.model.params.timeLimit = timeLimit
        self.model.params.threads = threads
        # self.model.params.outputFlag = 1
        self.model.params.LogToConsole = 1
        if not log_path:
            log_path = 'gurobi.log'
        self.model.params.logFile = log_path
        if MIPGap: self.model.params.MIPGap = MIPGap
        print('RUNNING MODEL')
        self.model.optimize()

        status = self.GUROBI_STATUS[self.model.status]
        if self.model.status == self.gurobipy.GRB.INFEASIBLE:
            raise NoSolutionsError(status)
        return status.lower(), self.model.objVal

    def solveAll(self, objective, var, method='min'):
        status, opt_value = self.solve(objective, method)
        solutions = [frozenset(a for a, y in var.items() if round(self.getValue(y)) > 0)]
        solutions += get_all_solutions(
            self, var, opt_value,
            candidates=solutions[0], iteration=0, mem=list()
        )
        solutions = [tuple(sorted(y for y in x)) for x in solutions]

        solutions = list(set(solutions))
        return status, opt_value, solutions

    @staticmethod
    def getValue(var):
        return var.x

    def changeUb(self, var, ub):
        var.ub = ub
        self.model.update()


class SCIP(Gurobi):
    def __init__(self, name):
        self.pyscipopt = importlib.import_module('pyscipopt')
        self.INF = 1e20
        self.model = self.pyscipopt.Model(name)
        super().__init__('SCIP')

    def update(self):
        pass

    def addConstr(self, *args, **kwargs):
        return self.model.addCons(*args, **kwargs)

    def addVar(self, *args, **kwargs):
        return self.model.addVar(*args, **kwargs)

    def quicksum(self, expr):
        return self.pyscipopt.quicksum(expr)

    def solve(self, objective=None, method='min'):
        if objective is not None:
            self.objective = objective
        self.model.setObjective(
            self.objective,
            'minimize' if method == 'min' else 'maximize'
        )
        # self.model.params.timeLimit = 60
        self.model.setRealParam('limits/time', 120)
        self.model.hideOutput()
        self.model.optimize()

        status = self.model.getStatus()
        if status == 'infeasible':
            raise NoSolutionsError(status)
        return status, self.model.getObjVal()

    def varName(self, var):
        return var.name

    def getValue(self, var):
        return self.model.getVal(var)

    def changeUb(self, var, ub):
        self.model.freeTransform()
        self.model.chgVarUb(var, ub)


class ShortReadModelTotalErrorDiscardObj(Gurobi):
    def __init__(self, stdev_coefficient=1, min_landmark_depth=.3, num_landmarks=5, num_landmark_groups=2, maxcopy=4, sequencing_error_rate=0.02, discard_penalty_multiplier=2):
        self.maxcopy=maxcopy
        self.stdev_coefficient = stdev_coefficient
        self.min_landmark_depth = min_landmark_depth
        self.num_landmark_groups = num_landmark_groups
        self.sequencing_error_rate = sequencing_error_rate
        self.discard_penalty_multiplier = discard_penalty_multiplier
        self.num_landmarks = num_landmarks
        super().__init__('ShortReadModelTotalErrorDiscardObjNew')


    def build(self, reads, candidates):
        objective = LinExpr()		    # objective function

        assignment = {}				    # {read.id: {c.id: D[cand][read] variable}}
        discard_error = LinExpr()
        already_seen_reads = set()
        self.depth_errors = {}
        self.candidates = candidates
        self.reads = reads

        log.info(f'Building model for {len(reads)} reads and {len(candidates)} candidates')

        num_candidates_processed = 0
        for allele_candidate in candidates:
            allele_candidate.copies = [self.addVar(vtype=gurobipy.GRB.BINARY, name=f'copies|{x}') for x in range(self.maxcopy+1)]
            allele_candidate.copy_multipler = LinExpr(self.quicksum([i*copy_var for i, copy_var in enumerate(allele_candidate.copies)]))
            #constraint (1) - One-hot encoding of all possible copy numbers for each allele
            self.addConstr(self.quicksum(allele_candidate.copies), GRB.EQUAL, 1)

            allele_candidate.assignment = {}

            for read in allele_candidate.reads: # all reads that have allele_candidate as a candidate
                if read.id not in already_seen_reads:
                    # make discard variables
                    read.discard_var = self.addVar(vtype=gurobipy.GRB.BINARY, name='discard|{}'.format(read.id))
                    read.discard_penalty = round(read.primary_mapping.query_alignment_length*self.sequencing_error_rate*self.discard_penalty_multiplier)
                    discard_error.add(read.discard_var, read.discard_penalty)
                    # reset attribute values
                    assignment[read.id] = {}
                    read.assignment = {}
                    read.assignment_expr = LinExpr()
                    already_seen_reads.add(read.id)

                # add read-candidate assignment variables (d_var)
                d_var = self.addVar(vtype=gurobipy.GRB.BINARY, name='D|{}|{}'.format(read.id, allele_candidate.id))
                try:
                    assignment[read.id][allele_candidate.id] = d_var
                except KeyError:
                    assignment[read.id] = {allele_candidate.id: d_var}
                allele_candidate.assignment[read.id] = d_var
                read.assignment[allele_candidate.id] = d_var
                objective.add(d_var, read.get_distance(allele_candidate))

            
            # # constraint (3): enforce called allele -> assigned reads
            # self.addConstr(self.quicksum(allele_candidate.assignment.values()), gurobipy.GRB.GREATER_EQUAL, allele_candidate.copies, name='Constr_3-{}'.format(allele_candidate.id))	# constraint 3

            # constraint (4): enforce assigned reads -> called allele
            self.addConstr(self.quicksum(allele_candidate.copies[1:])*len(reads), GRB.GREATER_EQUAL, self.quicksum(allele_candidate.assignment.values()))

            # calculate depth error
            allele_candidate.read_coverage = {}
            for i, group in enumerate(allele_candidate.landmark_groups):
                depth_error = LinExpr()
                lambdas = []
                for landmark in group:
                    reads_covering_pos = [r.assignment[allele_candidate.id] for r in landmark.reads]		# get d_var for reads covering pos		
                    allele_candidate.read_coverage[landmark.pos] = reads_covering_pos
                    
                    depth_diff = LinExpr()

                    for r in reads_covering_pos:
                        depth_diff.add(r)
                    
                    # abs value
                    depth_diff.add(allele_candidate.copy_multipler, -landmark.exp_cov)

                    diff_abs = self.addVar(lb=0, vtype=GRB.INTEGER)
                    self.addConstr(diff_abs + depth_diff >= 0, name=f'abs_val_1-{allele_candidate.id}')
                    self.addConstr(diff_abs - depth_diff >= 0, name=f'abs_val_2-{allele_candidate.id}')

                    depth_error.add(diff_abs)

                    self.addConstr(self.quicksum(reads_covering_pos) >= self.min_landmark_depth*landmark.exp_cov*allele_candidate.copy_multipler)
                    lambdas.append([self.stdev_coefficient*x for x in landmark.stdevs])
                
                lambdas = zip(*lambdas)

                self.addConstr(depth_error, gurobipy.GRB.LESS_EQUAL, (self.quicksum([c*sum(lam) for c, lam in zip(allele_candidate.copies, lambdas)])), name=f'depth_error-group_{i}-{allele_candidate.id}')

            allele_candidate.depth_error = depth_error
            self.depth_errors[allele_candidate.id] = depth_error
        
        for read in reads:
            # add constraint to prevent discarded reads from being assigned and enforce assignment to at most 1 copy
            try:
                self.addConstr(self.quicksum(read.assignment.values()), gurobipy.GRB.EQUAL, 1-read.discard_var, name='Constr_5-{}'.format(read.id))		# Constraint 5
            except AttributeError:
                # ignore reads with no candidates
                # log.warn(f'read {read.id} has an AttributeError')
                # print(read.id)
                continue

        
        self.objective = objective+discard_error
        self.discard_penalty = discard_error
