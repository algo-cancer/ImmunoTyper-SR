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

try:
    from builtins import object
except ImportError:
    from __builtin__ import object

import importlib

from common import log, SeqRecord, get_columns
import gurobipy
from gurobipy import LinExpr, GRB
from collections import Counter

def model(name, solver):
	def test_gurobi(name):
		try:
			model = Gurobi(name)
			log.warn('Using Gurobi')
		except ImportError:
			log.warn('Gurobi not found. Please install Gurobi and gurobipy Python package.')
			model = None
		return model

	# CPLEX does not support our quadratic models
	# def test_cplex(name):
	# 	try:
	# 		model = CPLEX(name)
	# 		log.warn('Using CPLEX')
	# 	except:
	# 		log.warn('CPLEX not found. Please install CPLEX and the following Python packages: cplex and docplex.')
	# 		model = None
	# 	return model

	def test_scip(name):
		try:
			model = SCIP(name)
			log.warn('Using SCIP')
		except ImportError:
			log.warn('SCIP not found. Please install SCIP and pyscipopt Python package.')
			model = None
		return model

	if solver == 'any':
		model = test_gurobi(name)
		# if model is None:
		# 	model = test_cplex(name)
		if model is None:
			model = test_scip(name)
		if model is None:
			raise Exception('No IP solver found. Aldy cannot solve any problems without matching IP solver. Please try installing Gurobi or SCIP.')
		return model
	else:
		fname = 'test_' + solver
		if fname in locals():
			return locals()[fname](name)
		else:
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



class Gurobi(object):
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

	def varName(self, var):
		return var.varName

	def abssum(self, vars, coeffs=None):
		total = 0
		for i, v in enumerate(vars):
			coeff = 1 if coeffs is None else coeffs[i]
			absvar = self.addVar()
			total += absvar * coeff
			self.addConstr(absvar + v >= 0)
			self.addConstr(absvar - v >= 0)
		return total

	def solve(self, objective=None, method='min', threads=6, timeLimit=4000, log_path=None): # ret obj val
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
		print 'RUNNING MODEL'
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

	def getValue(self, var):
		return var.x

	def changeUb(self, var, ub):
		var.ub = ub
		self.model.update()


class SCIP(Gurobi):
	def __init__(self, name):
		self.pyscipopt = importlib.import_module('pyscipopt')
		self.INF = 1e20
		self.model = self.pyscipopt.Model(name)

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


# CPLEX does not support our quadratic models.

# class CPLEX:
# 	def __init__(self, name):
# 		self.docplex = importlib.import_module('docplex.mp.model')
# 		self.model = self.docplex.Model(name=name)
# 		self.INF = self.model.infinity
# 		print self.INF

# 	def addConstr(self, *args, **kwargs):
# 		return self.model.add_constraint(*args, **kwargs)

# 	def addVar(self, *args, **kwargs):
# 		vtype = 'c'
# 		if 'vtype' in kwargs:
# 			vtype = 'b' if kwargs['vtype'] == 'B' else 'c'
# 			del kwargs['vtype']
# 		if vtype == 'b':
# 			v = self.model.binary_var(*args, **kwargs)
# 		else:
# 			v = self.model.continuous_var(*args, **kwargs)
# 		return v

# 	def quicksum(self, expr):
# 		return self.model.sum(expr)

# 	def update(self):
# 		pass

# 	def varName(self, var):
# 		return var.varName

# 	def abssum(self, vars, coeffs=None):
# 		return self.quicksum([self.model.abs(v) for v in vars])
		
# 	def solve(self, objective=None, method='min'): # ret obj val
# 		if objective is not None:
# 			self.objective = objective
# 		if method == 'min':
# 			self.model.minimize(self.objective)
# 		else:
# 			self.model.maximize(self.objective)
		
# 		# self.model.params.timeLimit = 60
# 		# self.model.params.threads = 2
# 		# self.model.params.outputFlag = 0
# 		# self.model.params.logFile = ''
		
# 		self.model.export_as_lp(path='./MODEL_DBG.LP')
# 		params = {}
# 		try:
# 			self.solution = self.model.solve(cplex_parameters=params, agent='local', log_output=False)
# 		except:
# 			self.model.export_as_lp(path='./MODEL_DBG.LP')
# 			self.solution = None
# 			exit(1)
# 		if not self.solution:
# 			raise NoSolutionsError(self.solution)
# 		return 'optimal', self.solution.get_objective_value()

# 	def solveAll(self, objective, var, method='min'):
# 		status, opt_value = self.solve(objective, method)
# 		solutions = [frozenset(a for a, y in var.iteritems() if round(self.getValue(y)) > 0)]
# 		solutions += get_all_solutions(
# 			self, var, opt_value,
# 			candidates=solutions[0], iteration=0, mem=list()
# 		)
# 		solutions = map(lambda x: tuple(sorted(y for y in x)), solutions)

# 		solutions = list(set(solutions))
# 		return status, opt_value, solutions

# 	def getValue(self, var):
# 		return self.solution.get_value(var)

# 	def changeUb(self, var, ub):
# 		var.ub = ub


class VarWrapper(object):
	count = None
	p_var = None
	total = None
	def __init__(self, var):
		self.var = var
	
	def X(self):
		return self.var.X



class BreakerModel(Gurobi):

	def __init__(self, expnum=None, minsize=None, expcov=None, maxcopy=None, min_cov=None, cov_weight=1000):
		self.expnum=expnum
		self.minsize=minsize
		self.expcov=expcov
		self.maxcopy=maxcopy
		self.cov_weight = cov_weight
		super(BreakerModel, self).__init__('Breaker')

	def __str__(self, expnum=None):
		result = ['BreakerModel var cov and var count']
		result.append(' '.join(['minnum', 'maxnum' 'minsize', 'expcov', 'maxcopy', 'min_cov', 'cov_weight']))
		result.append(' '.join([str(x) for x in [self.minnum, self.maxnum, self.minsize, self.expcov, self.maxcopy, self.cov_weight]]))
		return '\n'.join(result)

	def build(self, cluster, minnum=None, maxnum=None, minsize=None, expcov=None, maxcopy=None, min_cov=None, cov_weight=1000, log_path=None):
		log.debug(str(self))
		if not minnum:
			minnum=self.minnum
		if not maxnum:
			maxnum=self.minnum
		if not minsize:
			minsize=self.minsize
		if not expcov:
			expcov=self.expcov
		if not maxcopy:
			maxcopy=self.maxcopy
		if not cov_weight:
			cov_weight=self.cov_weight

		candidates = {}
		for c in cluster.candidates:
			for copy_num in range(maxcopy):
				copy = SeqRecord('{}_{}'.format(c.id, copy_num), '')
				copy.variants = c.variants
				candidates[copy.id] = copy

		## build structures to access data
		all_variants = set()
		
		for c in candidates.values():
			all_variants = all_variants.union(set(c.variants))
		candidate_variants = all_variants

		for r in cluster:
			all_variants = all_variants.union(set(r.variants))
		var_to_reads = dict([(v, {}) for v in all_variants])
		for r in cluster:
			for v in r.variants:
				try:
					var_to_reads[v][r.id] = r
				except KeyError as e:
					var_to_reads[v] = {r.id:r}

		## Build delta exists values
		copy_number = {}
		for c in candidates.values():
			delta_var = self.addVar(vtype=gurobipy.GRB.BINARY, name='D|{}|{}'.format(r.id, c.id))
			copy_number[c.id] = delta_var
			c.delta = delta_var

		## Constrain the total number of genes
		self.addConstr(self.quicksum([c.delta for c in candidates.values()]) >= minnum)
		self.addConstr(self.quicksum([c.delta for c in candidates.values()]) <= maxnum)

		## Build D[cand][read] assignment variable
		assignment = {}
		for r in cluster:
			assignment[r.id] = {}
			for c in candidates.values():
				d_var = self.addVar(vtype=gurobipy.GRB.BINARY, name='D|{}|{}'.format(r.id, c.id))
				if (r.candidate) == c.id:	## start with assignment of r to copy 0 of its prior candidate
					d_var.Start = 1
				else:
					d_var.Start = 0
				assignment[r.id][c.id] = d_var
				try:
					c.assignment[r.id] = d_var
				except AttributeError:
					c.assignment = {r.id: d_var}

				self.addConstr(c.delta >= d_var)

			r.assignment = assignment[r.id]
			self.addConstr(self.quicksum(r.assignment.values()), gurobipy.GRB.EQUAL, 1)

		## Constrain minimum size
		if minsize:
			for c in candidates.values():
				self.addConstr(self.quicksum(c.assignment.values())*c.delta >= minsize*c.delta)


		## build overall error 
		var_cov_save = {}
		var_cov_missing = {}
		var_cov_error = LinExpr()

		var_count_save = {}
		var_count_error= LinExpr()
		counter = 1
		for var, reads in var_to_reads.iteritems():
			if counter % 100 == 0: 
				log.debug('Finished adding {} of {} variants to model'.format(counter, len(var_to_reads)))
			counter += 1
			for c in candidates.values():
				
				sum_reads = LinExpr()
				for d in [r.assignment[c.id] for r in reads.values()]:
						sum_reads.add(d)

				## Variant coverage c, var building
				if var in c.variants:
					var_cov_expr = LinExpr()
					var_cov_abs  = self.addVar(vtype=gurobipy.GRB.INTEGER, name='abs_var_cov|{}|{}'.format(c.id, v))
					var_cov_data = VarWrapper(var_cov_abs)

					## save var_cov_abs
					try:
						var_cov_save[c.id][var] = var_cov_data
					except KeyError:
						var_cov_save[c.id] = {var: var_cov_data}
						c.var_cov = var_cov_save[c.id]
					
					var_cov_data.count = sum_reads

					var_cov_expr.add(sum_reads)
					var_cov_expr.add(-(c.delta*expcov))

					var_cov = self.addVar(vtype=gurobipy.GRB.INTEGER, lb=-GRB.INFINITY, name='var_cov|{}|{}'.format(c.id, v))
					self.addConstr((var_cov == var_cov_expr))
					self.addGenConstrAbs(var_cov_abs, var_cov)

					## add to error term 1:
					# var_cov_error.add(var_cov_abs/float(len(c.variants)))
					var_cov_error.add(var_cov_abs)

				elif var in candidate_variants.difference(c.variants):
					try:
						var_cov_missing[c.id][var] = sum_reads
					except KeyError:
						var_cov_missing[c.id] = {var: sum_reads}
						c.var_cov_missing = var_cov_missing[c.id]
					var_cov_error.add(sum_reads)


				## variant count error
				var_count_abs = self.addVar(vtype=gurobipy.GRB.INTEGER, name='abs_var_count|{}|{}'.format(c.id, v))
				var_count_data = VarWrapper(var_count_abs)
				var_count_expr = LinExpr()
				try:
					var_count_save[c.id][var] = var_count_data
				except KeyError:
					var_count_save[c.id] = {var: var_count_data}
					c.var_count = var_count_save[c.id]

				var_count_data.count = sum_reads

				var_count_expr.add(sum_reads)
				var_count_expr.add(-(c.delta*expcov))
				

				var_count = self.addVar(vtype=gurobipy.GRB.INTEGER, lb=-gurobipy.GRB.INFINITY, name='var_count|{}|{}'.format(c.id, v))
				self.addConstr((var_count == var_count_expr))
				self.addGenConstrAbs(var_count_abs, var_count)

				## make P variable
				p=None
				if min_cov:
					p = self.add_p_var(var_count_abs, sum_reads, min_cov)
					var_count_final = self.addVar(vtype=gurobipy.GRB.INTEGER)
					self.addConstr(var_count_final == var_count_abs*p)
				else:
					var_count_final = var_count_abs
				var_count_data.total = var_count_final


				var_count_expr = LinExpr()
				var_count_data.p_var = p

				## add to error term 2
				var_count_error.add(var_count_final)

		var_cov_error = var_cov_error *cov_weight
		self.objective = LinExpr(var_cov_error + var_count_error)

		self.var_count_error = var_count_error
		self.var_cov_error = var_cov_error

		self.candidates = candidates
		self.reads = cluster.reads
		self.var_to_reads = var_to_reads


	def get_assignments(self):
		result = []
		for c in self.candidates.values():
			if c.delta.X > 0:
				reads = []
				for r, v in c.assignment.iteritems():
					if v.X > 0:
						reads.append(r)
				result.append((c.id.split('_')[0], reads))
		return result, self.candidates



	def add_p_var(self, var_count_abs, count, min_cov):
		p = self.addVar(vtype=gurobipy.GRB.BINARY)
		p.constraint = self.addConstr(count/float(min_cov) - 0.01 > p)
		return p		


def print_model_assignments(m, print_all_var=False):
	print(print_assignments(m.candidates, m, print_all_var=print_all_var))

def print_assignments(candidates, model, print_all_var=False):
	result = ['\nCoding error: {}\nNon-coding error= {}\n'.format(model.var_cov_error.getValue(), model.var_count_error.getValue())]
	calls = []
	for c in candidates.values():
		if c.delta.X > 0:
			calls.append(c.id)
			result.append('----\n'+c.id+'\n----')
			result.append(str(c.delta.x))
			reads = []
			for r ,v in c.assignment.iteritems():
				if v.X > 0: 
					reads.append(r)
			result.append(get_columns(Counter(['_'.join(r.split('_')[:2]) for r in reads]).most_common()))
			var_cov_data = [('Variant', 'abs err', 'count reads')]
			try:
				for v, var_cov_abs in c.var_cov.iteritems():
					var_cov_data.append((v, var_cov_abs.X(), var_cov_abs.count.getValue()))
				result.append(get_columns(var_cov_data))
			except AttributeError:
				print('{} doesnt have variant coverage data?'.format(c.id))
			result.append('\nVariants not present in allele')
			var_missing_data = [('Variant', 'Count')]
			try:
				for v, val in c.var_cov_missing.iteritems():
					var_missing_data.append((v, val.getValue()))
				result.append(get_columns(var_missing_data))
			except AttributeError:
				print('{} doesnt have variant missing data?'.format(c.id))
			var_count_data = []
			
			try:
				for v, var_count_abs in c.var_count.iteritems():
					if var_count_abs.count.getValue() > 0:
						var_count_data.append((v, 
												var_count_abs.total.X, 
												var_count_abs.X(),
												var_count_abs.count.getValue(),
												'None' if not var_count_abs.p_var else var_count_abs.p_var.X,
												'None' if not var_count_abs.p_var else var_count_abs.p_var.constraint.RHS)
											)
				result.append('There are {} variants with an average coverage of {}\n'.format(len(var_count_data), round(sum([x[3] for x in var_count_data])/float(len([x for x in var_count_data if x[3] > 0])))))
			except AttributeError:
				print('{} doesnt have variant count data?'.format(c.id))
	
			if print_all_var:
				result.append(get_columns([('Variant', 'error', 'abs error' 'count', 'p var', 'p_constRHS')]+var_count_data))
	result = ['{} allele calls:'.format(len(calls))] + calls + result
	return '\n'.join(result)