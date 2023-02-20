import importlib
from abc import ABCMeta, abstractmethod
from typing import Literal

class NoSolutionsError(Exception):
    """Raised if a model is infeasible."""
    pass


class SolverMeta(metaclass=ABCMeta):
    @abstractmethod
    def solve(time_limit, threads, log_path) -> None:
        """_summary_

        Args:
            time_limit (int): Max time limit in seconds
            threads (int): Max thread usage
            log_path (str): Path to log file to write
        """
        raise NotImplementedError



class IlpSolverMeta(SolverMeta):
    """Meta class for ILP solvers. Child classes must implement all abstract classes"""
    SOLUTION_PRECISION = 1e-5
    VALID_VARIABLE_TYPES = set(['b', 'c', 'i'])
    

    def add_var(self, name: str, type: Literal['b', 'c', 'i'], lb: int=0, ub=None):
        """Adds variable to solver model. Calls self._add_var

        Args:
            name (str): Name of variable
            type (str): One of ('c', 'b', 'i')
            lb (int, optional): Lower bound of variable. Defaults to 0.
            ub (_type_, optional): Upper bound of variable. Defaults to None which will result in being set to infinity

        Raises:
            ValueError: If type is invalid

        Returns:
            None
        """
        type = type.lower()
        if type not in self.VALID_VARIABLE_TYPES:
            raise ValueError(f'Variable type {type} not one of valid options: {str(self.VALID_VARIABLE_TYPES)}')
        if ub is None:
            ub = self.INF
        return self._add_var(name, type, lb, ub)


    @abstractmethod
    def _add_var(self, name, type, lb, ub):
        '''Private method to actually create the variant - to be implemented in child class'''
        raise NotImplementedError


    @abstractmethod
    def solve(self, method: Literal['minimize', 'maximize'], 
              time_limit: int,
               threads: int, 
               log_path: str) -> None:
        """Solve the model using self.objective as the objective function

        Args:
            method (Literal[&#39;minimize&#39;, &#39;maximize&#39;]): Maximize or minimize the objective functions
            time_limit (int): Time limit in seconds
            threads (int): 
            log_path (str): Path to save solver log if applicable
        """


    @abstractmethod
    def add_constr(self, name: str, *args, **kwargs):
        """Add constraint. Args should include linear expression of variables.

        Args:
            name (str): Name of variable
        """
        raise NotImplementedError


    @abstractmethod
    def quicksum(self, expr: list):
        """Sum list of variables and coefficients and return linear expression.

        Args:
            expr (list): List of variables and coefficiants in form [var, -coef*var2,...]
        """
        raise NotImplementedError


    @staticmethod
    @abstractmethod
    def var_name(var):
        """Return variable name

        Args:
            var: Variable instance
        """
        raise NotImplementedError


    @abstractmethod
    def LinExpr(self, *args, **kwargs):
        """Return a linear expression instance
        """
        raise NotImplementedError


    @abstractmethod
    def get_value(self, var):
        """Return variable value

        Args:
            var: variable instance
        """
        raise NotImplementedError
    





class GurobiSolver(IlpSolverMeta):
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

    def add_constr(self, *args, **kwargs):
        return self.model.addConstr(*args, **kwargs)

    def addGenConstrAbs(self, *args, **kwargs):
        return self.model.addGenConstrAbs(*args, **kwargs)

    def _add_var(self, name: str, type: str, lb: int=0, ub=None, *args, **kwargs):    
        if type == "b":
            type = self.gurobipy.GRB.BINARY
        elif type == "i":
            type = self.gurobipy.GRB.INTEGER
        else:
            type = self.gurobipy.GRB.CONTINUOUS
        v = self.model.addVar(*args, vtype=type, lb=lb, ub=ub, name=name, **kwargs)
        return v

    def quicksum(self, expr):
        return self.gurobipy.quicksum(expr)

    def update(self):
        self.model.update()

    @staticmethod
    def var_name(var):
        return var.var_name
    
    def LinExpr(self, *args, **kwargs):
        return self.gurobipy.LinExpr(*args, **kwargs)

    def abssum(self, variables, coeffs=None):
        total = 0
        for i, v in enumerate(variables):
            coeff = 1 if coeffs is None else coeffs[i]
            absvar = self.add_var()
            total += absvar * coeff
            self.add_constr(absvar + v >= 0)
            self.add_constr(absvar - v >= 0)
        return total

    def solve(self, objective=None, method='min', threads=6, time_limit=4000, log_path=None, MIPGap=None): # ret obj val
        if objective is not None:
            self.objective = objective
        self.model.setObjective(
            self.objective,
            self.gurobipy.GRB.MINIMIZE if method == 'min' else self.gurobipy.GRB.MAXIMIZE
        )
        self.model.update()

        if time_limit: self.model.params.timeLimit = time_limit
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

    # def solveAll(self, objective, var, method='min'):
    #     status, opt_value = self.solve(objective, method)
    #     solutions = [frozenset(a for a, y in var.items() if round(self.get_value(y)) > 0)]
    #     solutions += get_all_solutions(
    #         self, var, opt_value,
    #         candidates=solutions[0], iteration=0, mem=list()
    #     )
    #     solutions = [tuple(sorted(y for y in x)) for x in solutions]

    #     solutions = list(set(solutions))
    #     return status, opt_value, solutions

    def get_value(self, var):
        try:
            if var.vtype == self.gurobipy.GRB.BINARY:
                return round(var.x) > 0
            elif var.vtype == self.gurobipy.GRB.INTEGER:
                return int(round(var.x))
            else:
                return var.x
        except AttributeError:
            return var.getValue()

    def changeUb(self, var, ub):
        var.ub = ub
        self.model.update()
