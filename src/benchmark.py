'''
Created on Feb 15, 2011

@author: ivihernandez
'''
import math
import random
from inspyred import ec


#from ecspy import benchmarks

class Benchmark(object):
    """Defines a global optimization benchmark problem.
    
    This abstract class defines the basic structure of a global
    optimization problem. Subclasses should implement the generator 
    and evaluator methods for a particular optimization problem, 
    which can be used with ECsPy evolutionary computations. 
    
    In addition to being used with evolutionary computations, subclasses
    of this are also callable. the arguments passed to such a call are
    combined into a list and passed as the single candidate to the 
    evaluator method. The single calculated fitness is returned.
    
    Public Attributes:
    
    - *dimensions* -- the number of inputs to the problem
    - *objectives* -- the number of outputs of the problem (default 1)
    - *bounder* -- the bounding function for the problem (default None)
    - *maximize* -- whether the problem is one of maximization (default True)
    
    """
    def __init__(self, dimensions, objectives=1):
        self.dimensions = dimensions
        self.objectives = objectives
        self.bounder = None
        self.maximize = True
        
    def __str__(self):
        if self.objectives > 1:
            return '%s (%d dimensions, %d objectives)' % (self.__class__.__name__, self.dimensions, self.objectives)
        else:
            return '%s (%d dimensions)' % (self.__class__.__name__, self.dimensions)
        
    def __repr__(self):
        return '%s' % self.__class__.__name__
    
    def generator(self, random, args):
        raise NotImplementedError
        
    def evaluator(self, candidates, args):
        raise NotImplementedError
        
    def __call__(self, *args):
        candidate = [a for a in args]
        fit = self.evaluator([candidate], {})
        return fit[0]