#standard imports
import math
import random
import sys
import copy
import numpy
import os
import collections
import time
#non standard imports
#import inspyred.ec
import inspyred
from inspyred.ec import terminators
from inspyred.ec import variators
from inspyred.ec import emo
#from inspyred.ec import Bounder

#import benchmark
#ivan's imports

import graphmlreader
import mopsda
import psda
import nsga2
import myutils
import candidatelist
#inspyred.benchmarks
#class FacilityLocationProblem(benchmark.Benchmark):
class FacilityLocationProblem():
    def __init__(self, inputs, prng, parameterReader, resultsFolder,embeddedEa):
        """
            @param inputs: list of files to work with.
            
            If there is only one param, the file is assumed to already contain
            the weighted distances 
            (http://people.brunel.ac.uk/~mastjjb/jeb/orlib/uncapinfo.html) 
            
            If there are two paramenters, the files are assumed to be in the GraphML format.
            In thi case, the file contains the x and y coordintes, and may contain
            additional attributes such as the demand or the cost
            (http://graphml.graphdrawing.org/primer/graphml-primer.html)
        """
        self.inputs = inputs #newline consider eliminating the non graphml case
        if len(inputs) == 2:
            fileNameDistribution = inputs[0]
            fileNameDemand = inputs[1]
            reader = graphmlreader.GraphMLReader(fileNameDistribution,fileNameDemand, resultsFolder)
            #self.distanceMap = {}
            self.readerType = "GRAPHML"
        else:
            raise Exception("wrong number of input files, provide distribution and demand centers") 
        
        self.distribution = reader.get_distribution()
        self.demand = reader.get_demand()
        self.maximize = True
        self.objectiveTypes = None
        self.debug = False
        self.evaluations = 0 #counter for the number of evaluations performed
        self.capacitated = None #capacitated or uncapacitated case?
        self.weighted = None #capacitated or uncapacitated case?
        self.minDistance = None #minimum of all the minimum distances
        self.maxDistance = None #maximum of all the minimum distances
        self.show = False #do not show the figures while running
        self.myFormat = "pdf"
        self.chromosomeType = None
        self.maxInterdiction = "get_maximum_distance_after_interdiction"
        self.medianInterdiction = "get_median_distance_after_interdiction"
        self.totalInterdiction = "get_total_distance_after_interdiction"
        self.set_objective_functions(parameterReader.get_objective_functions())
        self.set_capacitated(parameterReader.is_capacitated())
        self.set_weighted(parameterReader.is_weighted())
        self.set_objective_types(parameterReader.get_objective_functions_types())
        self.set_results_folder(resultsFolder)
        self.set_embedded_ea_algorithm(embeddedEa)
        self.set_prng(prng)
        if self.debug:
            print "in Facility:__init__"
    def get_distribution_centers(self):
        return self.distribution
    def get_demand_centers(self):
        return self.demand
    def get_objective_functions_types(self):
        return self.objectiveTypes
    def get_objective_functions(self):
        return self.objectiveFunctions
    def set_embedded_ea_population_size(self, pop):
        self.embedded_ea_population_size = pop
    def set_embedded_ea_generation_runs(self, gen):
        self.embedded_ea_generation_runs = gen
    def is_capacitated(self):
        """
            return True if running the capacitated case 
            (calculate_capacitated_distance).
            
            return False if running uncapacitated case
        """
        return self.capacitated
    
    def get_distance(self, demandCenter, distributionCenter):
        """return the weighted distance between two nodes
            @param: (a)=a node of the graph with attributes x and y (demandCenter)
                    (b)=a node of the graph with attributes x and y (distributionCenter)
            @return: weighted distance between a and b
        """ 
        a = self.demand.node[demandCenter]
        b = self.distribution.node[distributionCenter]
        distance = math.sqrt(  (a["x"]-b["x"])*(a["x"]-b["x"]) + (a["y"]-b["y"])*(a["y"]-b["y"]) ) 
        try:
            demand = a["weight"]
        except:
            demand = 1
        weightedDistance = demand * distance
         
            
        return  weightedDistance
    def set_prng(self, prng):
        """
            @param prng: the prng is the random number generator
            This is required, in order to call a EA from another EA
        """
        self.prng = prng
    def set_objective_parameters(self, objectiveParameters):
        """
            @arg objective_parameters: a dictionary with the
                name of the objective function followed
                by a list of parameters
        """
        self.objectiveParameters = objectiveParameters
    def set_objective_functions(self, objectiveFunctions):
        """
            #arg objective_function: list of objective functions to be called
        """
        self.objectiveFunctions = objectiveFunctions
        
        
    def set_results_folder(self, folderName):
        """
            @arg folderName: name of the folder were the results
                are going to be stored.
            @desc: This function allows the class Problem to know were
                put the results
        """
        self.resultsFolder = folderName
    def set_capacitated(self, capacitated):
        """
            Determine if the capacity of the dataset is going to be taken into account
        """
        self.capacitated = capacitated
    def set_weighted(self, weighted):
        """
            Determine if the weights of the dataset are going to be used or not
        """
        self.weighted = weighted
        if not self.weighted:
            #set the weight equal to zero if there are no weights
            for node in self.demand:
                self.demand.node[node]["weight"] = 1
    def get_results_folder(self ):
        """
            @desc: This function return the name of the folder
                were the results are going to be stored
        """
        return self.resultsFolder
    def get_prng(self):
        """
            @return : the random number generator 
        """
        return self.prng
    def set_embedded_ea_algorithm(self, algorithmName):
        """
            @param algorithmName: name of the algorithm to use in the 
            vulnerability function where an EA calls another EA.
            
            With this you can call PSDA from NSGA-II or from PSDA.
            Conversely, you can call NSGA-II from NSGA-II or from PSDA
            
            valid function names are: MOPSDA or NSGA-II
        """
        self.embeddedEaAlgorithm = algorithmName
    def get_embedded_ea_algorithm(self):
        """
            @return name: the name of the algorithm use in the routine
            get_vulnerability.
            
            That function allows an EA to call another EA.
            valid function names are: MOPSDA or NSGA-II
        """
        return self.embeddedEaAlgorithm
    
    def set_objective_types(self, objectiveTypes):
        """
            Set if each objective is going to be minimized or maximized
        """
        self.objectiveTypes = objectiveTypes
    
    def get_total_distance(self, ind):
        """
            Objective Function
            
            This functions calculates the distance between demand centers and 
            distribution centers. 
            The function knows wether the distribution centers have capacity or not.
            The function knows wether the demand centers have demand or not.
        """
        
        self.compute_distance(ind)
        return self.totalDistance    
    def get_maximum_distance(self, ind):
        """
            @param ind: chromosome
            @return: the maximum distance 
            Objective function
            
            This functions calculates the maximum distance between 
            all the distances, from demand centers to 
            distribution centers.
             
            The function calls count_distance, so that it knows: 
            The function knows wether the distribution centers have capacity or not.
            The function knows wether the demand centers have demand or not.
        """
        self.compute_distance(ind)
        return self.maxDistance
    def get_median_distance(self, ind):
        self.compute_distance(ind)
        
        ret = numpy.median(self.distances)
        return ret
    def get_embedded_ea_population_size(self):
        return self.embeddedEaPopulationSize
    def get_embedded_ea_generation_runs(self):
        return self.embeddedEaGenerationRuns
    
    
    
    def get_map(self,ind):
        """
            @desc: This function returns the mapping of which demand center
            goes to which distribution center:
            
            
            @input: the binary string representing the individual
            @output: a dictionary of the form: 
            Distribution Center -> List of Demand Centers assigned to it.
        """
        self.compute_distance(ind)            
        return self.map
    
    def evaluator(self, candidates, args):
        
        
        ret = None
        for name in self.objectiveFunctions:
            if name == "get_maximum_distance_after_interdiction":
                ret = self.maximum_distance_after_interdiction_evaluator(candidates, args)
            if name == "get_total_distance_after_interdiction":
                ret = self.total_distance_after_interdiction_evaluator(candidates, args)
            if name == "get_median_distance_after_interdiction":
                ret = self.median_distance_after_interdiction_evaluator(candidates, args)
        
        if ret == None:
            ret = self.simple_evaluator(candidates, args) 
        
        
        return ret
    
    
    
        
    def simple_evaluator(self, candidates, args):
        """
                This evaluator is called from evaluator.
                Is the function to be used in case none of the objective
                functions imply to call another EA.
                
                The evaluator of the Super Class: FacilityLocationProblem
                will call this function.
        """
        
        fitness = []
        for c in candidates:
            scores = []
            for name in self.objectiveFunctions:
                #no extra parameters
                methodToCall = getattr(self, name)
                result = methodToCall(c)
                scores.append(result)
            self.evaluations += 1
            
            pareto = emo.Pareto(scores, self.objectiveTypes)
            newdict = self.map.copy()
            c.set_map(newdict)
            fitness.append(pareto)
        
        return fitness
    
    def total_distance_after_interdiction_evaluator(self, candidates, args):
        """
            Perform basic computations needed
            for computing the maximum distance after interdiction.
            This evaluator is called from evaluator.
                
        """
        
        fitness = []
        newcandidates = [] 
        for c in candidates:
            
            scores = []
            for name in self.objectiveFunctions:
                #only obtain first two objective function's values
                #failure will be obtained afterwards
                if name == "get_total_distance_after_interdiction":
                        result = 0
                elif name == "get_interdicted_facilities":
                        result = 0
                else:
                    methodToCall = getattr(self, name)
                    
                    result = methodToCall(c)
                
                scores.append(result)
            self.evaluations += 1
            #add candidates and fitnes
            newdict = self.map.copy()
            c.set_map(newdict)
            #c.set_failures_chromosome()
            newcandidates.append(copy.deepcopy(c))
            fitness.append(emo.Pareto(scores, self.objectiveTypes), )
                
                
        #ONLY WAY to MODIFY a LIST inside the function, array splicing
        candidates[:] = newcandidates
        
        return fitness
    
    def median_distance_after_interdiction_evaluator(self, candidates, args):
        """
            Perform basic computations needed
            for computing the maximum distance after interdiction.
            This evaluator is called from evaluator.
                
        """
        fitness = []
        newCandidates = [] 
        for c in candidates:
            
            scores = []
            for name in self.objectiveFunctions:
                if name == "get_median_distance_after_interdiction":
                        result = 0
                elif name == "get_interdicted_facilities":
                        result = 0
                else:
                    methodToCall = getattr(self, name)
                    result = methodToCall(c)
                
                scores.append(result)
            self.evaluations += 1
            #add candidates and fitnes
            newCandidates.append(copy.deepcopy(c))
            fitness.append(emo.Pareto(scores, self.objectiveTypes), )
                
                
        #ONLY WAY to MODIFY a LIST inside the function, array splicing
        candidates[:] = newCandidates
        return fitness
        
        

    def maximum_distance_after_interdiction_evaluator(self, candidates, args):
        """
            Perform basic computations needed
            for computing the maximum distance after interdiction.
            This evaluator is called from evaluator.
                
        """
        fitness = []
        newCandidates = [] 
        for c in candidates:
            
            scores = []
            for name in self.objectiveFunctions:
                if name == "get_maximum_distance_after_interdiction":
                        result = 0
                elif name == "get_interdicted_facilities":
                        result = 0
                else:
                    methodToCall = getattr(self, name)
                    result = methodToCall(c)
                
                scores.append(result)
            self.evaluations += 1
            #add candidates and fitnes
            newCandidates.append(copy.deepcopy(c))
            fitness.append(emo.Pareto(scores, self.objective_types), )
                
                
        #ONLY WAY to MODIFY a LIST inside the function, array splicing
        candidates[:] = newCandidates
        return fitness
    
class UncapacitatedFacilityLocationProblem(FacilityLocationProblem):
    def __init__(self, inputs, prng, parameterReader, resultsFolder,embeddedEA):
        
        self.bounder = inspyred.ec.Bounder(0, 1)
        FacilityLocationProblem.__init__(self, inputs, prng, parameterReader, resultsFolder,embeddedEA)
        
        self.dimensions = self.distribution.number_of_nodes()
        
        
        
        self.chromosomeType = "binary"     #the chromosome is binary
        
    
    def generator(self, random, args):
        retval = [random.choice([0, 1]) for _ in xrange(self.dimensions )]
        return candidatelist.Candidatelist(retval)
    def get_deployed_facilities(self,elems):
        """
            Objective function
            
            count the number of bits that are on in the chromosome
            @param: elems = list of integers between 1 and 0
            @return: integer (number of ones in elem)
        """
        counter = 0
        for x in elems:
            if x == 1:
                counter = counter + 1
        if self.debug:
            print "in count_ones",counter
        return counter
    def compute_distance(self, ind):
        """
            
            
            @desc: Assign each demand center to its nearest distribution center (unconstrained).
            Demand can always be satisfied.
            
            @param ind: chromosome
            @param distanceFunction: either get_weighted_distance or get_distance
            
            @return: total distance from each demand center to its nearest distribution center
            
            @post-condition: a map of Distribution center -> list of demand centers is created.
        """
        """
        if len(self.distribution.nodes(data=True))==48:
        
            self.distribution.remove_node("42")
            self.distribution.remove_node("44")
            self.distribution.remove_node("45")
            self.distribution.remove_node('54')
            self.distribution.remove_node('55')
            
            self.distribution.add_node('32', {"x":32,"y":45, "weight":5})
            self.distribution.add_node('40', {"x":50,"y":40, "weight":4})
            self.distribution.add_node('41', {"x":23,"y":22, "weight":4})
            self.distribution.add_node('51', {"x":27,"y":5, "weight":3})
            self.distribution.add_node('52', {"x":52,"y":24, "weight":2})
        """
        
        distribution = self.distribution.nodes(data=True)#graph, whose nodes represent distribution centers)
        demand = self.demand.nodes(data=True)#graph, whose nodes represent demand centers).
        
        
        
        
        self.map = collections.defaultdict(set)#map of distribution center-> list of demand centers
        for node in self.distribution:
            self.map[node] = set()
        self.minDistance = sys.maxint #minimum of all the minimum distances
        self.maxDistance = -1 #maximum of all the minimum distances
        self.distances = [] #all the distances from each demand center to its closest DC
        distance = 0
        self.totalDistance = 0
        closestIndex = None
        #self.distanceMap = collections.defaultdict(float) #(DistributionCenter, DemandCenter) -> Distance
        
        ones = myutils.count_ones(ind)
        self.debug = False
        
        
        for demIndex, demValues in demand:    
            closestDistance = float(sys.maxint)
            counter = -1
            
            for distIndex, distValues in distribution:
                
                counter += 1
                
                if ind[counter] == 0:
                    continue
                else:
                    pass
                
                
                
                temp = self.get_distance(demIndex, distIndex)
                #self.distanceMap[(distIndex,demIndex)] = temp
                if temp < closestDistance:
                    closestDistance = temp
                    closestIndex = distIndex
             
            if closestIndex == None:
                #there were no Distribution Centers deployed
                distance = sys.maxint
                self.minDistance = sys.maxint
                self.maxDistance = sys.maxint
                self.distances.append(sys.maxint)
            else:
                self.map[closestIndex].add(demIndex) 
                  
                distance = distance + closestDistance
                self.minDistance = min(self.minDistance, closestDistance)
                self.maxDistance = max(self.maxDistance, closestDistance)
                self.distances.append(closestDistance)
        
        
        self.medianDistance = numpy.median(self.distances)   
        self.totalDistance = distance
        
        self.debug = False
        
        if self.debug:
            print "number of demand points",len(demand)
            print "number of distribution points",len(self.map)
            print "distance",distance
            print "keys in the map",len(self.map.keys())
            i = 1
            for key, value in self.map.iteritems():
                print i,")", key, value
                i = i + 1
            self.debug = False

class CapacitatedFacilityLocationProblem(FacilityLocationProblem):
    pass

class TotalDistanceAfterInterdiction(UncapacitatedFacilityLocationProblem):
    def __init__(self, parentObject, individual):
        self.bounder = inspyred.ec.Bounder(0, 1)
        self.individual = individual
        self.chromosome = individual.candidate
        self.distribution = parentObject.distribution
        self.demand = parentObject.demand
        self.dimensions = parentObject.dimensions
        self.objectiveFunctions = parentObject.objectiveFunctions
        #maximize: Distance after interdiction
        #minimize: interdicted facilities
        self.objectiveTypes = [True, False]
        self.maximize = True
        self.debug = False
        self.evaluations = 0 #counter for the number of evaluations performed
        self.capacitated = parentObject.capacitated #capacitated or uncapacitated case?
        self.weighted = parentObject.weighted #capacitated or uncapacitated case?
        self.minDistance = None #minimum of all the minimum distances
        self.maxDistance = None #maximum of all the minimum distances
        self.show = False #do not show the figures while running
        self.myFormat = "pdf"
        self.chromosomeType = None
        self.evaluations = 0
    def evaluator(self, candidates, args):
        fitness = []
        #newCandidates = []
        for c in candidates:
            self.evaluations += 1
            scores = []
            interdictionStrategy = c
            interdictedChromosome = myutils.turn_off(self.chromosome, interdictionStrategy)
            index = -1
            for name in self.objectiveFunctions:
                index = index + 1
                if name == "get_total_distance_after_interdiction":
                        result = self.get_total_distance(interdictedChromosome)
                elif name == "get_interdicted_facilities":
                        result = myutils.get_interdicted_facilities(self.chromosome, interdictionStrategy)
                else:
                    continue
                scores.append(result)
            #c.set_interdicted_facilities(interdictedChromosome)
            #newCandidates.append(copy.deepcopy(c))
            
            fitness.append(emo.Pareto(scores, self.objectiveTypes), )
            
        #candidates[:] = newCandidates
        return fitness
    
    
    
    
    



class InnerProblem:
    
    def __init__(self, externalProblem, parameterReader, externalEA):
        start = time.clock()
        self.prng = externalProblem.prng
        self.archive = []
        
        self.population = []
        initialPop = []
        initialFits = []
        objective_types = externalProblem.objectiveTypes
        self.resultsFolder = externalProblem.resultsFolder
        
        self.myFormat = externalProblem.myFormat
        self.show = externalProblem.show
        self.failureEAs = []
        for externalSolution in externalEA.archive:
            #print "externalSolution", externalSolution
            
            requiresInnerEA = False
            if "get_total_distance_after_interdiction" in externalProblem.objectiveFunctions:
                innerProblem = TotalDistanceAfterInterdiction(externalProblem, copy.deepcopy(externalSolution))
                requiresInnerEA = True
            
            else:
                raise(Exception("My Own Exception: Unsuported Failure Problem"))
            #choose inner algorithm
            if requiresInnerEA and externalProblem.get_embedded_ea_algorithm() == "MOPSDA":
                innerEA = psda.psda(externalProblem.prng, parameterReader.popSize, parameterReader.generations, innerProblem, [])
            elif requiresInnerEA and externalProblem.get_embedded_ea_algorithm() == "NSGA-II":
                innerEA = nsga2.nsga2(externalProblem.prng, parameterReader.popSize, parameterReader.generations, innerProblem, [])
            else:
                innerEA = innerProblem.get_pareto()
                #raise(Exception("My Own Exception: Unsuported inner Algorithm"))
            
            totalEvaluations = innerProblem.evaluations
            for internalSolution in innerEA.archive:
                #print "internalSolution", internalSolution
                elem = copy.deepcopy(externalSolution.candidate)
                
                elem.set_chromosome_after_interdiction(internalSolution.candidate)
                newdict = innerProblem.map.copy()
                elem.set_map(newdict) 
                
                initialPop.append(elem)
                scores = []
                scores.append(externalSolution.fitness[0])
                scores.append(externalSolution.fitness[1])
                scores.append(internalSolution.fitness[0])
                scores.append(internalSolution.fitness[1])
                fit = emo.Pareto(scores, objective_types)
                initialPop.append(elem)
                initialFits.append(fit)
                ind = inspyred.ec.Individual(elem)
                
                
                ind.fitness = fit
                self.population.append(ind)
                
        
            #print "----"
        """
        print "population"
        for x in self.population:
            print x
        """
        
        self.archive = myutils.my_best_archiver(random=self.prng, population=list(self.population), archive=list(self.archive), args=None)
        """
        print "archive"
        for x in self.archive:
            print x
        """ 
        
        
        end = time.clock()
        diff = (end - start)/60.0
        self.attributes = {}
        self.attributes["time"] = diff
        self.attributes["totalEvaluations"] = totalEvaluations
        
        
    def get_pareto(self):
        """
            Return a pareto set, whose individuals are the original chromosome
        """
        return self.archive
    

