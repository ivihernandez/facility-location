'''
Created on Aug 16, 2012

@author: ivihernandez
'''
#general imports
import random
import sys
import time
import os
import shutil
import copy
#non standard imports
import inspyred

#ivan's imports
import parameterreader
import facilitylocationproblem
import myutils
import psda
import mopsda
import nsga2

def multiple_runs(algorithmName, prng, problem, runs, pop, gen, step=None):
    archiver = inspyred.ec.archivers.best_archiver
    
    solutions = []#store several solutions
    total_population = []#store the population for the final pareto
    
    myArchive = []#store the pareto of the several runs
    
    
    maxTime = -1
    minTime = sys.maxint
    totalTime = 0
    for i in xrange(1, runs + 1):
        print algorithmName, "run #",i,
        start = time.clock()
        seeds = myutils.generate_seeds(problem, prng)
        if algorithmName == "MOPSDA":
            sol = psda.psda(prng=prng, popSize=pop, generations=gen,  problem=problem, seeds=seeds)
        elif algorithmName == "NSGA2":
            sol = nsga2.nsga2(prng=prng, popSize=pop, generations=gen,  problem=problem, seeds=seeds)
        else:
            raise(Exception("Unrecognized evolutionary algorithm. Try NSGA2 or MOPSDA"))
        #psda_sol = psda_basic(prng = prng, pop_size = pop, generations = gen,  problem = problem, initial_value = 0.5)
        end = time.clock()
        diff = end - start
        print diff
        #convert to minutes
        diff = diff/60.0
        totalTime += diff
        maxTime = max(maxTime, diff)
        minTime = min(minTime, diff)
        solutions.append(sol)
        for ind in sol.archive:
            total_population.append(ind)
        
    
    
    myArchive = archiver(random=prng, population=list(total_population), archive=list(myArchive), args=None)
    
    if algorithmName == 'MOPSDA':
        ea = mopsda.MOPSDA(prng)
    elif algorithmName == 'NSGA2':
        ea = inspyred.ec.emo.NSGA2(prng)
        
    ea.archive = myArchive 
    
    avgTime = totalTime/ float(runs) 
    totalEvaluations = problem.evaluations
    resultsFolder = problem.get_results_folder()
    
    
    attributes = {}
    attributes["totalTime"] = totalTime
    attributes["avgTime"] = avgTime
    attributes["maxTime"] = maxTime
    attributes["minTime"] = minTime
    attributes["totalEvaluations"] = totalEvaluations
    
    paretoOutputFilePath = os.path.join(resultsFolder,algorithmName+'-pareto.xml')
    #objs = problem.get_objective_functions()
    if step == None:
        step = ""
    else:
        step = str(step)
        paretoOutputFilePath = os.path.join(resultsFolder,algorithmName+"-"+step+'-pareto.xml')
    myutils.log_solution(paretoOutputFilePath, 
                         ea.archive, 
                         algorithmName, 
                         problem.get_objective_functions(),
                         problem.get_objective_functions_types(),
                         attributes)
    
    #sort population
    ea.archive.sort(cmp=myutils.individual_compare)
    ea.attributes = attributes
    
    return ea

def run_experiment(experimentFile, doPlot=True, prng=None ):
    if prng is None:
        prng = random.Random()
        myseed = time.time()
        myseed = 123
        prng.seed(myseed) 
    
    
    #read the parameters
    parameterReader = parameterreader.ParameterReader(experimentFile)
    parameterReader.show()
    
    #remove the folder with the results (if exists)
    resultsFolder = 'results '
    mydate = time.strftime("%a %d %b %Y %H %M %S", time.gmtime())
    resultsFolder += mydate
    
    try:
        shutil.rmtree(resultsFolder)
    except:
        pass
   
    #create a new folder for storing the results
    try:
        os.mkdir(resultsFolder)
    except:
        pass
    
    src = experimentFile
    fileName = os.path.basename(experimentFile)
    dst = os.path.join(resultsFolder, fileName) 
    shutil.copy(src, dst)
    #inputFolder = "inputs"
    
    pop = parameterReader.get_population_size()
    gen = parameterReader.get_generations()
    runs = parameterReader.get_runs()
    
    distributionCenters = parameterReader.datasetFiles[0]
    demandCenters = parameterReader.datasetFiles[1]
    args = [distributionCenters, demandCenters]
      
    
    
    innerArchive = None
    consideringFailures = False
    
    
    #NSGA-II
    if parameterReader.is_capacitated():
        problem = facilitylocationproblem.CapacitatedFacilityLocationProblem(args, prng, parameterReader, resultsFolder,"NSGA-II")
    else:
        problem = facilitylocationproblem.UncapacitatedFacilityLocationProblem(args, prng, parameterReader, resultsFolder,"NSGA-II")
    nsga2SolNoFailure = multiple_runs("NSGA2",prng, problem, runs, pop, gen, step=1)
    
    
    if "get_total_distance_after_interdiction" in problem.get_objective_functions():
        consideringFailures = True
    if "get_median_distance_after_interdiction" in problem.get_objective_functions():
        consideringFailures = True
    if "get_maximum_distance_after_interdiction" in problem.get_objective_functions():
        consideringFailures = True
    
    nsga2NoFailureTime = nsga2SolNoFailure.attributes["avgTime"]
    nsga2NoFailureEvaluations = nsga2SolNoFailure.attributes["totalEvaluations"]
    nsga2FailureTime = 0
    nsga2FailureEvaluations = 0
    if consideringFailures:
        inner   = facilitylocationproblem.InnerProblem(problem, parameterReader, nsga2SolNoFailure)
        innerArchive = inner.get_pareto()
        
        nsga2SolFailure = copy.deepcopy(nsga2SolNoFailure)
        nsga2SolFailure.archive = innerArchive
        paretoOutputFilePath = os.path.join(resultsFolder, "NSGA2-pareto-2.xml")
        myutils.log_solution(paretoOutputFilePath, 
                             nsga2SolFailure.archive, 
                             "NSGA2", 
                             problem.get_objective_functions(),
                             problem.get_objective_functions_types(),
                             nsga2SolFailure.attributes#nsga2SolFailure.attributes
                             )
        #nsga2FailureTime = nsga2SolFailure.attributes["avgTime"]
        #nsga2FailureEvaluations = nsga2SolFailure.attributes["totalEvaluations"]
        nsga2FailureTime = inner.attributes["time"]
        nsga2FailureEvaluations = inner.attributes["totalEvaluations"]
        
        
    #PSDA
    
    if parameterReader.is_capacitated():
        problem = facilitylocationproblem.CapacitatedFacilityLocationProblem(args, prng, parameterReader, resultsFolder,"MOPSDA")
    else:
        problem = facilitylocationproblem.UncapacitatedFacilityLocationProblem(args, prng, parameterReader, resultsFolder,"MOPSDA")
    
   
    psdaSolNoFailure = multiple_runs("MOPSDA",prng, problem, runs, pop, gen, step=1)
    psdaNoFailureTime = psdaSolNoFailure.attributes["avgTime"]
    psdaNoFailureEvaluations = psdaSolNoFailure.attributes["totalEvaluations"]
    psdaFailureTime = 0
    psdaFailureEvaluations = 0
     
    
    if consideringFailures:
        inner   = facilitylocationproblem.InnerProblem(problem, parameterReader, psdaSolNoFailure)
        innerArchive = inner.get_pareto()
        
        psdaSolFailure = copy.deepcopy(psdaSolNoFailure)
        psdaSolFailure.archive = innerArchive
        paretoOutputFilePath = os.path.join(resultsFolder, "MOPSDA-pareto-2.xml")
        myutils.log_solution(paretoOutputFilePath, 
                             psdaSolFailure.archive, 
                             "MOPSDA", 
                             problem.get_objective_functions(),
                             problem.get_objective_functions_types(),
                             psdaSolFailure.attributes
                             )
        
        psdaFailureTime = inner.attributes["time"]
        psdaFailureEvaluations = inner.attributes["totalEvaluations"]
    #combine solutions without failure
    combinedSolNoFailure = myutils.combine_paretos(prng, problem, [psdaSolNoFailure, nsga2SolNoFailure])
    paretoOutputFilePath = os.path.join(resultsFolder, "combined-pareto-1.xml")
    myutils.log_solution(paretoOutputFilePath, 
                         combinedSolNoFailure.archive, 
                         "MOPSDA and NSGA2", 
                         problem.get_objective_functions(),
                         problem.get_objective_functions_types(),
                         combinedSolNoFailure.attributes,
                         
                         )
    
    
    #combine PSDA and NSGA-II Solutions
    totalEvaluations = psdaFailureEvaluations + psdaNoFailureEvaluations + nsga2FailureEvaluations + nsga2NoFailureEvaluations
    totalTime = psdaFailureTime + runs * psdaNoFailureTime + nsga2FailureTime + runs * nsga2NoFailureTime
    attributes = {}
    attributes["totalCombinedTime"] = totalTime
    attributes["totalCombinedEvaluations"] = totalEvaluations
    attributes["seed"] = myseed
    if consideringFailures:
        combinedSol = myutils.combine_paretos_considering_failures(prng, problem, [psdaSolFailure, nsga2SolFailure])
    else:
        combinedSol = myutils.combine_paretos(prng, problem, [psdaSolNoFailure, nsga2SolNoFailure])
    paretoOutputFilePath = os.path.join(resultsFolder, "combined-pareto-2.xml")
    
    myutils.log_solution(paretoOutputFilePath, 
                         combinedSol.archive, 
                         "MOPSDA and NSGA2", 
                         problem.get_objective_functions(),
                         problem.get_objective_functions_types(),
                         combinedSolNoFailure.attributes,
                         attributes
                         )
    
    
