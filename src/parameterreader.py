'''
Created on Aug 12, 2011

@author: ivihernandez
'''
#standard imports
import collections
import xml.etree.cElementTree as ET
import sys
import os
#non standard imports

#ivan's imports

class ParameterReader:
    """
        This function is in charge of reading the parameters
        for running the experiments
        
    """
    def __init__(self, file):
        """
            Constructor 
            @param file: configuration file for the experiments to run 
        """
        self.datasetFiles = []
        self.objectiveTypes = []
        self.objectiveFunctions = []
        tree = ET.ElementTree(file=file)
        for elem in tree.iter():
            
            if elem.tag == 'inputFiles':
                self.numberOfInputFiles = int(elem.attrib['number'])
                for subelem in elem:
                    filePath = os.path.normpath(subelem.text)
                    if subelem.tag == 'demandCenters':
                        self.demandCentersFilePath = filePath
                    elif subelem.tag == 'distributionCenters':
                        self.distributionCentersFilePath = filePath
                    else:
                        self.matrixFilePath = filePath
                    
                    self.datasetFiles.append(filePath)
            if elem.tag == 'capacitated':
                if elem.text.lower() == 'false':
                    self.capacitated = False
                else:
                    self.capacitated = True
            if elem.tag == 'objectiveFunctions':
                self.numberOfObjectiveFunctions = int()
            if elem.tag == 'populationSize':
                self.popSize = int(elem.text)
            if elem.tag == 'generations':
                self.generations = int(elem.text)
            if elem.tag == 'runs':
                self.runs = int(elem.text)
            if elem.tag == 'objective':
                
                for subelem in elem:
                    #print "subelem",subelem.tag
                    if subelem.tag == 'type':
                        #print "the text",elem.text.lower()
                        if subelem.text.lower() == 'minimize':
                            self.objectiveTypes.append(False)
                        else:
                            self.objectiveTypes.append(True)
                    if subelem.tag == 'function':
                        self.objectiveFunctions.append(subelem.text)
        self.weighted = True
    
    def get_matrix_file_path(self):
        return self.matrixFilePath
    def get_demand_centers_file_path(self):
        return self.demandCentersFilePath
    def get_distribution_centers_file_path(self):
        return self.distributionCentersFilePath
    def get_input_files(self):
        return self.datasetFiles
    def get_number_of_input_files(self):
        return self.numberOfInputFiles
    def get_number_of_objective_functions(self):
        return len(self.objectiveFunctions)
    def get_objective_functions_types(self):
        return self.objectiveTypes
    def get_objective_functions(self):
        return self.objectiveFunctions
    def is_capacitated(self):
        return self.capacitated
    def get_population_size(self):
        return self.popSize
    def get_generations(self):
        return self.generations
    def get_runs(self):
        return self.runs
    def is_weighted(self):
        return self.weighted
    def get_objective_parameters(self):
        return self.objectiveParameters
    def show(self):
        print "number of input files ", self.get_number_of_input_files()
        print "files used as input files ", self.get_input_files()
        print "population ", self.get_population_size()
        print "generations ", self.get_generations()
        print "runs ", self.get_runs()
        print "capacitated case? ", self.is_capacitated()
        print "weighted case? ", self.is_weighted()
        print "number of objective functions ", self.get_number_of_objective_functions()
        print "objective functions types ", self.get_objective_functions_types()
        print "objective functions", self.get_objective_functions()
        #print self.objective_parameters
        """
        print "parameters per objective function ",
        for key, elem in self.objectiveParameters.iteritems():
            print key, elem
        print 
        """
if __name__ == '__main__':
    folder = r'C:\Users\ivihernandez\Documents\PhD-systems-engineering\projects\facility-location\src\tests-folder'
    fileName = 'swain-total-distance.xml'
    filePath = os.path.join(folder, fileName)
    parameterReader = ParameterReader(filePath)
    parameterReader.show()