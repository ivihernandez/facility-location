'''
Created on Feb 22, 2011

@author: ivihernandez
'''
#standard imports
import os

#non standard imports
import networkx as nx

#ivan imports
import myutils

def compare_nodes(a,b):
    (idA, attributesA) = a
    (idB, attributesB) = b
    #print idA, attributesA
    idA = int(idA)
    idB = int(idB)
    return idA - idB
    #return int(idA) - int(idB)

class GraphMLReader(object):
    '''
    classdocs
    '''
    """
    def write_mapping(self, graph, outputFile):
        distributionFile = open(outputFile, "w")
        distributionFile.write(myutils.get_xml_header())
        size = graph.number_of_nodes()
        mystr = '<pairing points=\"' + str(size) + '\">' + os.linesep
        distributionFile.write(mystr)
        for index, node in enumerate(graph.nodes(data=True)):
            (myid, _) = node
            mystr = '<point id=\"' + str(myid) + '\" '
            mystr += ' index=\"' + str(index) + '\">'  
            mystr += os.linesep
             
            mystr += '</point>' + os.linesep
            distributionFile.write(mystr) 
        mystr = "</pairing>" + os.linesep
        distributionFile.write(mystr)
        distributionFile.close()
    """ 
    def __init__(self, fileNameDistribution, fileNameDemand, resultsFolder=None):
        '''
        Constructor
        '''
        
        self.distribution = nx.read_graphml(fileNameDistribution, node_type=int)
        self.demand = nx.read_graphml(fileNameDemand, node_type=int)
        
        
        
        #write the mapping of distribution
        baseName = os.path.basename(fileNameDistribution)
        distributionFileMapName = baseName[0:-4] + "-distribution-mapping.xml"
        distributionFilePath = os.path.join(resultsFolder, distributionFileMapName)
        #self.write_mapping(self.distribution, distributionFilePath)
        #write the mapping of demands
        baseName = os.path.basename(fileNameDemand)
        demandFileMapName = baseName[0:-4] + "-demand-mapping.xml"
        demandFilePath = os.path.join(resultsFolder, demandFileMapName)
        """
        if resultsFolder != None:
            self.write_mapping(self.distribution, distributionFilePath)
            self.write_mapping(self.demand, demandFilePath)
        """
        
         
        
    def get_distribution(self):
        return self.distribution
    def get_demand(self):
        return self.demand
    def get(self):
        return (self.distribution, self.demand)
        
if __name__ == '__main__':
    print 'Reader started'
    distributionFileFolder = r'C:\Users\ivihernandez\Documents\PhD-systems-engineering\projects\facility-location\src\inputs'
    distributionFileName = 'swain-blue-yellow.xml'
    distributionFilePath = os.path.join(distributionFileFolder, distributionFileName)
    demandFilePath = distributionFilePath
    reader = GraphMLReader(distributionFilePath, demandFilePath, "temp")
    
    