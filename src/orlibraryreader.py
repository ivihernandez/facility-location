'''
Created on Jul 22, 2011

@author: ivihernandez
'''

import networkx as nx

class ORLibraryReader(object):
    """
        This reader can retrieve this type of documents:
        
    """
    def __init__(self, fileName):
        file = open(fileName, "r")
        line = file.readline()
        values = line.strip().split()
        m = int(values[0])   #number of potential facilities
        n = int(values[1])   #number of demand centers
        
        self.distributionNodes  = nx.Graph()
        self.demandNodes = nx.Graph()
        #map where distances are stored
        #key = (demandCenter, distributionCenter)
        #value = distance
        self.distanceMap = {} 
        for i in xrange(m):
            line = file.readline()
            values = line.strip().split()
            capacity = values[0]
            cost = float(values[1])
            thedict = {}
            thedict["cost"] = cost
            thedict["capacity"] = capacity
            id = i
            
            self.distributionNodes.add_node(id, thedict)
        for j in xrange(n):
            line = file.readline()
            values = line.strip().split()
            demand = int(values[0])
            
            id = j
            thedict = {}
            thedict["demand"] = demand
            
            itemsRead = 0
            while True:
                if itemsRead >= m:
                    break
                line = file.readline()
                values = line.strip().split()
                for i in xrange(len(values)):
                    self.distanceMap[(j, itemsRead)] = float(values[i])
                    
                    self.demandNodes.add_node( id, thedict)
                    itemsRead += 1
            
        print self.demandNodes
        print self.distributionNodes
    def get_distribution(self):
        """
            retur the list of distribution nodes (no coordinates included)
        """
        return self.distributionNodes
    def get_demand(self):
        """
            retur the list of demand nodes (no coordinates included)
        """
        return self.demandNodes
    def get_distance_map(self):
        """
            return the matrix with the distances from demandCenters to distributionCenters
        """
        return self.distanceMap
            