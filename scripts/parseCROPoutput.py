import sys
import csv
import os
from os.path import basename

"""
Author: Henry Yichao Dong
Email : yichaodo@usc.edu
This script parses the CROP output files and extracts the
cluster information in a clean format as a Cluster object

The specific output can be customized by printing different member fields
of the cluster objects stored in the clusterList variable
"""

class Cluster:
    def __init__(self, repName, clusterSize, repSequence = ""):
        self.name = repName
        self.size = clusterSize
        self.seq = repSequence
        self.stdDev = 0.0
        self.members = [] #a list of member ids

    def print(self):
        print(self.name + ' : ' + str(self.size) )
        print(self.seq)
        print(self.members)


def main():
    numCmdArgs = len(sys.argv)
    if numCmdArgs != 2:
        print ('Correct Usage:')
        print ('python3 parseCROPoutput.py input.fasta')
        sys.exit(0)

    baseFname= sys.argv[1]
    totalNumClusters = 0;
    clusterList = [] #stores all of the Cluster objects
    clusterFname = baseFname + '.cluster'
    clusterFastaFname = clusterFname + '.fasta'
    clusterListFname = clusterFname + '.list'

    #First process the .cluster file to get the cluster sizes
    with open(clusterFname) as clusterFile:
        totalNumClusters = int(next(clusterFile))
        for line in clusterFile:
            if not line:
                continue
            lineValues = line.split('\t')
            cluster = Cluster (lineValues[0], int(lineValues[1]))
            cluster.stdDev = float(lineValues[2])
            clusterList.append(cluster)
    assert len(clusterList) == totalNumClusters

    #Then process the .cluster.fasta file to include the rep seq in each cluster
    with open(clusterFastaFname) as clusterFastaFile:
        idx = 0
        for line in clusterFastaFile:
            if line[0] == '>':
                assert line[1:].strip('\n') == clusterList[idx].name
            else:
                clusterList[idx].seq = line.strip('\r\n ')
                idx += 1

    #Lastly process the .cluster.list file to add cluster members IDs
    with open(clusterListFname) as clusterListFile:
        idx = 0
        for line in clusterListFile:
            members = line.split('\t')[1].split(',')
            members[-1] = members[-1].strip()
            clusterList[idx].members = members
            assert clusterList[idx].size == len(members)
            idx += 1

    #sort clusters by decreasing size
    clusterList.sort(key=lambda cluster: cluster.size, reverse=True)

    #print out the info of each cluster to STDOUT
    for cluster in clusterList:
        cluster.print()
        print()

    totalNumSeqs = 0
    for cluster in clusterList:
        totalNumSeqs += cluster.size
    print('Total Number of Sequences = ' + str(totalNumSeqs) )
    print('Total Number of Clusters  = ' + str(totalNumClusters) + '\n')

if __name__ == "__main__":
    main()
