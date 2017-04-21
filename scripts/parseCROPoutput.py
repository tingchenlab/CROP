import sys
import csv
import os
from os.path import basename

"""
Author: Henry Yichao Dong
Email : yichaodo@usc.edu
This script parses the CROP output files and extracts the
cluster information in a clean format as a Cluster object

Usage:
python3 parseCROPoutput.py input.fasta

This scripts outputs two files, one with a simpler format, one with more details

The simpler output file should look like this:

1000 reads, 100 clusters identified
8 5 1,2,3,8,10
(followed by 99 other lines)
Line one is the summary
Starting from line 2, each line shows
centerID clusterSize member1ID,member2ID,member3ID,...

The detailed output file (has the additional center sequence info) should look like this:

1000 reads, 100 clusters identified
8 5 1,2,3,8,10
CTACTGGGATACGTCGGTTATCATC
(Followed by another 198 lines, since every cluster is represented by 2 lines)

More specific outputs can be customized by printing different member fields
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
        print(self.name + ' ' + str(self.size) )
        print(self.seq)
        print(self.members)


def outputSimpleFile(simpleFname, clusterList):
    totalNumSeqs = 0
    for cluster in clusterList:
        totalNumSeqs += cluster.size

    with open(simpleFname, "w") as outFile:
        outFile.write(str(totalNumSeqs) + ' reads, ' + str(len(clusterList)) + ' clusters identified\n')
        for cluster in clusterList:
            outFile.write(cluster.name + ' ' + str(cluster.size) + ' ')
            membersStr = cluster.members[0]
            for i in range (1, len(cluster.members), 1):
                membersStr += ',' + cluster.members[i]
            outFile.write(membersStr + '\n')


def outputDetailFile(detailFname, clusterList):
    totalNumSeqs = 0
    for cluster in clusterList:
        totalNumSeqs += cluster.size

    with open(detailFname, "w") as outFile:
        outFile.write(str(totalNumSeqs) + ' reads, ' + str(len(clusterList)) + ' clusters identified\n')
        for cluster in clusterList:
            outFile.write(cluster.name + ' ' + str(cluster.size) + ' ')
            membersStr = cluster.members[0]
            for i in range (1, len(cluster.members), 1):
                membersStr += ',' + cluster.members[i]
            outFile.write(membersStr + '\n')
            outFile.write(cluster.seq + '\n')


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

    #outputs two files, one has simple output and one with more details
    simpleFname = baseFname + ".simple.txt";
    detailFname = baseFname + ".detail.txt";
    outputSimpleFile(simpleFname, clusterList)
    outputDetailFile(detailFname, clusterList)


if __name__ == "__main__":
    main()
