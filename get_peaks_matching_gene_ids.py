import sys
import math


peakAnnotFile = sys.argv[1]
geneListFile = sys.argv[2]
outFile = sys.argv[3]


distTssFilter = math.inf
if len(sys.argv) == 5:
    distTssFilter = int(sys.argv[4])


geneSet = set()
with open(geneListFile, 'r') as f:
    for line in f:
        geneSet.add(line.strip())

matchingPeaks = []        
with open(peakAnnotFile, 'r') as f:
    f.readline()
    for line in f:
        peakData = line.strip().split("\t")
        peakID = peakData[3]
        geneName = ""
        distTss = sys.maxsize
        if len(peakData)==6:
            geneName = peakData[5]
            distTss = int(peakData[4])
        if geneName in geneSet and abs(distTss) < distTssFilter:
            matchingPeaks.append(peakData[3])
            
with open(outFile, 'w') as f:
    for matchingPeak in matchingPeaks:
        f.write(matchingPeak)
        f.write("\n")