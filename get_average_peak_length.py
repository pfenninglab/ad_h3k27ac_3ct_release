import argparse

def getData(line):
    split = line.split()
    chr = split[0]
    start = int(split[1])
    end = int(split[2])
    
    return (chr, start, end)

def calculatePeakLength(peakFile):
    inf = open(peakFile, 'r')
            
    totalPeakLength = 0
    numPeaks = 0
    
    for line in inf:
        (chr, start, end) = getData(line)
        totalPeakLength += end-start
        numPeaks += 1
    
    avgPeakLength = 1.0*totalPeakLength/numPeaks
    inf.close()
    return avgPeakLength, totalPeakLength


if __name__=="__main__":
    parser = argparse.ArgumentParser(description = "Calculate average region length in a bed file")
    parser.add_argument('-i', '--peak-file', help =  "Input file containing peaks", required=True)
    args = parser.parse_args()    
    print calculatePeakLength(args.peak_file)
