import argparse

def getData(line):
    split = line.split()
    chr = split[0]
    start = int(split[1])
    end = int(split[2])
    
    return (chr, start, end)

def namePeaks(peakFile, out, prefix, padding):
    inf = open(peakFile, 'r')
    outf = open(out, 'w')
    
    count = 0
    for line in inf:
        count+=1
        (chr, peakStart, peakEnd) = getData(line)
        outf.write("\t".join([chr, str(peakStart), str(peakEnd), prefix+str(count).zfill(padding)]))
        outf.write("\n")
        
    
    inf.close()
    outf.close()

if __name__=="__main__":
    parser = argparse.ArgumentParser(description = "Script to name peaks in a bed file")
    parser.add_argument('-i', '--peak-file', help = "input peak bed file", required = True)
    parser.add_argument('-o', '--out', help = "output peak bed file", required = True)
    parser.add_argument('-p', '--prefix', help = "prefix of the label", required = True)
    parser.add_argument('-w', '--padding', help = "length of the number for label", required = True, type = int)
    
    args = parser.parse_args()
    
    namePeaks(args.peak_file, args.out, args.prefix, args.padding)
    
