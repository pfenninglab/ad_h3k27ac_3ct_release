import gzip

if __name__=="__main__":
    inputFile = "ucsc_table_browser_ucsc_known_genes_all_genes.txt.gz"
    with gzip.open(inputFile, 'rt') as f:
        line = f.readline()
        for line in f:
            data = line.strip().split("\t")
            [chrom, strand, start, stop, exonCount, exonStarts, exonEnds, name] = data

            assert(exonStarts[-1]==",")
            assert(exonEnds[-1]==",")
            exonStarts = exonStarts[0:len(exonStarts)-1]
            exonEnds = exonEnds[0:len(exonEnds)-1]

            exonStarts = [int(val) for val in exonStarts.split(",")]
            exonEnds = [int(val) for val in exonEnds.split(",")]
            
            exonLengths = [exonEnds[i] - exonStarts[i] for i in range(len(exonEnds))]

            blockStarts = [val - int(start) for val in exonStarts]
            
            blockStarts = [str(val) for val in blockStarts]
            exonLengths = [str(val) for val in exonLengths]            
            print("\t".join([chrom, start, stop, name, "0", strand, start, stop, "0,0,0", exonCount, ",".join(exonLengths), ",".join(blockStarts)]))
            
