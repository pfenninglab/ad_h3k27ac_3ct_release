library(optparse) 


option_list = list( 
                    make_option(c("-m", "--matrix_counts"), type="character", 
                                                    help="file containing count matrix"),
                    make_option(c("-i", "--info_sample"), type="character", 
                                                    help="file containing sample info"),
                    make_option(c("-p", "--peak_info"), type="character", 
                                                    help="file containing peak info"),                                                     
                    make_option(c("-d", "--design"), type="character", default=NULL, 
                                                    help="design string for analysis"),
                    make_option(c("-v", "--variable_contrast"), type="character", default=NULL, 
                                                    help="output contrast variable"),
                    make_option(c("-u", "--contrast_numerator"), type="character", default=NULL, 
                                                    help="contrast numerator"),                                                        
                    make_option(c("-l", "--contrast_denominator"), type="character", default=NULL, 
                                                    help="contrast denominator"),
                    make_option(c("-q", "--padj_cutoff"), type="double", default=0.05, 
                                                    help="cutoff to use for adjusted p-value"),
                    make_option(c("-e", "--proj_id"), type="character", default=NULL,
                                                    help = "sample to keep"),                               
                    make_option(c("-f", "--dont_collapse"), action="store_true", default=FALSE,
                                                    help = "whether to collapse replicates"),  
                    make_option(c("-c", "--cell_type"), type="character", default="all", 
                                                    help="cell types to include in analysis"),
                    make_option(c("-b", "--brain_region"), type="character", default="all", 
                                                    help="brain regions to include in analysis"),
                    make_option(c("-s", "--sex"), type="integer", default=2, 
                                                    help="cell types to include in analysis (0 for female, 1 for male, 2 for both"),   
                    make_option(c("-g", "--case_control"), type="integer", default=2, 
                                                    help="whether to include cases or controls or both (0 for controls, 1 for cases, 2 for both"),
                    make_option(c("-n", "--nsc_cutoff"), type="double", default=1.0, 
                                                    help="cutoff to use for NSC"),
                    make_option(c("-r", "--rsc_cutoff"), type="double", default=0.4, 
                                                    help="cutoff to use for RSC"),
                    make_option(c("-a", "--pbc_cutoff"), type="double", default=0.5, 
                                                    help="cutoff to use for PBC"),
                    make_option(c("-o", "--out_prefix"), type="character", default=NULL, 
                                                    help="output prefix for 3 output files"),    
                    make_option(c("-t", "--transform"), action="store_true", default=FALSE, 
                                                    help="perform rlog transformation and plot pca")
                                                  
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print(opt)

if (is.null(opt$matrix_counts)){
  print_help(opt_parser)
  stop("At least one argument must be supplied ", call.=FALSE)
}


countsFile <- opt$matrix_counts
sampleInfoFile <- opt$info_sample
peakIDsFile <- opt$peak_info
designFormula <- opt$design
contrastVariable <- opt$variable_contrast
contrastNumerator <- opt$contrast_numerator
contrastDenominator <- opt$contrast_denominator
padjCutoff <- opt$padj_cutoff
sampleToKeep <- opt$proj_id
dontCollapse <- opt$dont_collapse
cellType <- opt$cell_type
brainRegion <- opt$brain_region
sex <- opt$sex
caseControl <- opt$case_control
nscCutoff <- opt$nsc_cutoff
rscCutoff <- opt$rsc_cutoff
pbcCutoff <- opt$pbc_cutoff
outPrefix <- opt$out_prefix
performRLog <- opt$transform

countMatrix <- read.table(countsFile, header = FALSE, sep = "\t", skip = 0)
sampleInfo <- read.table(sampleInfoFile, header = TRUE, row.names=1, sep = "\t", skip=0)
peakInfo <- read.table(peakIDsFile, header = FALSE, row.names=1, sep = "\t", skip=0)

#making sure the row names are valid variable names, essential for row names and column names to match
rownames(sampleInfo) <- make.names(rownames(sampleInfo))
all(rownames(sampleInfo) == colnames(countMatrix))
colnames(countMatrix) <- rownames(sampleInfo)
all(rownames(sampleInfo) == colnames(countMatrix))

sampleInfo$proj_id <- factor(sampleInfo$proj_id)
sampleInfo$cell_type <- factor(sampleInfo$cell_type)
sampleInfo$brain_region <- factor(sampleInfo$brain_region)
sampleInfo$attempt <- factor(sampleInfo$attempt)
sampleInfo$replicate <- factor(sampleInfo$replicate)
sampleInfo$binary_amyloid <- factor(sampleInfo$binary_amyloid)
sampleInfo$msex <- factor(sampleInfo$msex)
sampleInfo$is_microglia <- factor(as.integer(sampleInfo$cell_type == "Microglia"))
sampleInfo$is_glia <- factor(as.integer(sampleInfo$cell_type == "Glia"))
sampleInfo$is_neuron <- factor(as.integer(sampleInfo$cell_type == "Neuron"))
sampleInfo$amyloid_sqrt <- sqrt(sampleInfo$amyloid)
sampleInfo$apoe <- factor(sampleInfo$apoe)

#Filtering samples based on input criteria
if (sex!=2){
    sampleInfo <- sampleInfo[sampleInfo$msex==sex,] 
}

if (caseControl!=2){
    sampleInfo <- sampleInfo[sampleInfo$binary_amyloid==caseControl,]
}


sampleInfo <- sampleInfo[sampleInfo$NSC>nscCutoff,]
sampleInfo <- sampleInfo[sampleInfo$RSC>rscCutoff,]
sampleInfo <- sampleInfo[sampleInfo$PBC1>pbcCutoff,]
sampleInfo <- sampleInfo[sampleInfo$PBC2>pbcCutoff,]


if (cellType=="n"){
    sampleInfo <- sampleInfo[sampleInfo$is_neuron==1,]
} else if(cellType=="g"){
    sampleInfo <- sampleInfo[sampleInfo$is_glia==1,]
} else if(cellType=="m"){
    sampleInfo <- sampleInfo[sampleInfo$is_microglia==1,]
} else if(cellType=="nn"){
    sampleInfo <- sampleInfo[sampleInfo$is_neuron==0,]
}

if (brainRegion=="h"){
    sampleInfo <- sampleInfo[sampleInfo$brain_region=="HPC",]
} else if (brainRegion=="d"){
    sampleInfo <- sampleInfo[sampleInfo$brain_region=="DLPFC",]
}

if (!is.null(sampleToKeep)){
    sampleInfo <- sampleInfo[sampleInfo$proj_id==sampleToKeep,]
}

#Subsampling count matrix based on selected samples
countMatrix <- countMatrix[,rownames(sampleInfo)]
all(rownames(sampleInfo) == colnames(countMatrix))



library(DESeq2)
#initializing deseq matrix

print("Number of remaining samples")
print(nrow(sampleInfo))
print(sampleInfo$cell_type)
print(sampleInfo$proj_id)
print(sampleInfo$msex)
print(sampleInfo$binary_amyloid)
print(designFormula)
dds <-  DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = sampleInfo,
                              design = as.formula(paste("~", designFormula, sep=" ")))

#collapsing technical replicates                             
ddsCollapsed <- collapseReplicates( dds,
                                   groupby = make.names(paste(dds$proj_id,dds$brain_region, dds$cell_type,sep="_")))

                                   
# fitting deseq model
if(dontCollapse){
    dds <- DESeq(dds)
} else {                                   
    dds <- DESeq(ddsCollapsed)
}

print("Available results:")
print(resultsNames(dds))


if(!is.null(contrastVariable)){
    if(is.null(contrastNumerator)){
        res <- results(dds, alpha=padjCutoff, name=contrastVariable)
    } else {
        res <- results(dds, alpha=padjCutoff, contrast=c(contrastVariable,contrastNumerator,contrastDenominator))
    }
} else {
    useContrast <- TRUE
    cat("Results by name (y/n)?:")
    userPrompt <- readLines(file("stdin"), n = 1, ok = FALSE)
    if (userPrompt %in% c("y", "Y", "yes", "YES")) {
        useContrast <- FALSE
    }

    if (useContrast) {
    
        cat("Enter contrast name:")
        contrastVariable <- readLines(file("stdin"), n = 1, ok = FALSE)
        cat("Enter contrast numerator:")
        contrastNumerator <- readLines(file("stdin"), n = 1, ok = FALSE)
        cat("Enter contrast denominator:")
        contrastDenominator <- readLines(file("stdin"), n = 1, ok = FALSE)
        #extracting results
        res <- results(dds, alpha=padjCutoff, contrast=c(contrastVariable,contrastNumerator,contrastDenominator))
    } else {
        cat("Enter result name:")
        contrastVariable <- readLines(file("stdin"), n = 1, ok = FALSE)
        res <- results(dds, alpha=padjCutoff, name=contrastVariable)
    }
}


#printing results summary
print(summary(res))

#plotting MA plot
jpeg(paste(outPrefix,"MA.jpg",sep="_"))
plotMA(res, alpha=padjCutoff)
dev.off()

#plotting histogram of p-values
jpeg(paste(outPrefix,"pvalhist.jpg",sep="_"))
hist( res$pvalue, breaks=20, col="grey" )
dev.off()

rownames(res) <- rownames(peakInfo)
res_sorted <- res[order(res$padj),]
res_sorted_remove_na <- na.omit(res_sorted)


peaks_up <- res_sorted_remove_na[res_sorted_remove_na$log2FoldChange>0 & res_sorted_remove_na$padj < padjCutoff,]
peaks_down <- res_sorted_remove_na[res_sorted_remove_na$log2FoldChange<0 & res_sorted_remove_na$padj < padjCutoff,]


write.table(peakInfo[rownames(peaks_up),], file=paste(outPrefix,"up.txt",sep="_"), sep="\t",col.names=FALSE, quote=FALSE)
write.table(peakInfo[rownames(peaks_down),], file=paste(outPrefix,"down.txt",sep="_"), sep="\t",col.names=FALSE, quote=FALSE)
write.table(res_sorted_remove_na, file=paste(outPrefix,"results_sorted.txt",sep="_"), sep="\t", quote=FALSE)


# perform rlog transform and plot PCA if specified


if (performRLog) {
    print("Performing rlog transformation")
    rld <- rlog(dds)
    print("Finished rlog transformation")
    
    sampleDists <- dist(t(assay(rld)))

    sampleDistMatrix <- as.matrix( sampleDists )
    rownames(sampleDistMatrix) <- paste( rld$cell_type, rld$proj_id, sep="_" )
    colnames(sampleDistMatrix) <- paste( rld$cell_type, rld$proj_id, sep="_" )
    library( "gplots" )
    library( "RColorBrewer" )
    colours = colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

    jpeg(paste(outPrefix,"dendrogram.jpg",sep="_"))    
    heatmap.2( sampleDistMatrix,
              trace="none",
              col=colours,
              cexRow=0.5,
              cexCol=0.5,
              margins=c(8,8))
    dev.off()
    
    jpeg(paste(outPrefix,"pca.jpg",sep="_"))
    plotPCA(rld, intgroup = c("cell_type"))
    dev.off()
}

warnings()



