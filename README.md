# Analysis scripts for H3K27ac data from sorted nuclei of hippocampus and prefrontal cortex of Alzheimer's human postmortem brains

This is a release repository for code accompanying the paper titled "Cell type-specific histone acetylation profiling of Alzheimer's Disease and integration with genetics" 


## Plotting sample info jitter plots and heatmap

TODO

## Individual peak calling for all dlPFC and hippocampus samples using the AQUAS ChIP-seq workflow

The following commands will process hippocampus samples individually (starting from fastq stage)
```
sbatch hip_cases_individual_peak_calling.sb
sbatch hip_controls_individual_peak_calling.sb
```

The following command will process dlPFC samples individually (starting from bwa aligned bam stage - hg19)
```
sbatch dlpfc_chipseq_script.sb
```

The above commands will generate 78 peak calls corresponding to peaks active in each individual sample - (16 hippocampus samples + 10 dlPFC samples) x (3 cell type populations). Since, ChIPs and Input experiments were conducted in duplicate, the two ChIP and Input replicates were together input to the pipeline to generate these peak calls. These commands also compute quality control metrics such as NSC, RSC, FRiP for each individual sample. A table of these QC metrics is provided in the main paper for each sample.


## Calling reproducible peaks in each brain region separately for each cell type and subjects with and without Amyloid pathology

```
sbatch run_dlpfc_cases_by_cell_type.sb
sbatch run_dlpfc_controls_by_cell_type.sb
sbatch run_hip_cases_by_cell_type.sb
sbatch hip_controls_by_cell_type.sb
```

The above commands will generate 12 peak calls corresponding to peaks active in each cell type for the two different brain regions and subjects with and without amyloid pathology - (3 cell types x 2 brain regions x 2 amyloid/no amyloid)

## Merging and combining peaks to generate a consensus peak set

```
zcat /hpc/peaks/controls/Neuron/peak/macs2/overlap/optimal_set/Neuron_ppr.naive_overlap.filt.narrowPeak.gz /hpc/peaks/controls/Microglia/peak/macs2/overlap/optimal_set/Microglia_ppr.naive_overlap.filt.narrowPeak.gz /hpc/peaks/controls/Glia/peak/macs2/overlap/optimal_set/Glia_ppr.naive_overlap.filt.narrowPeak.gz /hpc/peaks/cases/Neuron/peak/macs2/overlap/optimal_set/Neuron_ppr.naive_overlap.filt.narrowPeak.gz /hpc/peaks/cases/Microglia/peak/macs2/overlap/optimal_set/Microglia_ppr.naive_overlap.filt.narrowPeak.gz /hpc/peaks/cases/Glia/peak/macs2/overlap/optimal_set/Glia_ppr.naive_overlap.filt.narrowPeak.gz /dlpfc/peaks/controls/Neuron/peak/macs2/overlap/optimal_set/Neuron_ppr.naive_overlap.filt.narrowPeak.gz /dlpfc/peaks/controls/Microglia/peak/macs2/overlap/optimal_set/Microglia_ppr.naive_overlap.filt.narrowPeak.gz /dlpfc/peaks/controls/Glia/peak/macs2/overlap/optimal_set/Glia_ppr.naive_overlap.filt.narrowPeak.gz /dlpfc/peaks/cases/Neuron/peak/macs2/overlap/optimal_set/Neuron_ppr.naive_overlap.filt.narrowPeak.gz /dlpfc/peaks/cases/Microglia/peak/macs2/overlap/optimal_set/Microglia_ppr.naive_overlap.filt.narrowPeak.gz /dlpfc/peaks/cases/Glia/peak/macs2/overlap/optimal_set/Glia_ppr.naive_overlap.filt.narrowPeak.gz | awk -vFS='\t' -vOFS='\t' '{print $1,$2,$3}' | sort -k1,1 -k2,2n | bedtools merge -i stdin -d 200 > dlpfc_hpc_combined_set_combined_200.bed
```

## Naming and numbering peaks

Function to name and number peaks. Takes an input prefix for the peak names and adds numbers automatically. Ordering of output is same as input.

```
name_peaks.py [-h] -i PEAK_FILE -o OUT -p PREFIX -w PADDING
```

Argument definitions:
```
  -h, --help            show this help message and exit
  -i PEAK_FILE, --peak-file PEAK_FILE
                        input peak bed file
  -o OUT, --out OUT     output peak bed file
  -p PREFIX, --prefix PREFIX
                        prefix of the label
  -w PADDING, --padding PADDING
                        length of the number for label
```

To name the consensus peak set, I use the following command:

```
python name_peaks.py -i dlpfc_hpc_combined_set_combined_200.bed -o dlpfc_hpc_combined_set_combined_200_name.bed -p LHTPFCHPC -w 15
```

To create a saf file (used for feature counting) from the above generated bed file, I use the following command:

```
awk -vFS='\t' -vOFS='\t' '{print $4,$1,$2+1,$3,"+"}' dlpfc_hpc_combined_set_combined_200_name.bed > dlpfc_hpc_combined_set_combined_200_name.saf
```

## Feature counts to compute read counts for every sample in the consensus peak set

```
sbatch feature_counts_dlpfc_control_samples.sb 
sbatch feature_counts_dlpfc_case_samples.sb 
sbatch feature_counts_hpc_control_samples.sb 
sbatch feature_counts_hpc_case_samples.sb 
```

To create a combined matrix with all the counts, I use the following command:

```
paste /counts/hpc_cases/*.countSimp.txt /counts/hpc_controls/*.countSimp.txt /counts/dlpfc_cases/*.countSimp.txt /counts/dlpfc_controls/*.countSimp.txt | tail -n +2 > /counts/all_samples.countSimp.no_header.txt
```

## DESeq2 wrapper

To easily call DESeq2 from the command line for running multiple DESeq2 experiments, I created a command line wrapper around DESeq2 that uses optparse.

```
Usage: run_deseq_all_samples.R [options]


Options:
	-m MATRIX_COUNTS, --matrix_counts=MATRIX_COUNTS
		file containing count matrix

	-i INFO_SAMPLE, --info_sample=INFO_SAMPLE
		file containing sample info

	-p PEAK_INFO, --peak_info=PEAK_INFO
		file containing peak info

	-d DESIGN, --design=DESIGN
		design string for analysis

	-v VARIABLE_CONTRAST, --variable_contrast=VARIABLE_CONTRAST
		output contrast variable

	-u CONTRAST_NUMERATOR, --contrast_numerator=CONTRAST_NUMERATOR
		contrast numerator

	-l CONTRAST_DENOMINATOR, --contrast_denominator=CONTRAST_DENOMINATOR
		contrast denominator

	-q PADJ_CUTOFF, --padj_cutoff=PADJ_CUTOFF
		cutoff to use for adjusted p-value

	-e PROJ_ID, --proj_id=PROJ_ID
		sample to keep

	-f, --dont_collapse
		whether to collapse replicates

	-c CELL_TYPE, --cell_type=CELL_TYPE
		cell types to include in analysis

	-b BRAIN_REGION, --brain_region=BRAIN_REGION
		brain regions to include in analysis

	-s SEX, --sex=SEX
		cell types to include in analysis (0 for female, 1 for male, 2 for both

	-g CASE_CONTROL, --case_control=CASE_CONTROL
		whether to include cases or controls or both (0 for controls, 1 for cases, 2 for both

	-n NSC_CUTOFF, --nsc_cutoff=NSC_CUTOFF
		cutoff to use for NSC

	-r RSC_CUTOFF, --rsc_cutoff=RSC_CUTOFF
		cutoff to use for RSC

	-a PBC_CUTOFF, --pbc_cutoff=PBC_CUTOFF
		cutoff to use for PBC

	-o OUT_PREFIX, --out_prefix=OUT_PREFIX
		output prefix for 3 output files

	-t, --transform
		perform rlog transformation and plot pca

	-h, --help
		Show this help message and exit
```

## Cell type-specific hyperacetylated peaks and Abeta associated DAR calling

```
source run_deseq_all_samples.sh
```
The above shell script calls the DESeq2 wrapper script 48 different times with different arguments. In the first three calls, I call the wrapper to contrast cell types against each other to identify cell type-specific peaks for Neuron, Microglia, and Glia. In the next 45 commands, I identify amyloid beta correlated DARs in different brain regions, cell types and sexes. These 45 contrasts correspond to the contrasts defined in Figure 3 and Supplementary Table 4 of the main manuscript text.

## Calling cell type reproducible peaks
I use the following script and command to call cell type reproducible peaks across both profiled brain regions by feeding in pooled alignments from hippocampus controls and dlPFC controls and assessing reproducibility between the two brain regions using the AQUAS ChIP-seq workflow.

```
source cell_type_reproducible_dlpfc_hpc_combined.sh
```

## S-LDSC analysis
Create merged background set for cell type-specific differential peaks
```
cat /deseq_analysis/cell_type_deseq/allbr_cell_type_effect_for_neuron_up.bed /deseq_analysis/cell_type_deseq/allbr_cell_type_effect_for_microglia_up.bed /deseq_analysis/cell_type_deseq/allbr_cell_type_effect_for_glia_up.bed | sort -k1,1 -k2,2n > bedtools merge -i stdin > merged_set_of_differential_cell_type_specific_peaks.bed
```

Running S-LDSC analysis
```
source run_ldsc_cell_type_specific.sh
```

## S-LDSC analysis (cell type reproducible peaks)
```
zcat /peaks/cell_type_peaks_hpc_dlpfc_reproducible/dlpfc_hpc_combined_controls_Neuron/peak/macs2/overlap/optimal_set/dlpfc_hpc_combined_controls_Neuron_ppr.naive_overlap.filt.narrowPeak.gz /peaks/cell_type_peaks_hpc_dlpfc_reproducible/dlpfc_hpc_combined_controls_Microglia/peak/macs2/overlap/optimal_set/dlpfc_hpc_combined_controls_Microglia_ppr.naive_overlap.filt.narrowPeak.gz /peaks/cell_type_peaks_hpc_dlpfc_reproducible/dlpfc_hpc_combined_controls_Glia/peak/macs2/overlap/optimal_set/dlpfc_hpc_combined_controls_Glia_ppr.naive_overlap.filt.narrowPeak.gz | sort -k1,1 -k2,2n | bedtools merge -i stdin > dlpfc_hpc_combined_controls_cell_type_reproducible.bed

```

```
source run_ldsc_cell_type_reproducible.sh
```

## Permutation test analysis (is the code for this ok to release???)

TODO

## Plotting GWAS enrichment results

TODO

## HOMER annotation
```
annotatePeaks.pl dlpfc_hpc_combined_set_combined_200_name.bed hg19 > dlpfc_hpc_combined_set_combined_200_name_homer_gene_annot.txt
annotatePeaks.pl /deseq_analysis/binary_amyloid_deseq/glia_dlpfc_binary_amyloid_up.bed hg19 > /deseq_analysis/binary_amyloid_deseq/glia_dlpfc_binary_amyloid_up_homer_gene_annot.txt
annotatePeaks.pl /deseq_analysis/sex_specific_amyloid_deseq/glia_hpc_sex_0_specific_binary_amyloid_down.bed hg19 > /deseq_analysis/sex_specific_amyloid_deseq/glia_hpc_sex_0_specific_binary_amyloid_down_homer_gene_annot.txt
```

```
awk -vFS='\t' -vOFS='\t' '{print $2,$3,$4,$1,$10,$16}' dlpfc_hpc_combined_set_combined_200_name_homer_gene_annot.txt > dlpfc_hpc_combined_set_combined_200_name_homer_gene_annot_disttss_gene_name_column.txt
```

## Analysis of Habib et al marker genes to validate the sorting of neuronal, microglial and glial populations

Script to extract peak names for peaks annotated to a given gene list:
```
Usage: get_peaks_matching_gene_ids.py [peak_annot_file] [gene_list_file] [output_file] [filter_disttss]

peak_annot_file - tab delimited file with columns representing chr, start, stop, peak name, distance to TSS, Nearest Gene Name
gene_list_file - list of gene names, one name per line
output_file - list of peak names that are output, one name per line
filter_disttss - Optional. Distance cutoff for filtering out peaks. Peaks with greater than filter_disttss distance from nearest gene will not be included in the output even if the gene is present in the input gene list. Default value: infinity to include all peaks annotated to genes in input set.
```

The following shell script will call the above python script to extract the peak names for peaks annotated to the 15 Habib et al cell type marker gene sets. Then, it will extract the cell type-specificity H3K27ac log2fc values from the 3 cell type contrasts in DESeq2 for each of the 15 "marker peak" sets.
```
source habib_markers_extract_l2fcs.sh
```

The following python script reads in the previously extracted cell type-specificity log2fc values and conducts a t-test to check for H3K27ac enrichment at the 15 "marker peak" sets in the 3 cell type populations we profiled. It then plots a heatmap with the t-test effect size along with p-values. It also performs multiple hypothesis correction using the Bonferroni adjustment across the 45 different tests (3 FANS populations x 15 Habib et al cell type clusters). This plot is presented in Supplementary Figure 6b in the main paper.

```
python habib_markers_plot_combined_l2fc.py
```

For restricting the analysis to peaks <5kb from nearest TSS, I use the following commands:

```
source habib_markers_disttss_filter_extract_l2fcs.sh
python habib_markers_disttss_filter_plot_combined_l2fc.py
```

These will generate the plot in Supplementary Figure 6a in the main paper.

## Habib markers individual sample analysis

TODO

## Plotting peak specificity at individual AD risk loci

TODO

## Plotting distance to TSS of Abeta associated DARs

TODO

## VST generation for full peak set for all samples

TODO

## VST generation for female hippocampus OEG, male hippocampus OEG, and dlPFC OEG samples

TODO

## Post-hoc analysis, box plots of covariates

TODO

## GO results plotting

TODO

## snRNA-seq Mathys et al analysis

TODO

## Averaging bigwigs

TODO

## pyGenomeTracks visualizations for Abeta associated DARs

TODO

# Misc

Below is a general purpose script for processing peak data. It uses argparse and takes in command line arguments.


## Calculating Average Peak Length

Function to calculate average peak length in a narrowPeak bed file. Takes an input narrowPeak file (or possibly other bed format files) and prints the average region length.

```
get_average_peak_length.py [-h] -i PEAK_FILE
```

Argument definitions:
```
  -h, --help            show this help message and exit
  -i PEAK_FILE, --peak-file PEAK_FILE
                        Input file containing peaks
```
