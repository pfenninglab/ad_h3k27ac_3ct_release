import matplotlib
matplotlib.use('agg')
matplotlib.rcParams['svg.fonttype'] = 'none'
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from scipy import stats, misc
import numpy as np


def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero.

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower offset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax / (vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highest point in the colormap's range.
          Defaults to 1.0 (no upper offset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap


h3k27ac_cell_type_names = ["neuron", "microglia", "glia"]
marker_names = ["NeuN+", "Pu.1+", "NeuN-/Pu.1-"]
brain_regions = ["allbr"]
clusters = ["exPFC1", "exPFC2", "exCA1", "exCA3", "GABA1", "GABA2", "exDG", "MG", "ODC1", "ODC2", "OPC", "ASC1", "ASC2", "NSC", "END"]
numClusters = len(clusters)
totGreaterThanZeroTests = len(clusters)*len(h3k27ac_cell_type_names)*len(brain_regions)
totCellTypeTests = misc.comb(len(h3k27ac_cell_type_names),2)*len(clusters)*len(brain_regions)
padj_threshold = 0.05
padj_stringent = 0.01
logfc_cutoff = 0.5
numColsPlot = 5
numRowsPlot = 3
h1 = 0.2
h2 = 0.05
prefix = "/habib_markers_analysis/"
suffix = "l2fc.txt"



tStatMatrix = dict()
pValMatrix = dict()
l2fcMeansMatrix = dict()

plotColor = dict()

nameMapping = {
                "exPFC1": "Excitatory PFC 1",
                "exPFC2": "Excitatory PFC 2",
                "exCA1": "Excitatory HPC CA1",
                "exCA3": "Excitatory HPC CA3",
                "GABA1": "GABAergic 1",
                "GABA2": "GABAergic 2",
                "exDG": "Excitatory Dentate Gyrus",
                "MG": "Microglia",
                "ODC1": "Oligodendrocyte 1",
                "ODC2": "Oligodendrocyte 2",
                "OPC": "Olig Precursor",
                "ASC1": "Astrocyte 1",
                "ASC2": "Astrocyte 2",
                "NSC": "Neural Stem Cell",
                "END": "Endothelial"
              }

for cluster in ["exPFC1", "exPFC2", "exCA1", "exCA3", "GABA1", "GABA2", "exDG"]:
    plotColor[cluster] = (0.65,0.16,0.16,0.3)
    
for cluster in ["MG"]:
    plotColor[cluster] = (0,0,1,0.3)
    
for cluster in ["ODC1", "ODC2"]:
    plotColor[cluster] = (1,1,0,0.3)
    
for cluster in ["OPC", "ASC1", "ASC2", "NSC", "END"]:
    plotColor[cluster] = (0,0.5,0,0.3)

for brain_region in brain_regions:
    tStatMatrix[brain_region] = []
    pValMatrix[brain_region] = []
    l2fcMeansMatrix[brain_region] = []
    plt.figure(figsize=(20, 12))
    for i, cluster in enumerate(clusters):
        allL2fcData=[]
        t_stats = []
        p_vals = []
        for h3k27ac_cell_type_name in h3k27ac_cell_type_names:
            l2fcFile = prefix + "_".join([brain_region, cluster, h3k27ac_cell_type_name, suffix])
            l2fcData = []
            with open(l2fcFile, 'r') as f:
                for line in f:
                    l2fc = float(line.strip())
                    l2fcData.append(l2fc)
            allL2fcData.append(l2fcData)

            # test whether mean is close to 0
        t_stats = []
        p_vals = []
        l2fcMeans = []
        for x in range(len(allL2fcData)):
            t_stat, p_val = stats.ttest_1samp(allL2fcData[x], logfc_cutoff)
            l2fcMean = np.mean(allL2fcData[x])
            #one-tailed p value (divide by 2)
            p_val = p_val/2
            t_stats.append(t_stat)
            p_vals.append(p_val)
            l2fcMeans.append(l2fcMean)
        
        t_stats = np.array(t_stats)
        p_vals = np.array(p_vals)
        l2fcMeans = np.array(l2fcMeans)
        
        l2fcMeansMatrix[brain_region].append(l2fcMeans)
        tStatMatrix[brain_region].append(t_stats)
        
        ax = plt.subplot(numRowsPlot, numColsPlot, i+1, facecolor=plotColor[cluster])
        ax.set_title(nameMapping[cluster], fontsize=18)
        ax.set_ylabel('log2FC (cell type)', fontsize=15)
        ax.tick_params(axis='x', labelsize=15)
        ax.boxplot(allL2fcData, labels=marker_names)
        greater_zero_adj_p_vals = p_vals*totGreaterThanZeroTests
        pValMatrix[brain_region].append(greater_zero_adj_p_vals)
        for j, adj_p_val in enumerate(greater_zero_adj_p_vals):
            if adj_p_val <= padj_threshold:
                x1 = j-h1+1
                x2 = j-h1-h2+1
                y1 = logfc_cutoff
                y2 = np.median(allL2fcData[j])
                if y2>y1:
                    ax.plot([x1,x2,x2,x1], [y1,y1,y2,y2], c='k')
                    sigString = "*"
                    if adj_p_val <= padj_stringent:
                        sigString = "**"
                    ax.text(x2-h2,(y1+y2)/2.0, sigString, ha='right', va='center', color='k')
        
        # test each cell type against the other
        for j in range(len(allL2fcData)):
            for k in range(j+1, len(allL2fcData)):
                if (np.median(allL2fcData[j])>0 and greater_zero_adj_p_vals[j]<padj_threshold) and (np.median(allL2fcData[k])>0 and greater_zero_adj_p_vals[k]<padj_threshold):
                        t_stat, p_val = stats.ttest_ind(allL2fcData[j], allL2fcData[k], equal_var=False)
                        adj_p_val = p_val*totCellTypeTests
                        if adj_p_val <= padj_threshold:
                            x1 = j+1+h2
                            x2 = k+1-h2
                            y1 = np.amax([np.amax(val) for val in allL2fcData])+(k-j)*h1
                            y2 = y1+h2
                            ax.plot([x1,x1,x2,x2], [y1,y2,y2,y1], c='k')
                            sigString = "*"
                            if adj_p_val <= padj_stringent:
                                sigString = "**"
                            ax.text((x1+x2)/2.0, y2, sigString, ha='center', va='bottom', color='k')
        
        
        ax.axhline(y=0.0, linestyle='dashed', c='k')
        ax.set_ylim([-3,3])
    
    plt.tight_layout()
    plt.savefig(prefix + "_".join([brain_region, "l2fchist.svg"]))
    plt.close()

for brain_region in tStatMatrix:
	print(brain_region)
	print(tStatMatrix[brain_region])
	print(pValMatrix[brain_region])


for brain_region in l2fcMeansMatrix:
    currMatrix = np.stack(l2fcMeansMatrix[brain_region]).T
    print(currMatrix.shape)
    currPVals = np.stack(pValMatrix[brain_region]).T
    plt.figure(figsize=(20,12))
    fig, ax = plt.subplots()
    
    vmax = np.amax(currMatrix)
    vmin = np.amin(currMatrix)
    vmid = 1 - vmax / (vmax + abs(vmin))
    new_cmap = shiftedColorMap(matplotlib.cm.bwr, midpoint=vmid, name='shrunk')
    
    
    im = ax.imshow(currMatrix, cmap=new_cmap)

    ax.set_xticks(np.arange(currMatrix.shape[1]))
    ax.set_yticks(np.arange(currMatrix.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels([nameMapping[cluster] for cluster in clusters])
    ax.set_yticklabels(marker_names)
                   

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=False, bottom=False,
                    left=False, right=False,
                    labeltop=True, labelbottom=False)    
    
    
    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=90, ha="left",
             rotation_mode="anchor")

    for i in range(currPVals.shape[0]):
        for j in range(currPVals.shape[1]):
            val = np.round(-np.log10(min(currPVals[i,j],1)),decimals=0)
            if not val==float("inf"):
                val = int(val)
            if currMatrix[i,j] > 0:
                text = ax.text(j, i, val,
                    ha="center", va="center", color="k", fontsize=12)
         

    cbar = plt.colorbar(im, orientation="horizontal", shrink=0.50)
    cbar.ax.xaxis.set_ticks_position('top')
    cbar.ax.xaxis.set_label_position('bottom')
    cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation=45, ha='left')
    cbar.ax.set_xlabel("mean(log2FC)", va="top")
    
    plt.savefig(prefix + "_".join([brain_region,"meanl2fc","heatmap",".svg"]))
    plt.close()
    
for brain_region in tStatMatrix:
    currMatrix = np.stack(tStatMatrix[brain_region]).T
    currPVals = np.stack(pValMatrix[brain_region]).T
    plt.figure(figsize=(20,12))
    fig, ax = plt.subplots()
    
    vmax = np.amax(currMatrix)
    vmin = np.amin(currMatrix)
    vmid = 1 - vmax / (vmax + abs(vmin))
    new_cmap = shiftedColorMap(matplotlib.cm.bwr, midpoint=vmid, name='shrunk')
    
    
    im = ax.imshow(currMatrix, cmap=new_cmap)

    ax.set_xticks(np.arange(currMatrix.shape[1]))
    ax.set_yticks(np.arange(currMatrix.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels([nameMapping[cluster] for cluster in clusters])
    ax.set_yticklabels(marker_names)
                   

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=False, bottom=False,
                    left=False, right=False,
                    labeltop=True, labelbottom=False)    
    
    
    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=90, ha="left",
             rotation_mode="anchor")

    for i in range(currPVals.shape[0]):
        for j in range(currPVals.shape[1]):
            val = np.round(-np.log10(min(currPVals[i,j],1)),decimals=0)
            if not val==float("inf"):
                val = int(val)
            if currMatrix[i,j] > 0:
                text = ax.text(j, i, val,
                    ha="center", va="center", color="k", fontsize=12)
         

    cbar = plt.colorbar(im, orientation="horizontal", shrink=0.50)
    cbar.ax.xaxis.set_ticks_position('top')
    cbar.ax.xaxis.set_label_position('bottom')
    cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation=45, ha='left')
    cbar.ax.set_xlabel("t-statistic", va="top")
    
    plt.savefig(prefix + "_".join([brain_region,"tstat","heatmap",str(logfc_cutoff)+".svg"]))
    plt.close()
    print(brain_region)
    print(currPVals)
    print(l2fcMeansMatrix[brain_region])
