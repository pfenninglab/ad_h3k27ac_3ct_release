import matplotlib
matplotlib.use('agg')
matplotlib.rcParams['svg.fonttype'] = 'none'
import matplotlib.pyplot as plt
import numpy as np


new_bar_width = 0.1
data = np.array([[1604, 627, 52231],
                 [277,273, 175865],
                 [81,129,120305],
                 [0,0,3611]])
data_normed = data / data.sum(axis=0)
colorLabels = ["#ff8c00ff", "#990033ff", 'white']
labels = ["Glia Female Hip. Hypoacetylated", "Glia dlPFC Hyperacetylated", "Full Peak Set"]
plt.figure()
for i in range(data.shape[1]):
    plt.bar(np.arange(data.shape[0])+i*new_bar_width, data_normed[:,i], new_bar_width, color=colorLabels[i], edgecolor='k', label=labels[i])
plt.xlabel('Absolute distance to TSS (kb)')
plt.ylabel('Peaks')
plt.xticks(np.arange(data.shape[0])+new_bar_width, ('0 to 5', '5 to 50', '50 to 500', '> 500'))
plt.yticks()
plt.legend(bbox_to_anchor=(0.5,0.0), loc='upper center')
plt.savefig('peaks_sets_disttss.svg')


