import sys
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import scanpy as sc

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42

samples = np.array(('DOW8726783', 'DOW8726784', 'DOW8726785', 'DOW8726786', 'DOW8726787', 'DOW8726788', 'DOW8726789', 'DOW8726790'))

for i in range(len(samples)):
    print(i)
    output_dir = '/home/jovyan/AA_DOW/scrublet/'
    input_dir = '/home/jovyan/data/AA_DOW/AA_DOW/cellranger302_count_33315_AA_' + samples[i] + '_GRCh38-3_0_0-premrna/'
    adata = sc.read_10x_h5('/home/jovyan/data/AA_DOW/AA_DOW/cellranger302_count_33315_AA_' + samples[i] + '_GRCh38-3_0_0-premrna/output_filtered.h5', genome='background_removed')
    counts_matrix = adata.X
    genes = np.array(scr.load_genes(input_dir + 'raw_feature_bc_matrix/features.tsv', delimiter='\t', column=1))

    print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
    print('Number of genes in gene list: {}'.format(len(genes)))

    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                              min_cells=3, 
                                                              min_gene_variability_pctl=85, 
                                                              n_prin_comps=30)

    predicted_doublets = predicted_doublets*1
    predicted_doublets = predicted_doublets.astype(int)
    detected_doublets_rate = round(scrub.detected_doublet_rate_, 4)
    overall_doublets_rate = round(scrub.overall_doublet_rate_, 4)

    np.savetxt(output_dir + '/' + samples[i] + '_' + 'doublets_scores.txt', doublet_scores)   
    np.savetxt(output_dir +  '/' + samples[i] + '_' + 'predicted_doublets.txt', predicted_doublets)                              
    with open(output_dir +  '/' + samples[i] + '_' + 'detected_doublets_rate.txt', 'w') as f:
      f.write('%f' % detected_doublets_rate)  

    with open(output_dir + '/' + samples[i] + '_' +  'overall_doublets_rate.txt', 'w') as f:
      f.write('%f' % overall_doublets_rate)

    f = scrub.plot_histogram()
    f[0].savefig(output_dir + '/' + samples[i] + '_' + "doubletScore_histogram.pdf", bbox_inches='tight')
    
