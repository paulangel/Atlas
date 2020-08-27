import numpy as np
import pandas as pd
import sklearn
import gc
import functions
import scipy

data           = pd.read_csv('/Users/pwangel/Downloads/myeloid_atlas_expression_v7.1.tsv', sep='\t', index_col=0)
previous_genes = pd.read_csv('/Users/pwangel/PlotlyWorkspace/combine_data/blood/outputs_for_front_end/iMac_genes.tsv', sep='\t', index_col=0)
annotations    = pd.read_csv('/Users/pwangel/PlotlyWorkspace/combine_data/blood/outputs_for_front_end/iMac_annotations.tsv', sep='\t', index_col=0)
annotations['Dataset'] = [i.split(';') for i in annotations.index.values.astype(str)]

annotations = annotations.loc[annotations.Dataset!='6638'] #### Do this if you want to get rid of dataset 6638
data        = data[annotations.index]

data = functions.transform_to_percentile(data)

genes = functions.calculate_platform_dependence(data, annotations)

functions.plot_atlas(data, genes, annoations, '/Users/pwangel/Downloads')


#KW_Htest_results = functions.plot_KW_Htest(data, annotations, genes, '/Users/pwangel/Downloads/')

#functions.plot_gene_platform_dependence(data, annotations, genes, '/Users/pwangel/Downloads/')

#H_index_list, retained_genes_list = functions.resample_clustering(data, annotations, resample_strategy='bootstrap', n_resamples=2, n_clusters_list=[3,4])

