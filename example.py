import numpy as np
import pandas as pd
import sklearn
import gc
import functions

blood_data        = pd.read_csv('/Users/pwangel/PlotlyWorkspace/combine_data/blood/outputs_for_front_end/blood_atlas_expression.tsv', sep='\t', index_col=0)
previous_genes    = pd.read_csv('/Users/pwangel/PlotlyWorkspace/combine_data/blood/outputs_for_front_end/blood_atlas_genes.tsv', sep='\t', index_col=0)
blood_annotations = pd.read_csv('/Users/pwangel/PlotlyWorkspace/combine_data/blood/outputs_for_front_end/blood_atlas_annotations.tsv', sep='\t', index_col=0)

#genes = functions.calculate_platform_dependence(blood_data, blood_annotations)

H_index_list, retained_genes_list = functions.resample_clustering(blood_data, blood_annotations, resample_strategy='bootstrap', n_resamples=2, n_clusters_list=[3,4])

