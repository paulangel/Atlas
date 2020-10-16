import numpy as np
import pandas as pd
import sklearn, sys, gc, scipy
import functions

blood_atlas_colours = pd.read_csv('/Users/pwangel/Atlas/data/blood_atlas_colours.tsv', sep='\t').set_index('Sample Source')
blood_atlas_colours = {key:value[0] for key, value in zip(blood_atlas_colours.index.values, blood_atlas_colours.values)}

yidis_data           = pd.read_csv('/Users/pwangel/Downloads/imac_probit.tsv', sep='\t', index_col=0)
yidis_data.columns   = [i.replace('.', ';') for i in yidis_data.columns.values]
yidis_annotations   = pd.read_csv('/Users/pwangel/Downloads/yidis_annotations.tsv', sep='\t', index_col=0)
yidis_annotations.index  = [i.replace('.', ';') for i in yidis_annotations.index.values]

#data           = pd.read_csv('/Users/pwangel/Downloads/myeloid_atlas_expression_v7.1.tsv', sep='\t', index_col=0)
#annotations    = pd.read_csv('/Users/pwangel/PlotlyWorkspace/combine_data/blood/outputs_for_front_end/i_annotations.tsv', sep='\t', index_col=0)
data           = pd.read_csv('/Users/pwangel/PlotlyWorkspace/combine_data/blood/outputs_for_front_end/blood_atlas_expression.tsv', sep='\t', index_col=0)
annotations    = pd.read_csv('/Users/pwangel/PlotlyWorkspace/combine_data/blood/outputs_for_front_end/blood_atlas_annotations.tsv', sep='\t', index_col=0)

previous_genes = pd.read_csv('/Users/pwangel/PlotlyWorkspace/combine_data/blood/outputs_for_front_end/iMac_genes.tsv', sep='\t', index_col=0)
previous_genes.index = [i.replace('-', '_') for i in previous_genes.index.values]
previous_genes.index = [i.replace('.', '_') for i in previous_genes.index.values]

#yidis_genes    = pd.read_csv('/Users/pwangel/Downloads/variance partition result.txt', sep='\t', index_col=0)
yidis_genes    = pd.read_csv('/Users/pwangel/Downloads/yidi_var_part_noweights.tsv', sep='\t', index_col=0)
yidis_genes.index = [i.replace('-', '_') for i in yidis_genes.index.values]
yidis_genes.index  = [i.replace('.', '_') for i in yidis_genes.index.values]

data        = data[annotations.index]
yidis_data  = yidis_data[yidis_annotations.index]

old_genes   = yidis_data.index.values
yidis_data.index  = [i.replace('.', '_') for i in yidis_data.index.values]
yidis_data.index  = [i.replace('-', '_') for i in yidis_data.index.values]

print(stop)

#data = functions.transform_to_percentile(data)

sys.path.append('/Users/pwangel/Gene_Analysis/combine_data')
import mega_functions

print("Calculating platform dependence")
python_platform_varPart = functions.calculate_platform_dependence(yidis_data, yidis_annotations)
R_platform_varPart        = mega_functions.cut_genes_that_depend_on_platform(yidis_data, yidis_annotations, 0.2, False)

import matplotlib
import matplotlib.pyplot as pyplot

pyplot.scatter(python_platform_varPart.Platform_VarFraction, yidis_genes.Batch.values, s=10) 
pyplot.xlabel("R Variance Partition Platform Only")
pyplot.ylabel("R Yidis Platform and Celltype")

pyplot.show()

pyplot.scatter(python_platform_varPart.Platform_VarFraction, celltype_varPart.Platform_VarFraction, s=10) 
pyplot.xlabel("R Variance Partition Platform Only")
pyplot.ylabel("Python Platform and Celltype")

pyplot.show()

pyplot.scatter(yidis_genes.Batch.values, celltype_varPart.Platform_VarFraction, s=10) 
pyplot.xlabel("R Platform and Celltype")
pyplot.ylabel("Python Platform and Celltype")

pyplot.show()

print(stop)

celltype_varPart = functions.calculate_celltype_dependence(yidis_data, yidis_annotations)
celltype_varPart.to_csv('/Users/pwangel/Downloads/python_celltype_platform_probit.tsv', sep='\t')

#functions.plot_pca(data, annotations, genes, labels=['Dataset', 'celltype', 'Platform_Category'], colour_dict=blood_atlas_colours)

#KW_Htest_results = functions.plot_KW_Htest(data, annotations, genes, '/Users/pwangel/Downloads/')

#functions.plot_gene_platform_dependence(data, annotations, genes, '/Users/pwangel/Downloads/')

#H_index_list, retained_genes_list = functions.resample_clustering(data, annotations, resample_strategy='bootstrap', n_resamples=2, n_clusters_list=[3,4])

