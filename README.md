# atlas

Simple code for generating transcriptomic atlases 

Consists of a collection of python functions in the functions.py file, and an example of how to use them, in the example.ipynb. 

Within the functions.py file there are the following functions:

transform_to_percentile: Transforms expression values in a dataframe to percentile values.

calculate_platform_dependence: Applies the variance modelling to determine the fraction of each genes variance due to the platform effect.

resample_clustering: Applies resampling to determine the stability of the atlas. Either bootstrap of leave-one-dataset-out resampling. Can take a long time as at each resampling iteration the variance modelling is redone.

calc_H_index: Using in the resampling_clustering function as the measure of cluster stability.

plot_KW_Htest: Scans through a range of platform variance threshold and at each threshold uses the Kruskal Wallis H test to measure the dependence of each principal component upon platform. 

KruskalWallisHTest: the Kruskal Wallis H test function.

plot_gene_platform_dependence_distribution: Plots the distribution of platform variance.

plot_gene_platform_dependence_cumulative_distribution: Plots the cumulative distribution of platform variance (i.e. how many genes pass a given threshold).

plot_pca: Plots the 3d PCA of the atlas.
