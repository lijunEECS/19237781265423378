# File structure

*code_2Dtest* : 2D demo version. Run code in this folder to generate figures for 2D demostration.

*code* : code used to solve 360D simulation problem.

# HSCS code structure

*HSCS.m* : main program, consist of (1)presampling, (2)kmeans clustering and (3)mixture importance sampling.

## kmeans cluster (including calculate min-norm points) uses the following functions:

*kmeans_main.m* : main function to do k-means clustering on failed samples.

*Wkmeans.m* : to do weighted kmeans clustering.

*Cluster_norm.m* : to calculate the norm used in kmeans. (cosine distance in current version)

*DBI/DVI.m* : Davies-Bouldin Index/ Dunn Validity Index. Two indices used to evaluate cluster results.

*sampleWeight.m* : to calculate the weight of each sample based on pdf.

*reindex.m* : to map the indices to [1,M] after removing some empty clusters.

*min_norm.m* : to capture the min norm of a cluster center. (refer to HSCS paper for details)

*isFailure.m* : run SPICE simulation to check if a sample is failed or not.

## mixture importance sampling uses the following functions:

*generateMISSamples.m* : generate a batch of samples based on cluster results.

*weight.m* : calculate the likelihood ratio.

*mis_main.m* : main function to do mixture importance sampling based on cluster results.

*getMeanSigma_cp22nmdemo.m* : generate mean and sigma values.
