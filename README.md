# Single Cell Non-negative Matrix Factorization (scNMF) toolkit

scNMF is a toolkit for compressing single cell datasets, factorizing these compressed spaces, and visualizing the results.

# Divisive clustering
`FindClusters.Divisive`
The master function.

Divisive bipartite spectral clustering is a method for compressing information in single cell datasets.

# Non-negative matrix factorization
`RunNMF`

Currently only dense matrices are supported. Sparse matrices may be supported in the future to avoid transferring a dense matrix into memory, but NMF will not be faster since zeros are not treated as missing values.

# Visualization
Tools in development.

# Random thoughts
- Speed of factorization for sparse matrices has been dramatically increased by regarding zeros as missing values. However, this necessitates a calculation of loss only on positive values in the input matrix, and thus a very different result is obtained. In single cell data, zeros are not equivalent to missing values and therefore they must be weighted in the iterative NNLS steps.
