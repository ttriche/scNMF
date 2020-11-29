# Single Cell Non-negative Matrix Factorization (scNMF) toolkit

scNMF is a toolkit for compression of single cell datasets (divisive clustering), fast factorization of these compressed spaces (NMF), and high-level visualization of the results.

Vignettes:
 - Divisive cell clustering (basic version, advanced version)
 - Dimensional reduction with scNMF (basic version, advanced version)
 - Programs of genetic coactivation with deep scNMF (basic version, advanced version)
 - Exploring geractions with deep scNMF (basic version, advanced version)

## 1. Divisive clustering
`scNMF::FindClusters.Divisive()` is a very fast divisive clustering method which operates on sparse single cell matrices in Seurat objects. It can be used for fast clustering as an alternative to existing graph-based approaches.

Divisive bipartite spectral clustering is a method for compressing information in single cell datasets.

## 2. Non-negative matrix factorization
`scNMF::nnmf()` is the fastest NMF algorithm we are aware of. It builds off off `NNLM::nnmf()`, but does not support non-random initialization, regularization, or masking, and uses only sequential coordinate descent to solve least squares against a mean squared error loss function.

### Software and Implementation

###
*Sparse matrix support:* Currently `nnmf` does not currently support sparse matrices. This is planned for future development. However, the only advantage of sparse support is that matrices will not need to be densely mapped to memory, there will be no speed gains in the NMF solver itself because zeros are not treated as missing values.

*Zeros not treated as missing values:* Algorithms which are faster than `scNMF::nnmf()` likely treat zeros as missing values (i.e. `rsparse::WRMF`). When zeros are treated as missing values, the model becomes imputation-greedy and seeks to maximize prediction of only positive values in thei nput. This ignores the fact that many values in single cell count matrices are truly zero, despite some amount of zero-inflation. Thus, while zeros may be treated as missing values in recommender systems (i.e. Netflix prize, movie recommendations), they should not be treated as missing values in single cell experiments.

*Cross-validation:* 

## 3. Visualization
Tools in development.
