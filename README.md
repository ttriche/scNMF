# scNMF: single cell non-negative matrix factorization toolkit

scNMF is a toolkit for compression of single cell datasets (divisive clustering) and fast factorization of these compressed spaces (NMF). NMF on single cell data can learn the architecture of gene coactivation programs that yield observed transcriptional states, and can be used to visualize cell clusters based on these coactivated gene programs.

## Vignettes
1. Divisive cell clustering (basic version, advanced version)
2. Dimensional reduction with NMF for cell visualization (basic version, advanced version)
3. Deep NMF to learn programs of differential genetic coactivation (basic version, advanced version)
// Later? Exploring genetic interactions with systematic deep scNMF (basic version, advanced version)

## Main parts of publication:
Part 1. Fast and robust NMF on single-cell datasets after compression by divisive clustering
 - A cross-validation objective for robustness: bipartite matching between factorizations of non-overlapping datasets
 - Robustness of factorization of raw Seurat example datasets with various NMF algorithms (NNLM::nnmf, scNMF, CoGAPS, scater RunNMF)
 - Robustness of factorization of Seurat datasets compressed by divisive clustering
 
Part 2. NMF captures cellular identity better than PCA
 - Cells can be projected onto a lower-dimensional manifold from their NMF factor mappings. NMF coordinates resolve true cell identities better than PCA (Seurat CITE-seq dataset)
 - Divisive clustering resolves cellular identities better than PCA graph-based Louvain clustering
 - Divisive cellular clustering tree from RNA-seq data with overlaid Seurat CITE-seq protein expression
 - Louvain cellular clustering tree (hierarchical clustering of determined clusters) with overlaid Seurat CITE-seq protein expression

Part 3. Deep NMF captures genetic coactivation programs
 - Framework for deep NMF
 - Deep NMF on Mouse Organogenesis Cell Atlas
 - Deep NMF on Human Cell Landscape
 - Gene set enrichment in NMF factors (heatmap)
 - LDAvis of selected factors
 - Hierarchical clustering tree of an individual gene (i.e. FZD4) showing all factors in which it expresses

Functions yet to write:
 - Unrooted igraph network plot of cell clustering graph














 *To do:*
 Vignette 1: Cluster cells in the Seurat hcabm40k dataset as basic introduction, visualize on UMAP coordinates. Visualize cells on UMAP coordinates. Visualize cells with igraph network (write visualization tool). Demonstrate use of different stopping critera for divisive cell clustering (i.e. Newman-Girvan modularity, cosine distance, min cells only)
 Vignette 1 advanced: Cluster cells in Seurat pbmc CITE-seq data to show that clustering captures cell clusters marked by protein genes better than Louvain clustering
 
Vignette 2: Run NMF on a clustering from vignette 1 on hcabm40k, average expression, run NMF with function that maps dimensional reduction to all cells and saves the slot in Seurat object. Use NMF to project UMAP/tSNE coordinates and compare to PCA coordinates. Compare cellular distributions of clustering identities determined by divisive clustering vs. graph-based clustering. UMAP plots make clear that NMF resolves cellular populations better, transient populations linked by an "ithsmus" between clusters.

Vignette 3: Run divisive clustering on hcabm40k, run NMF, run NMF2 on NMF$W, run NMF3 on NMF$W2 examine Gene Ontology term enrichment in W1, W2, and W3. Use LDAVis to visualize genes in various factors. Measure enrichment of learned factors across cell types, use Seurat to visualize enrichment patterns.

Vignette 4 - maybe move to a later publication? Run divisive clustering on multiple datasets (MOCA, gastrulation, MCA, Tabula Muris, Tabula Muris Senis), combine to run NMF, save data from second-generation NMF in package as data("mouse_gene_programs"). Annotate with paracrine coupling. Explore WNT and FZD matching, coupling across sets. Filter by membrane-integral (receptors), nuclear-intrinsic (transcription factors), secreted, annotate by autocrine or paracrine. Examine signaling pathway structures.

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
