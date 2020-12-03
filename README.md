# scNMF: single cell non-negative matrix factorization toolkit

scNMF is a toolkit for compression of single cell datasets (divisive clustering) and fast factorization of these compressed spaces (NMF). NMF on single cell data can learn the architecture of gene coactivation programs that yield observed transcriptional states, and can be used to visualize cells on UMAP/tSNE coordinates based on these coactivated gene programs.

## Vignettes
1. Divisive cell clustering 
2. Dimensional reduction with NMF for cell visualization
3. Deep NMF to learn programs of differential genetic coactivation 
4. Benchmarking scNMF against existing NMF methods (figure 1) 
// Later? Exploring genetic interactions with systematic deep scNMF

# Manuscript concept

### Title: Deep NMF of Single Cell Transcriptomes Learns Gene Coactivation Programs

#### Figure 1. Fast and robust NMF on single-cell datasets after compression by divisive clustering (hcabm40k dataset)
<ol>
 <li>Schematic of genes x cells -> compressed genes x cell clusters -> NMF -> A) Deep NMF for gene set learning or B) cell visualization</li>
 <li>Number of clusters generated from a range of distance-based stopping criteria, modularity-based stopping criteria, or min cells stopping criteria</li>
 <li>Robustness of various NMF methods on raw data. CoGAPS/scNMF/scater/nnmf</li>
 <li>Robustness of various NMF methods on compressed data. CoGAPS/scNMF/scater/nnmf</li>
 <li>Runtime for each method</li>
</ol>

#### Figure 2. Visualization of cell clusters on NMF coordinates
<ol>
 <li>Schematic: Projection of cells onto a lower-dimensional manifold from NMF factor mappings with NNLS</li>
  <li>tSNE plot of cells from Figure 1 on PCA coordinates with Louvain clustering (Seurat DimPlot)</li>
  <li>tSNE plot of cells from Figure 1 on PCA coordinates with divisive clustering (Seurat DimPlot)</li>
  <li>tSNE plot of cells from Figure 1 on NMF coordinates with Louvain clustering (Seurat DimPlot)</li>
  <li>tSNE plot of cells from Figure 1 on NMF coordinates with divisive clustering (Seurat DimPlot)</li>
</ol>

#### Figure 3. Deep NMF captures genetic coactivation programs
<ol>
 <li>Schematic: Framework for deep NMF</li>
 <li>Deep NMF on Mouse Organogenesis Cell Atlas at each time point</li>
 <li>Gene ontology term enrichment in deep NMF factors (fgsea, gheatmap)</li>
 <li>fgsea enrichment plot for selected factor (fgsea)</li>
 <li>LDAvis of specific factors (LDAvis)</li>
 <li>Plot factors on tSNE plot with NMF coordinates (i.e. FZD4/Lrp5/Tspan12), cell clusters as assigned by MOCA authors, labeled with MOCA labels</li>
</ol>

#### Figure 4. WNT pathway transcriptional architecture as learned by Deep NMF
<ol>
 <li>Unrooted graph of hierarchical clustering of FZD+ vectors showing FZD1-10 bias as pie charts.</li>
 <li>Cell type enrichment of each vector</li>
 <li>Transmembrane gene associations with FZDs in each factor (GPCRs, LRPs, kinases, other)</li>
 <li>Transcription factors marking each factor</li>
 <li>Lrp1/Lrp4/Lrp5/Lrp6/Lrp8/Vldlr/Ldlr coexpression (cosine similarity) bar chart with each FZD</li>
 <li>Co-incidence of Ras/Mapk/Egfr/Erk with Wnt, shown as a bar chart</li>
</ol>

#### Data availability:
<ul>
 <li>Compressions of the MOCA dataset available as a sparse matrix</li>
 <li>hcabm40k dataset is available through the `SeuratData` package</li>
 <li>Cell annotations for the MOCA dataset are available as `cell_annotate.csv` on the MOCA figshare page</li>
 <li>The R package is publicly available on GitHub and will be submitted to Bioconductor at the time of submission</li>
</ul>

<hr>
<br>
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
