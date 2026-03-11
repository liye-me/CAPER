# CAPER
![CAPER](/figure/CAPER.png)
CAPER is primarily built on a matrix factorization framework that explicitly disentangles 
meaningful biological effects from unwanted variation. CAPER directly generates 
interpretable latent factors and a batch-corrected expression matrix in which 
the biological signal of interest is preserved and isolated, thereby enabling more 
reliable downstream analyses in complex multi-condition single-cell datasets.

# Installation
CAPER is implemented in C++, compiled using Rcpp, and distributed as an R package.

```
library(devtools)
install_github("liye-me/CAPER")
```

# Usage
The main function is CAPER_rcpp. You can find the instructions and an example by '?CAPER_rcpp'.

## Example
We apply CAPER to 4 scRNA-seq datasets from the peripheral blood mononuclear cells (PBMCs) of 2 patients with systemic lupus erythematosus (SLE), which included 2 IFN-β-stimulated samples (stim; N = 1,109 cells) and 2 unstimulated controls (control; N = 1,055 cells).

### Preprocessing Pipeline
```
library(Seurat)
library(sva)
object <- readRDS("/path/to/your/CAPER/inst/extdata/exampledata.rds")
objectlist <- SplitObject(object,split.by = "condition")
objectlist <- lapply(X = objectlist, FUN = function(x) {
x <- NormalizeData(x)
x <- ScaleData(x)
x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 500)
})
object_c <- objectlist[[1]]
object_s <- objectlist[[2]]
hvg.c <- head(VariableFeatures(object_c), n = 500)
hvg.s <- head(VariableFeatures(object_s), n = 500)
hvg.union <- union(x = hvg.c, y = hvg.s)
X_c <- object_c@assays$RNA@scale.data[hvg.union,]
X_s <- object_s@assays$RNA@scale.data[hvg.union,]
X_c <- ComBat(dat=as.matrix(X_c), batch=as.character(object_c$batch), mod=NULL, par.prior=TRUE, prior.plots=FALSE)
X_s <- ComBat(dat=as.matrix(X_s), batch=as.character(object_s$batch), mod=NULL, par.prior=TRUE, prior.plots=FALSE)
```

### CAPER Model
##### Inputs
- X_c: A gene-by-cell matrix for condition c. (G × n^c) → as.matrix(X_c)
- X_s: A gene-by-cell matrix for condition s. (G × n^s) → as.matrix(X_s)
- k1: The number of global latent dimensions.
- k2: The number of condition-specific latent dimensions.
- sigma_c: The variance scale for condition c. (mean of per-gene variances)
- sigma_s: The variance scale for condition s. (mean of per-gene variances)
- max_iter: The maximum number of EM iterations to perform.
- epsilon: The convergence threshold for the EM algorithm.
- verbose: If TRUE, progress information will be printed during execution (default is FALSE).

##### Outputs
A list of model outputs:
- Lambda: Cell-by-latent loading matrix Λ ((n^c+n^s) × (k1+2*k2))
- E_Z: Posterior mean E[Z|X] ((k1+2*k2) × G), where each column corresponds to a gene g and Z_g = [z_g; z_g^c; z_g^s]
- psi_diag: Diagonal of Ψ (noise covariance), length (n^c+n^s)
```
sigma_c <- mean(apply(X_c, 1, var))
sigma_s <- mean(apply(X_s, 1, var))
library(CAPER)
res <- CAPER_rcpp(X_c = as.matrix(X_c), X_s = as.matrix(X_s), k1 = 30, k2 = 10, sigma_c = sigma_c, sigma_s = sigma_s, max_iter = 10000, epsilon = 0.001, verbose = FALSE)
```
### Gene Expression Matrix Reconstruction and Loading Matrix Extraction
```
k1 <- 30

# Reconstruct (cells × genes)
Xhat <- res$Lambda[, 1:k1, drop = FALSE] %*% res$E_Z[1:k1, , drop = FALSE]
rownames(Xhat) <- c(colnames(X_c), colnames(X_s))
colnames(Xhat) <- hvg.union

# Create Seurat object (genes × cells)
Xhat_t <- t(Xhat)
correct <- CreateSeuratObject(counts = Xhat_t)
correct <- SetAssayData(correct, slot = "data", new.data = Xhat_t)
correct <- SetAssayData(correct, slot = "scale.data", new.data = Xhat_t)

correct <- FindVariableFeatures(correct)
correct$condition <- c(as.character(object_c$condition), as.character(object_s$condition))
correct$batch     <- c(as.character(object_c$batch),     as.character(object_s$batch))
correct$celltype  <- c(as.character(object_c$celltype),  as.character(object_s$celltype))

X_emb_matrix <- res$Lambda[, 1:k1, drop = FALSE]
rownames(X_emb_matrix) <- c(colnames(X_c), colnames(X_s))

correct[["emb"]] <- CreateDimReducObject(
  embeddings = X_emb_matrix,
  key = "X_",
  assay = DefaultAssay(correct)
)
```

## Our group

 <https://sqsun.github.io/>.