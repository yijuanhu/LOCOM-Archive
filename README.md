# LOCOM

A logistic regression model for testing differential abundance in compositional microbiome data

The LOCOM package provides the locom function for testing differential abundance of individual taxa (or OTUs) that is based on the log ratio of relative abundances between each taxon and a reference taxon. It also provides a global test of whether there is any differential abundance in a microbiome community. The tests accommodate continuous and discrete (binary) traits of interest and allows adjustment of confounders.

## Installation

```r
install.packages(c("Rcpp", "RcppArmadillo", "metap", "utils", "psych", "permute", "parallel", "devtools", "BiocManager"))
BiocManager::install("BiocParallel")
devtools::install_github("yijuanhu/LOCOM")
```

## An example

Apply LOCOM to a dataset from the study of smoking effects on the microbiome of the human upper respiratory tract. (Charlson et al., 2010). 

```r
# loading data
data("throat.otu.table")
data("throat.meta")
data("throat.otu.taxonomy")
otu.table <- data.matrix(throat.otu.table)
Y <- ifelse(throat.meta$SmokingStatus == "NonSmoker", 0, 1)
C <- data.matrix(model.matrix(Y ~ throat.meta$Sex + throat.meta$AntibioticUsePast3Months_TimeFromAntibioticUsage - 1))[, -1]

# filtering out three samples with antibiotic use
filter.out.sam <- which(C[, 3] == 0)
otu.table.filter <- otu.table[-filter.out.sam,]
Y <- Y[-filter.out.sam]
C <- C[-filter.out.sam,]

# filtering out rare OTUs
prop.presence.thresh <- 0.2
prop.presence <- colMeans(otu.table.filter > 0)
filter.out.otu <- which(prop.presence < prop.presence.thresh)
if (length(filter.out.otu) > 0) {
    otu.table.filter <- otu.table.filter[, -filter.out.otu]
    prop.presence <- prop.presence[-filter.out.otu]
}

# running locom
res <- locom(otu.table = otu.table.filter, Y = Y, C = C[, 1], fdr.nominal = 0.1, seed = 1)

# summarizing results
# ordering the detected OTUs by their p-values. If no OTU is detected, we can still provide a summary table for
# the top (e.g., 10) OTUs by re-defining o = order(res$p.otu)[1:10]
w <- match(res$detected.otu, colnames(res$p.otu))
o <- w[order(res$p.otu[w])]

summary.table <- data.frame(otu.name = colnames(res$p.otu)[o],
                            mean.freq = colMeans(otu.table.filter/rowSums(otu.table.filter))[o],
                            prop.presence = prop.presence[o],
                            p.value = signif(res$p.otu[o], 3),
                            q.value = signif(res$q.otu[o], 3),
                            effect.size = signif(res$effect.size[o], 3),
                            otu.tax = throat.otu.taxonomy[as.numeric(colnames(res$p.otu)[o]) + 1],
                            row.names = NULL)
summary.table
```

