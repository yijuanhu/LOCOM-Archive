#' A logistic regression model for testing differential abundance in compositional microbiome data (LOCOM)
#'
#' This function allows you to test
#' (1). whether any OTU (or taxon) is associated with the trait of interest with FDR control, based on the log ratio of relative abundances between pairs of taxa, and
#' (2). whehter the whole community is associated with the trait (a global test).
#' The tests accommodate both continuous and discrete (binary) traits and allows adjustment of confounders.
#'
#' This function uses a sequential stopping criterion (similar to that of Sandve et al. 2011) for the permutation procedure,
#' which stops when all taxon-level tests have either reached the pre-specified
#' number (default 100) of rejections or yielded a q-value (by the Benjamini-Hochberg [BH] procedure) that is below the
#' nominal FDR level (default 0.2). The permutation procedure is always terminated if a pre-specified maximum number (see description of \code{n.perm.max}) of
#' permutations have been generated. The global test is based on all permutation replicates when the procedure stops/terminates.
#'
#'
#' @param otu.table the OTU table (or taxa count table) in which the rows correspond to samples and the columns correspond to OTUs (taxa).
#' @param Y the trait of interest.
#' @param C the other (confounding) covariates to be adjusted.
#' @param fdr.nominal the nominal FDR. The default is 0.2.
#' @param seed a user-supplied integer seed for the random number generator in the
#'   permutation procedure. The default is NULL, which means that an integer seed will be
#'   generated internally and randomly. In either case, the integer seed will be stored
#'   in the output object in case the user wants to reproduce the permutation replicates.
#' @param prev.cut a real value between 0 and 1; taxa with prevalence (percentage of nonzeros) less than prev.cut are excluded. The default is 0.2.
#' @param n.perm.max the maximum number of permutations. The default is NULL, in which case \code{n.otu} * \code{n.rej.stop} * (1/\code{fdr.nominal})
#'   is used where \code{n.otu} is the total number of OTUs (that have non-zero counts in at least one sample).
#' @param n.rej.stop the minimum number of rejections (i.e., the permutation test
#'   statistic exceeds the observed test statistic) to obtain before stopping the permutation procedure. The
#'   default is 100.
#' @param n.cores the number of cores to be used for parallel computing. The default is 1.
#' @param adjustment method to adjust p-value: Benjamini-Hochberg procedure (BH) or Sandve's adjustment 
#' @param ref.otu reference OTU. The default is NULL, which means the most abundant OTU will be chosen as the reference OTU to fit the model.
#' @return A list consisting of
#' \itemize{
#'   \item effect.size - effect size at each OTU, i.e., beta_j,1 - median_j'=1,...J(beta_j',1)
#'   \item p.otu - p-values for OTU-specific tests
#'   \item q.otu - q-values (adjusted p-values by the BH procedure) for OTU-specific tests
#'   \item detected.otu - detected OTUs (using the column names of the OTU table) at the nominal FDR
#'   \item p.global - p-value for the global test
#'   \item n.perm.completed - number of permutations completed
#'   \item seed - the seed used to generate the permutation replicates
#' }
#' @export
#' @examples
#' # loading data
#'
#' data("throat.otu.table")
#' data("throat.meta")
#' data("throat.otu.taxonomy")
#'
#' otu.table <- data.matrix(throat.otu.table)
#'
#' Y <- ifelse(throat.meta$SmokingStatus == "NonSmoker", 0, 1)
#' C <- data.matrix(model.matrix(Y ~ throat.meta$Sex + throat.meta$AntibioticUsePast3Months_TimeFromAntibioticUsage - 1))[, -1]
#'
#' # filtering out three samples with antibiotic use
#'
#' filter.out.sam <- which(C[, 3] == 0)
#' otu.table.filter <- otu.table[-filter.out.sam,]
#' Y <- Y[-filter.out.sam]
#' C <- C[-filter.out.sam,]
#'
#' # filtering out rare OTUs
#'
#' prop.presence.thresh <- 0.2
#' prop.presence <- colMeans(otu.table.filter > 0)
#' filter.out.otu <- which(prop.presence < prop.presence.thresh)
#' if (length(filter.out.otu) > 0) {
#'     otu.table.filter <- otu.table.filter[, -filter.out.otu]
#'     prop.presence <- prop.presence[-filter.out.otu]
#' }
#'
#' # running locom
#' 
#' res <- locom(otu.table = otu.table.filter, Y = Y, C = C[, 1], fdr.nominal = 0.1, seed = 1, adjustment = "Sandev", n.cores = 4)
#'
#' # summarizing results
#' # ordering the detected OTUs by their p-values. If no OTU is detected, we can still provide a summary table for
#' # the top (e.g., 10) OTUs by re-defining o = order(res$p.otu)[1:10]
#'
#' w <- match(res$detected.otu, colnames(res$p.otu))
#' o <- w[order(res$p.otu[w])]
#'
#' summary.table <- data.frame(otu.name = colnames(res$p.otu)[o],
#'                   mean.freq = colMeans(otu.table.filter/rowSums(otu.table.filter))[o],
#'                   prop.presence = prop.presence[o],
#'                   p.value = signif(res$p.otu[o], 3),
#'                   q.value = signif(res$q.otu[o], 3),
#'                   effect.size = signif(res$effect.size[o], 3),
#'                   otu.tax = throat.otu.taxonomy[as.numeric(colnames(res$p.otu)[o]) + 1],
#'                   row.names = NULL)
#' summary.table
locom <- function(otu.table, Y, C = NULL,
                  fdr.nominal = 0.2, seed = NULL, prev.cut = 0.2,
                  n.perm.max = NULL, n.rej.stop = 100, n.cores = 4, ref.otu = NULL,
                  adjustment = "BH"){
    
    # remove zero OTUs
    
    w = which(colSums(otu.table>0)>= prev.cut * nrow(otu.table))
    if (length(w) < ncol(otu.table)) {
        cat(paste(ncol(otu.table)-length(w), 'OTU(s) with fewer than', prev.cut * nrow(otu.table), 'in all samples are removed', sep=" "), "\n")
        otu.table = otu.table[,w,drop=FALSE]
    }
    
    n.sam <- nrow(otu.table)
    n.otu <- ncol(otu.table)
    
    # find reference OTU
    
    if (is.null(ref.otu)) {
        mean.freq <- colMeans(otu.table/rowSums(otu.table))
        ref.otu <- which.max(mean.freq)
    }
    otu.table[which(otu.table[,ref.otu]==0), ref.otu] <- 1 # replace 0 by 1
    
    # -----------------------
    # Observed statistic
    # -----------------------
    
    if (is.null(C) == TRUE) {
        X <- cbind(Y, 1)
        Yr <- Y
    } else {
        X <- cbind(Y, 1, C)
        Yr <- resid(lm(Y ~ C))
    }
    
    freq.table.ref <- otu.table/(otu.table + otu.table[,ref.otu])
    n.X <- ncol(X)
    XX <- CalculateXX(X)
    weight = (otu.table + otu.table[,ref.otu])
    tol = 1e-6
    iter_max = 50
    Firth_thresh = 0.4
    prop.presence <- colMeans(otu.table>0)
    
    res.obs <- Newton(freq_table = freq.table.ref, X = X, XX = XX, 
                      beta_init = array(0, dim = c(n.X, n.otu)), weight = weight,
                      tol = tol, iter_max = iter_max, Firth_thresh = Firth_thresh,
                      prop_presence = prop.presence)
    
    beta <- res.obs[1,]
    beta <- beta - median(beta)
    
    shrinkage <- ifelse(is.null(C), 0.5, 0.75)
    beta_init = rbind(0, shrinkage*res.obs[2:n.X,])
    
    # -----------------------
    # Permutation
    # -----------------------
    
    if (is.null(seed)) seed = sample(1:10^6, 1)
    set.seed(seed)
    
    if (is.null(n.perm.max)) n.perm.max = n.otu * n.rej.stop * (1/fdr.nominal)
    
    beta.perm <- array(NA, dim = c(n.perm.max, n.otu))
    n.rej.otu.left <- rep(0, n.otu)
    n.rej.otu.right <- rep(0, n.otu)
    
    n.perm.block <- 1000
    n.block <- n.perm.max/n.perm.block
    n.perm.core <- n.perm.block/n.cores
    
    
    parallel.perm <- function(i) {
        
        perm_Newton(freq_table = freq.table.ref, Yr = Yr, X = X, XX = XX, 
                    beta_init = beta_init, weight = weight,
                    perm = perm.mat[, (i*n.perm.core + 1):((i+1)*n.perm.core)],
                    tol = tol, iter_max = iter_max, Firth_thresh = Firth_thresh, 
                    prop_presence = prop.presence)
    } # parallel.perm

    
    n.perm.completed <- 0
    
    for(i.block in 1:n.block){
        
        perm.mat <- t(shuffleSet(n.sam, n.perm.block)) - 1
        
        cat("permutations:", n.perm.completed + 1, "\n")
        
        if (n.cores > 1) {
            
            if (Sys.info()[['sysname']] == 'Windows') {
                parallel.stat = bplapply(0:(n.cores - 1), parallel.perm, BPPARAM = MulticoreParam(workers=n.cores))
            } else {
                parallel.stat = mclapply(0:(n.cores - 1), parallel.perm, mc.cores = n.cores)
            }
            res.perm <- do.call(rbind, parallel.stat)
            
        } else {
            
            res.perm <- perm_Newton(freq_table = freq.table.ref, Yr = Yr, X = X, XX = XX,
                                    beta_init = beta_init, weight = weight,
                                    perm = perm.mat,
                                    tol = tol, iter_max = iter_max, Firth_thresh = Firth_thresh,
                                    prop_presence = prop.presence)
        }
        
        n.perm.completed <- i.block*n.perm.block
        n.perm.completed.inv <- 1 / (n.perm.completed+1)
        
        w <- ((i.block-1)*n.perm.block+1):n.perm.completed
        beta.perm[w,] <- res.perm
        beta.perm[w,] <- beta.perm[w,] - apply(beta.perm[w,], 1, median)
        
        n.rej.otu.left <- n.rej.otu.left + rowSums(beta>=t(beta.perm[w,]))
        n.rej.otu.right <- n.rej.otu.right + rowSums(beta<=t(beta.perm[w,]))
        n.rej.otu <- 2*pmin(n.rej.otu.left+1, n.rej.otu.right+1)
        p.otu <- n.rej.otu * n.perm.completed.inv
        q.otu <- fdr.Sandev(p.otu)
        
        if (all(q.otu <= fdr.nominal | n.rej.otu >= n.rej.stop)) break
        
    } # permutation
    
    if(adjustment == "BH"){
        q.otu <- p.adjust(p.otu, method = "BH")
    }
    detected.otu <- colnames(otu.table)[which(q.otu < fdr.nominal)]
    
    # ------------------------
    # Global p-value
    # ------------------------
    
    beta.all <- rbind(beta, beta.perm)
    
    p.otu1 <- apply( beta.all, 2, function(x) 0.5*(2*pmin(rank(x), n.perm.completed +2-rank(x))*2-3)*n.perm.completed.inv )
    
    stat.global <- rowSums(1/p.otu1)
    p.global <- (sum(stat.global[1] <= stat.global[-1]) + 1) * n.perm.completed.inv
    
    
    otu.names <- colnames(otu.table)
    beta <- matrix(beta, nrow = 1)
    p.otu <- matrix(p.otu, nrow = 1)
    q.otu <- matrix(q.otu, nrow = 1)
    colnames(beta) <- otu.names
    colnames(p.otu) <- otu.names
    colnames(q.otu) <- otu.names
    
    
    return(list(effect.size = beta,
                p.otu = p.otu,
                q.otu = q.otu,
                detected.otu = detected.otu,
                p.global = p.global,
                n.perm.completed = n.perm.completed,
                seed = seed))
    
} # locom


fdr.Sandev = function(p.otu) {
    m = length(p.otu)
    p.otu.sort = sort(p.otu)
    n.otu.detected = seq(1, m)
    pi0 = min(1, 2/m*sum(p.otu))

    qval.sort = m * pi0 * p.otu.sort / n.otu.detected
    j.min.q = 1
    while (j.min.q < m) {
        min.q = min( qval.sort[j.min.q:m] )
        new.j.min.q = (j.min.q-1) + max( which(qval.sort[j.min.q:m]==min.q) )
        qval.sort[j.min.q:new.j.min.q] = qval.sort[new.j.min.q]
        j.min.q = new.j.min.q+1
    }
    mat = match(p.otu, p.otu.sort)
    qval.orig = qval.sort[mat]
    results = qval.orig

    return(results)

} # fdr.Sandev




