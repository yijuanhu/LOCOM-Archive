scenario <- "binary_Y_only"      # binary_Y_only, continuous_Y_only, binary_Y_binary_C, binary_Y_continuous_C, continuous_Y_binary_C, continuous_Y_continuous_C
causal.type <- 1   
n.sam <- 100        
beta <- 3         
have.bias <- FALSE
n.cores <- 1

n0.sam <- n.sam * 0.5
n1.sam <- n.sam * 0.5
C <- NULL
C.name <- "C"
disp <- 0.02        
a0 <- (1-disp)/disp
depth.mu <- 10000 
depth.lower <- 2000     
fdr.target <- 0.2 
filter.thresh <- 0.2 
n.sim <- 1

depth.sd <- depth.mu/3 

library(dirmult)
library(vegan)
library(Rcpp)
library(RcppArmadillo)
library(R.utils)
library(metap)
library(permute)
library(LOCOM)
source("ancom_v2.1.R")

library(ALDEx2)
library(ANCOMBC)
library(Wrench)
library(dacomp)
library(exactRankTests)
library(phyloseq)


pi <- read.table("input_throat/fit_dirmult_pi.txt", header=FALSE, as.is=TRUE)[,1]
n.otus <- length(pi)
ref.otu <- which.max(pi) # 468


if (scenario == "binary_Y_only" | scenario == "continuous_Y_only"){  # no confounder
    
    betaC <- NULL
  
    if (causal.type == 1){
        otu.Y.only <- c(which(pi >= 0.005))[1:20]
    }
    if (causal.type == 2){
        otu.Y.only <- order(pi, decreasing = TRUE)[1:5]
    }
    if (causal.type == 3){
        set.seed(0)
        otu.Y.only <- sort(sample(1:n.otus, 500))
    }
    if (causal.type == 4){
        otu.Y.only <- c(which(pi >= 0.001 & pi <= 0.002))[1:20]
    }
    otu.Y.C <- NULL
    otu.C.only <- NULL
    causal.otus <- c(otu.Y.only)
  
} else {
  
    betaC <- 2 # confounding effect
  
    if (causal.type == 1){
        w <- which(pi >= 0.002)
        otu.Y.only <- w[1:10]
        otu.C.only <- w[21:30]
        otu.Y.C <- w[11:20]
    } 
    if (causal.type == 2){
        o <- order(pi, decreasing = TRUE)
        otu.Y.only <- o[1:3]
        otu.C.only <- o[31:33]
        otu.Y.C <- o[4:5]
    }
    causal.otus <- c(otu.Y.only, otu.Y.C)
}
noncausal.otus <- setdiff(1:n.otus, causal.otus)

if (causal.type %in% c(3, 4)){
    set.seed(0)
    beta.otu <- runif(length(otu.Y.only), 1/beta, beta)
} else {
    beta.otu <- rep(beta, length(otu.Y.only))
}
beta.otu.log <- log(beta.otu)

if (have.bias == TRUE){
    set.seed(0)
    bias.factor.log <- rep(NA, n.otus)
    bias.factor.log[otu.Y.only] <- rnorm(length(otu.Y.only), 1, 0.8)
    bias.factor.log[-otu.Y.only] <- rnorm(n.otus-length(otu.Y.only), 0, 0.8)
    bias.factor <- exp(bias.factor.log)
}


summarize_otu_results <- function(qvalue, causal.otus, noncausal.otus, fdr.target=0.2) {
  
    otu.detected = sort(which(qvalue < fdr.target))
    n.otu = length(otu.detected)
  
    if (n.otu > 0) {
        sen = sum(otu.detected %in% causal.otus)/length(causal.otus)
        sep = 1 - sum(otu.detected %in% noncausal.otus)/length(noncausal.otus)
        fdr = n.otu - sum(otu.detected %in% causal.otus)
        fdr = fdr/n.otu
    } else {
        sen = 0
        sep = 1
        fdr = 0
    }
  
    out = list(n.otu=n.otu, sen=sen, sep=sep, fdr=fdr)
  
    return(out)
}


otu.all <- list()
p.global.all <- list()

for (i.seed in 1:n.sim) {
  
    cat("i.seed", i.seed, "\n")
  
    set.seed(i.seed)
  
    # -----------------------------------
    # Simulating data
    # -----------------------------------
  
    if (scenario == "binary_Y_only"){
        Y <- c(rep(0, n0.sam), rep(1, n1.sam))
    } else if (scenario == "continuous_Y_only"){
        Y <- runif(n.sam, -1, 1)
    } else if (scenario == "binary_Y_binary_C"){
        Y <- c(rep(0, n0.sam), rep(1, n1.sam))
        rho <- 0.2
        C <- c(rbinom(n0.sam, 1, rho), rbinom(n1.sam, 1, 1 - rho))
    } else if (scenario == "binary_Y_continuous_C"){
        Y <- c(rep(0, n0.sam), rep(1, n1.sam))
        C <- c(runif(n0.sam, -1, 1), runif(n1.sam, 0, 2))
    } else if (scenario == "continuous_Y_binary_C"){
        C <- c(rep(0, n0.sam), rep(1, n1.sam))
        Y <- c(runif(n0.sam, -1, 1), runif(n1.sam, 0, 2))
    } else if (scenario == "continuous_Y_continuous_C"){
        Y <- runif(n.sam, -1, 1)
        Y_indep <- runif(n.sam, -1, 1)
        rho <- 0.5
        C <- rho * Y + sqrt(1-rho^2) * Y_indep 
    }
  
    depth.sim <- rnorm(n.sam, depth.mu, depth.sd)
    depth.sim[depth.sim < depth.lower] <- depth.lower
    depth.sim <- round(depth.sim)
  
    otu.table.sim = rdirichlet(n = n.sam, a0*pi)
  
    # Introducing experimental bias
  
    if (have.bias == TRUE){
        otu.table.sim <- t(t(otu.table.sim) * bias.factor)
    }
  
    # Introducing the effects of Y
  
    otu.table.sim[, otu.Y.only] <- otu.table.sim[, otu.Y.only] * exp(Y %*% t(beta.otu.log))
  
    # Introducing the effects of C
  
    if (!is.null(betaC)){
        otu.table.sim[, otu.Y.C] <- otu.table.sim[, otu.Y.C] * beta^Y * betaC^C
        otu.table.sim[, otu.C.only] <- otu.table.sim[, otu.C.only] * betaC^C
    }
  
    # Normalizing
  
    otu.table.sim <- otu.table.sim/rowSums(otu.table.sim)
  
    # Generating count data
  
    for (i in 1:n.sam) {
        otu.table.sim[i,] <- rmultinom(1, depth.sim[i], otu.table.sim[i,])
    }
    colnames(otu.table.sim) <- c(1:n.otus)
  
  
    # -----------------------------------
    #  Filtering otus
    # -----------------------------------
  
    prop.presence <- colMeans(otu.table.sim > 0)
    otus.keep <- which(prop.presence >= filter.thresh)
    
    otu.table.sim.filter <- otu.table.sim[, otus.keep]
    ref.otu.filter <- which(colnames(otu.table.sim.filter) == ref.otu)
    causal.otus.filter = which(otus.keep %in% causal.otus)
    noncausal.otus.filter = setdiff(1:length(otus.keep), causal.otus.filter)
    
    cat(n.otus-length(otus.keep), "rare otus are filtered out", "\n")
    cat(length(causal.otus.filter), "causal otus remained after filtering\n") 
  
  
    ###########
    # LOCOM
    ###########
    
    res.locom <- locom(otu.table = otu.table.sim.filter, Y = Y, C = C, 
                       ref.otu = ref.otu.filter, fdr.nominal = fdr.target, 
                       n.perm.max = 50000, n.rej.stop = 100, n.cores = n.cores)	
    
    otu.locom <- summarize_otu_results(qvalue = res.locom$q.otu, causal.otus = causal.otus.filter, noncausal.otus = noncausal.otus.filter)
  
  
    #############
    # PERMANOVA
    #############
    
    otu.table.pc <- otu.table.sim.filter + 0.5  # pseudocount 0.5
    otu.table.clr <-  log(otu.table.pc) - apply(otu.table.pc, 1, function(x) mean(log(x)))
    otu.dist <- vegdist(otu.table.clr, method="euclidean")
    if (is.null(C)){
        res.adonis2 <- adonis2(otu.dist ~ Y, data = data.frame(Y=Y), permutations = 9999)
    } else{
        res.adonis2 <- adonis2(otu.dist ~ C + Y, data = data.frame(Y=Y, C=C), permutations = 9999)
    }
  
  
    ###########
    # ALDEx2
    ###########
    
    otu.table.aldex2 <- t(otu.table.sim.filter)
    if (scenario == "binary_Y_only"){
        aldex.sample <- aldex.clr(otu.table.aldex2, as.character(Y), mc.samples=128, verbose=TRUE)
        res.aldex2 <- aldex.ttest(aldex.sample, paired.test=FALSE, hist.plot = FALSE)
        aldex2.qvalue <- res.aldex2$we.eBH
    } else {
        if (is.null(C)){
            cov <- data.frame("Y" = Y)
            mm <- model.matrix(~ Y, cov)
        } else {
            cov <- data.frame("Y" = Y, C)
            colnames(cov)[2:ncol(cov)] <- C.name
            mm <- model.matrix( as.formula(paste(" ~ Y + ", paste(C.name, collapse = " + "), sep = "")), cov)
        }
        aldex.sample <- aldex.clr(otu.table.aldex2, mm, mc.samples=128, verbose=TRUE)
        res.aldex2 <- aldex.glm(aldex.sample)
        aldex2.qvalue <- res.aldex2$`model.Y Pr(>|t|).BH`
    }
    otu.aldex2 <- summarize_otu_results(qvalue = aldex2.qvalue, causal.otus = causal.otus.filter, noncausal.otus = noncausal.otus.filter)
    
    
    if (scenario == "binary_Y_only" | scenario == "binary_Y_binary_C" | scenario == "binary_Y_continuous_C"){
      
        ###########
        # ANCOM
        ###########
        
        otu.table.ancom <- t(otu.table.sim.filter)
        sam.ID <- paste0("ID", c(1:n.sam))
        colnames(otu.table.ancom) <- sam.ID
        
        if (is.null(C)){
            meta_data <- data.frame(sam.ID = sam.ID, Y = Y)
            adj_formula <- NULL
        } else {
            meta_data <- data.frame(sam.ID = sam.ID, Y = Y, C)    
            colnames(meta_data)[3:ncol(meta_data)] <- C.name
            adj_formula = paste(C.name, collapse = " + ")
        }
        
        prep <- feature_table_pre_process(feature_table = otu.table.ancom, meta_data = meta_data, sample_var = "sam.ID", 
                                          group_var = NULL, out_cut = 0.05, zero_cut = 1, lib_cut = 1000, neg_lb = FALSE)
        res.ancom <- ANCOM(prep$feature_table, prep$meta_data, prep$structure_zeros, main_var = "Y", p_adj_method = "BH", fdr.target, adj_formula = adj_formula, rand_formula = NULL)
      
        detected.otu.ancom <- which(res.ancom$out[, paste("detected_", 0.9, sep = "")]== TRUE)
        ancom.qvalue <- rep(1, length(otus.keep))
        ancom.qvalue[detected.otu.ancom] <- 0
        
        otu.ancom <- summarize_otu_results(qvalue = ancom.qvalue, causal.otus = causal.otus.filter, noncausal.otus = noncausal.otus.filter)
        
        ###########
        # ANCOM-BC
        ###########
        
        otu.table.ancombc <- t(otu.table.sim.filter)
        colnames(otu.table.ancombc) <- sam.ID
        
        if (is.null(C)){
            meta.data <- data.frame(group = Y, row.names = sam.ID, stringsAsFactors = FALSE)
            formula <- "group"
        } else {
            meta.data <- data.frame(group = Y, C, row.names = sam.ID, stringsAsFactors = FALSE)
            colnames(meta.data)[2:ncol(meta.data)] <- C.name
            formula <- paste("group + ",  paste(C.name, collapse = " + "), sep = "")
        }
        
        tax.table <- matrix(sample(letters, 7*length(otus.keep), replace = TRUE), nrow = length(otus.keep), ncol = 7)
        rownames(tax.table) <- rownames(otu.table.ancombc)
        colnames(tax.table) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
        
        otu.phylo <- otu_table(otu.table.ancombc, taxa_are_rows = TRUE)
        meta.phylo <- sample_data(meta.data)
        tax.phylo <- tax_table(tax.table)
        phylo.obj <- phyloseq(otu.phylo, meta.phylo, tax.phylo)
        
        res.ancombc <- ancombc(phyloseq = phylo.obj, formula = formula,
                               p_adj_method = "holm", prv_cut = 0, lib_cut = 1000, # 0.2: LOCOM filter; 0.1: ANCOM filter
                               group = "group", struc_zero = TRUE, neg_lb = FALSE,
                               tol = 1e-5, max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)
        
        otu.ancombc <- summarize_otu_results(qvalue = res.ancombc$res$q_val[,1], causal.otus = causal.otus.filter, noncausal.otus = noncausal.otus.filter)
    }
    
    if (scenario == "binary_Y_only" | scenario == "continuous_Y_only"){
      
        ###########
        # DACOMP
        ###########
        
        selected.ref.otus = dacomp.select_references(X = otu.table.sim.filter, median_SD_threshold = 1.3, verbose = F)
        test.method = ifelse(length(unique(Y)) == 2, DACOMP.TEST.NAME.WILCOXON, DACOMP.TEST.NAME.SPEARMAN)
        
        res.dacomp = dacomp.test(X = otu.table.sim.filter, y = Y, ind_reference_taxa = selected.ref.otus, 
                                 test = test.method, verbose = F, q = 0.05)
        dacomp.qvalue <- p.adjust(res.dacomp$p.values.test, method = 'BH')
        
        otu.dacomp <- summarize_otu_results(qvalue = dacomp.qvalue, causal.otus = causal.otus.filter, noncausal.otus = noncausal.otus.filter)
    }
    
    if (scenario == "binary_Y_only"){
      
        ###########
        # WRENCH
        ###########
        
        res.wrench <- wrench(t(otu.table.sim.filter), condition = Y, z.adj = TRUE)
        otu.table.wrench <- otu.table.sim.filter*res.wrench$nf
        wrench.pvalue <- apply(otu.table.wrench, 2, function(x) wilcox.exact(x[Y==1], x[Y==0])$p.value)
        wrench.qvalue <- p.adjust(wrench.pvalue, "BH")
        otu.wrench <- summarize_otu_results(qvalue = wrench.qvalue, causal.otus = causal.otus.filter, noncausal.otus = noncausal.otus.filter)
        
        #################
        # Pseudocount 0.5
        #################
        
        otu.table.pc <- otu.table.sim.filter + 0.5 
        otu.table.pc.alr <- log(otu.table.pc/otu.table.pc[,ref.otu.filter])
        pc.pvalue <- apply(otu.table.pc.alr, 2, function(x) wilcox.exact(x[Y==1], x[Y==0])$p.value)
        pc.qvalue <- p.adjust(pc.pvalue, "BH")
        otu.pseudo.half <- summarize_otu_results(qvalue = pc.qvalue, causal.otus = causal.otus.filter, noncausal.otus = noncausal.otus.filter)
        
        ###########
        # TSS
        ###########
        
        otu.table.tss <- otu.table.sim.filter/rowSums(otu.table.sim.filter)
        tss.pvalue <- apply(otu.table.tss, 2, function(x) wilcox.exact(x[Y==1], x[Y==0])$p.value)
        tss.qvalue <- p.adjust(tss.pvalue, "BH")
        otu.tss <- summarize_otu_results(qvalue = tss.qvalue, causal.otus = causal.otus.filter, noncausal.otus = noncausal.otus.filter)
    }
    
    
    # -----------------------------
    # Combining results
    # -----------------------------
    
    p.global.all[[i.seed]] <- data.frame("LOCOM" = res.locom$p.global, "PERMANOVA-half" = res.adonis2["Y", 5])
    
    if (scenario == "binary_Y_only") {
        otu.all[[i.seed]] <- rbind(unlist(otu.locom), unlist(otu.aldex2), unlist(otu.ancom), unlist(otu.ancombc), unlist(otu.dacomp), unlist(otu.wrench), unlist(otu.tss), unlist(otu.pseudo.half))
        rownames(otu.all[[i.seed]]) <- c("LOCOM", "ALDEx2", "ANCOM", "ANCOM-BC", "DACOMP", "WRENCH", "Wilcox-TSS", "Wilcox−alr−half")
    } else if (scenario == "binary_Y_binary_C" | scenario == "binary_Y_continuous_C") {
        otu.all[[i.seed]] <- rbind(unlist(otu.locom), unlist(otu.aldex2), unlist(otu.ancom), unlist(otu.ancombc))
        rownames(otu.all[[i.seed]]) <- c("LOCOM", "ALDEx2", "ANCOM", "ANCOM-BC")
    } else if (scenario == "continuous_Y_only") {
        otu.all[[i.seed]] <- rbind(unlist(otu.locom), unlist(otu.aldex2), unlist(otu.dacomp))
        rownames(otu.all[[i.seed]]) <- c("LOCOM", "ALDEx2", "DACOMP")
    } else if (scenario == "continuous_Y_binary_C" | scenario == "continuous_Y_continuous_C"){
        otu.all[[i.seed]] <- rbind(unlist(otu.locom), unlist(otu.aldex2))
        rownames(otu.all[[i.seed]]) <- c("LOCOM", "ALDEx2")
    }
    colnames(otu.all[[i.seed]]) <- c("n.otu", "sen", "sep", "fdr")
  
}

# Output

p.global.summary <- Reduce("+", lapply(p.global.all, function(x) as.numeric(x < 0.05)))/length(p.global.all)
filename.p.global <- paste("S_", scenario, "_cau", causal.type, "_lib", library.mu, "_n", n.sam, "_filter", filter.thresh, "_fdr", fdr.target, "_beta", beta, "_global.txt", sep = "")
write.table(p.global.summary, file = filename.p.global)

otu.summary <- Reduce("+", otu.all)/length(otu.all)
filename.otu <- paste("S_", scenario, "_cau", causal.type, "_lib", library.mu, "_n", n.sam, "_filter", filter.thresh, "_fdr", fdr.target, "_beta", beta, "_otu.txt", sep = "")
write.table(otu.summary , file = filename.otu)
