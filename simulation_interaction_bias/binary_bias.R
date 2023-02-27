options(echo=TRUE) # if you want to see commands in output file
options(digits = 10)
args=(commandArgs(TRUE))
print(args)

if(length(args)==0){
  print("No arguments supplied")
  
  test_global <- 1    # test global hypothesis or not?
  test_otu <- 1       # test individual otus or not?
  
  Y_type <- 1         # 1: binary Y (testing variable), 0: continuous Y
  causal_type <- 1   # 1: spike.col=468 (freq=0.10); 2: spike.col=673 (freq=0.05); 3: spike.col=391 (freq=0.01)
  library_mu <- 10000 # mean of library size
  disp <- 0.02        # dispersion parameter
  n_sam <- 100        # total sample size
  alpha <- 5
  gamma <- 2
  freq_thres <- 0
  censor_thres <- 0.8
  ss <- 3
  refer_causal <- 0
  M <- 1
  random <- 0
  phi <- 1
  
  n_sim <- 1
  i_seed <- 11
  
} else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

library(psych)
library(dirmult) # simPop
library(Rcpp)
library(RcppArmadillo)
library(R.utils) 
library(metap)
library(permute)
library(LOCOM)

summarize_otu_results <- function(otuname, qvalue, causal_otus, non_causal_otus, fdr_target=0.1) {
  otu_detected = otuname[which(qvalue < fdr_target)]
  n_otu = length(otu_detected)
  sen = sum(otu_detected %in% causal_otus)/length(causal_otus)
  sep = 1 - sum(otu_detected %in% non_causal_otus)/length(non_causal_otus)
  fdr = n_otu - sum(otu_detected %in% causal_otus)
  (fdr = fdr/n_otu)
  
  out = list(n_otu = n_otu, sen = sen, sep = sep, fdr = fdr, otu_detected = otu_detected)
  return(out)
}

filterOTU <- function(otu.table, freq.table, refer.col, censor.thres, freq.thres, eps = 1){
  # Args:
  #  otu.table: simulated otu table
  #  refer.col: chosen referenec otu idx
  #  censor.thres: threshold used to presence and absence
  #  freq.thre: threshold for mean frequency
  #  eps: pseudo count used when there is any 0 in the reference col
  
  n.sam <- nrow(otu.table)
  n.otus <- ncol(otu.table)
  censoring.rate <- apply(otu.table, 2, function(x) sum(x==0)/n.sam) 
  filter.idx <- which(colMeans(freq.table) < freq.thres | censoring.rate > censor.thres)
  cat("There are total", length(filter.idx), "otu will be filter out", "\n")
  
  if(length(filter.idx)>0){
    otuname <- c(1:n.otus)
    names(otuname) <- c(1:n.otus)
    otuname <- otuname[-filter.idx]
    refer.col <- which(otuname == refer.col)
    otu.table.filter <- otu.table[, -filter.idx]
  } else {
    otu.table.filter <- otu.table
  }
  otu.table.filter[which(otu.table.filter[,refer.col]==0), refer.col] <- eps  
  
  
  return(list(otu.table.filter = otu.table.filter, refer.col = refer.col, 
              otuname = otuname, filter.idx = filter.idx, censoring.rate=censoring.rate))
}

library_sd = library_mu/ss # standard error in library size  #can use normal to simulate library_size
lib_sizer_lower = 2000     # lower bound for library size

for (i_seed in (1+(k-1)* K) : (k* K)) {
  cat("i_seed",i_seed,"\n")
  
  # ----------------------------
  # Set frequency value
  # ----------------------------
  pi <- read.table("./input_throat/fit_dirmult_pi.txt", header=FALSE, as.is=TRUE)[,1]
  n_otus = length(pi)
  
  if(causal_type == 1){
    random.col <- which(pi >= 0.005)[1:20]
    spike.col <- random.col
  } 
  
  if(causal_type == 2){
    random.col <- order(pi, decreasing = TRUE)[1:5]
    spike.col <- random.col
  }
  
  causal_otus_sim = spike.col
  n_otus_causal_sim = length(causal_otus_sim)
  non_causal_otus_sim = setdiff(1:n_otus, causal_otus_sim)
  n_otus_noncausal_sim = length(non_causal_otus_sim)
  
  causal.idx <- rep(0, n_otus)
  causal.idx[causal_otus_sim] <- 1
  # -----------------------------
  
  inter.col <- c(1:n_otus)
  
  # ------------------------
  # simulate bias factor
  # ------------------------
  set.seed(0)
  half.otu <- sample(c(1:n_otus), n_otus/2)
  set.seed(i_seed)
  gamma_otu <- rnorm(n_otus, 0, 0.8)
  gamma_otu_sign <- -sign(t(matrix(rep(gamma_otu, n_otus), nrow = n_otus)))
  if(M == 1){
    epsilon <- 2 * matrix(rnorm(n_otus * length(inter.col), 0.5, 0.1), nrow = length(inter.col))
    gamma_otu_int <- phi * abs(gamma_otu) * epsilon * gamma_otu_sign # 
  }
  if(M == 2){
    epsilon <- 2 * matrix(rnorm(n_otus * length(inter.col), 0.5, 0.1), nrow = length(inter.col))
    epsilon[causal_otus_sim, ] <- 2 * matrix(rbeta(length(causal_otus_sim) * length(inter.col), 0.5, 0.5), ncol = length(inter.col))
    gamma_otu_int <- phi * abs(gamma_otu) * epsilon * gamma_otu_sign # 
  }
  if(M == 3){
    epsilon <- 2 * matrix(rnorm(n_otus * length(inter.col), 0.5, 0.1), nrow = length(inter.col))
    epsilon[half.otu, ] <- 2 * matrix(rbeta(length(half.otu) * length(inter.col), 0.5, 0.5), ncol = length(inter.col))
    gamma_otu_int <- phi * abs(gamma_otu) * epsilon * gamma_otu_sign # 
  }

  diag(gamma_otu_int) <- 0
  
  # -----------------------------
  # Simulate Data
  # -----------------------------
  frac <- 0.5
  set.seed(i_seed)
  refer.col <- integer(0)
  while(length(refer.col) == 0){
    pi0 <- pi
    pi1 <- pi0
    n0_sam <- n_sam * frac
    n1_sam <- n_sam * (1 - frac)
    Y <- c(rep(0, n0_sam), rep(1, n1_sam))
    
    # ------------------------
    # simulation library size
    # ------------------------
    library_size <- rnorm(n_sam, library_mu, library_sd)
    library_size[library_size < lib_sizer_lower] <- lib_sizer_lower
    library_size <- round(library_size)
    
    # ------------------------
    # simulate data
    # ------------------------
    alpha_log <- log(alpha)
    alpha.vector <- rep(0, n_otus)
    alpha.vector[causal_otus_sim] <- alpha_log
    
    a0 <- (1- disp)/disp
    otu_table_sim_0 = dirmult::rdirichlet(n = n0_sam, a0 * pi0)
    otu_table_sim_1 = dirmult::rdirichlet(n = n1_sam, a0 * pi0)
    
    bias_factor_log_mat <- matrix(rep(gamma_otu, n_sam), ncol = n_otus, byrow = TRUE)
    bias_factor_log_inter_mat <- array(0, dim = c(n_sam, n_otus))
    
    # -------------------------
    # introduce bias
    # ------------------------- 
    otu_table_sim <- rbind(otu_table_sim_0, otu_table_sim_1)
    otu_table_sim[ ,spike.col] <- otu_table_sim[ ,spike.col] * exp(alpha_log * Y)
    otu_table_sim <- otu_table_sim/rowSums(otu_table_sim)
    bias_factor_log_inter_mat_normalized <- (otu_table_sim)[ ,inter.col]  %*%  t(gamma_otu_int)
    bias_factor_log_inter_mat <- bias_factor_log_inter_mat_normalized
    
    bias_factor <- exp(bias_factor_log_mat + bias_factor_log_inter_mat)  
    otu_table_sim <- rbind(otu_table_sim_0, otu_table_sim_1) * bias_factor
    
    # --------------------------
    # introduce trait effect
    # --------------------------
    otu_table_sim[ ,spike.col] <- otu_table_sim[ ,spike.col] * exp(alpha_log * Y)
    otu_table_sim <- otu_table_sim/rowSums(otu_table_sim)
    
    otu_table_sim_count <- array(0, dim = c(n_sam, n_otus))
    for (i in 1:n_sam) {
      otu_table_sim_count[i, ] <- rmultinom(1, library_size[i], otu_table_sim[i, ])
    }
    otu_table_sim <- otu_table_sim_count
    freq_table_sim <- otu_table_sim/library_size
    # -----------------------------
    
    # -----------------------------
    # choose reference column
    # -----------------------------
    refer.col.causal <- which.max(apply(freq_table_sim[, causal_otus_sim], 2, mean))
    refer.col.noncausal <- which.max(apply(freq_table_sim[, non_causal_otus_sim], 2, mean))
    refer.col.causal <- causal_otus_sim[refer.col.causal]
    refer.col.noncausal  <- non_causal_otus_sim[refer.col.noncausal]
    
    # -----------------------------
    # Filter Data
    # -----------------------------
    colnames(otu_table_sim) <- c(1:ncol(otu_table_sim))
    if(refer_causal == 0){
      filterData <- filterOTU(otu.table = otu_table_sim, freq.table = freq_table_sim, refer.col = refer.col.noncausal, 
                              censor.thres = censor_thres, freq.thres = freq_thres, eps = 1)
    } else {
      filterData <- filterOTU(otu.table = otu_table_sim, freq.table = freq_table_sim, refer.col = refer.col.causal, 
                              censor.thres = censor_thres, freq.thres = freq_thres, eps = 1)    
    }

    otu.table.filter <-filterData$otu.table.filter
    refer.col <- filterData$refer.col
    otuname <- filterData$otuname
    filter.idx <- filterData$filter.idx
    censoring.rate <- filterData$censoring.rate
  }
  
  causal_otus_sim_filter <- causal_otus_sim[!causal_otus_sim %in% filter.idx]
  non_causal_otus_sim_filter <- non_causal_otus_sim[!non_causal_otus_sim %in% filter.idx]
  cat("number of causal otu:", length(causal_otus_sim_filter), "\n")
  
  # -----------------------------
  # LOCOM estimation
  # -----------------------------
  fdr_target <- 0.2
  ptm <- proc.time()
  res <- locom(otu.table = otu.table.filter, Y = Y, C = NULL,
               ref.otu = refer.col, prev.cut = 0, 
               fdr.nominal = fdr_target, seed = NULL,
               n.perm.max = 50000, n.rej.stop = 100, n.cores = 1)
  proc.time() - ptm
  
  res.summary <- summarize_otu_results(otuname = colnames(otu.table.filter),
                                       qvalue = res$q.otu,
                                       causal_otus = causal_otus_sim_filter,
                                       non_causal_otus  = non_causal_otus_sim_filter,
                                       fdr_target = fdr_target)
  
  n_otu <-  res.summary$n_otu
  sen_otu <-  res.summary$sen
  sep_otu <-  res.summary$sep
  fdr_otu <-  res.summary$fdr
  mat <- c(n_otu, sen_otu, sep_otu, fdr_otu)
  p_global_mat <- res$p.global
  
}
quit(save = "no", status = 0, runLast = TRUE)
