###########################
## Initialize  Functions ##
###########################

rm(list=ls())
setwd("...")
source("bcf_iv.R")
library(doParallel)
library(foreach)
library(dplyr)

## Generate Data Functions
# n: number of data points
# p: number of covariates 
# rho: correlation within the covariates
# i = effect_size
# nsim: number of datasets created
# compliance: compliance rate

generate_data <- function(n, p, rho, effect_size, rules, compliance){
  # Generate Variables
  mu = rep(0, p)
  null = 0
  Sigma = matrix(rho, nrow = p, ncol = p) + diag(p)*(1-rho)
  rawvars = mvrnorm(n=n, mu=mu, Sigma=Sigma)
  pvars = pnorm(rawvars)
  binomvars = qbinom(pvars, 1, 0.5) 
  X = binomvars
  x1 = X[,1]
  x2 = X[,2]
  x3 = X[,3]
  x4 = X[,4]
  x5 = X[,5]
  x6 = X[,6]
  x7 = X[,7]
  x8 = X[,8]
  x9 = X[,9]
  x10 = X[,10]
  X <- cbind(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
  
  # Generate unit level observed exposure
  w1 <- rbinom(n, 1, compliance)  
  w0 <- numeric(n) 
  
  # Generate unit level potential outcome
  y0 <- rnorm(n)  
  y1 <- numeric(n)  
  
  # Generate Heterogeneity
  if(rules==2){
    y1[x1==0 & x2==0] <- y0[x1==0 & x2==0] + w1[x1==0 & x2==0] * effect_size
    y1[x1==0 & x2==1] <- y0[x1==0 & x2==1] + w1[x1==0 & x2==1] * null
    y1[x1==1 & x2==0] <- y0[x1==1 & x2==0] + w1[x1==1 & x2==0] * null
    y1[x1==1 & x2==1] <- y0[x1==1 & x2==1] + w1[x1==1 & x2==1] * -effect_size
  }
  if(rules==4){
    y1[x1==0 & x2==0] <- y0[x1==0 & x2==0] + w1[x1==0 & x2==0] * i
    y1[x1==0 & x2==1] <- y0[x1==0 & x2==1] + w1[x1==0 & x2==1] * i*0.5
    y1[x1==1 & x2==0] <- y0[x1==1 & x2==0] + w1[x1==1 & x2==0] * -i*0.5
    y1[x1==1 & x2==1] <- y0[x1==1 & x2==1] + w1[x1==1 & x2==1] * -i
  }
  
  # Generate Random Instrument
  z <- rbinom(n, 1, 0.5) 
  
  #  Unit level observed exposure and observed response
  w <- z * w1 + (1-z) * w0  
  y <- z * y1 + (1-z) * y0 
  
  # Observed data
  dataset <- as.data.frame(cbind(y, z, w, X))
  return(dataset)
}

generate_data_confounding <- function(n, p, rho, effect_size, compliance){
  # Generate Variables
  mu = rep(0, p)
  null = 0
  Sigma = matrix(rho, nrow = p, ncol = p) + diag(p)*(1-rho)
  rawvars = mvrnorm(n=n, mu=mu, Sigma=Sigma)
  pvars = pnorm(rawvars)
  binomvars = qbinom(pvars, 1, 0.5) 
  X = binomvars
  x1 = X[,1]
  x2 = X[,2]
  x3 = X[,3]
  x4 = X[,4]
  x5 = X[,5]
  x6 = X[,6]
  x7 = X[,7]
  x8 = X[,8]
  x9 = X[,9]
  x10 = X[,10]
  X <- cbind(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
  
  # Generate unit level observed exposure
  w1 <- rbinom(n, 1, compliance)  
  w0 <- numeric(n) 
  
  # Generate unit level potential outcome
  y0 <- rnorm(n, mean = x3 + x4, sd=1)  
  y1 <- numeric(n) 
  
  # Generate Heterogeneity
  y1[x1==0 & x2==0] <- y0[x1==0 & x2==0] + w1[x1==0 & x2==0] * i
  y1[x1==0 & x2==1] <- y0[x1==0 & x2==1] + w1[x1==0 & x2==1] * null
  y1[x1==1 & x2==0] <- y0[x1==1 & x2==0] + w1[x1==1 & x2==0] * null
  y1[x1==1 & x2==1] <- y0[x1==1 & x2==1] + w1[x1==1 & x2==1] * -i
  
  # Generate Instrument
  logit.prob = -1 + x3 - x4
  prob = exp(logit.prob)/(1+exp(logit.prob))
  
  # Generate Treatment Indicator
  z = rbinom(n, 1, prob=prob) 
  
  #  Unit level observed exposure and observed response
  w <- z * w1 + (1-z) * w0
  y <- z * y1 + (1-z) * y0
  
  # Observed data
  dataset <- as.data.frame(cbind(y, z, w, X))
  return(dataset)
}



f_score <- function(TP, FP, FN){
  if ( TP + 0.5*(FP + FN) == 0){
    F_score = 0
  } else{
    F_score = (TP/(TP + 0.5*(FP + FN)))
  }
}

fpr <- function(FP, TN){
  if ( (FP + TN) == 0){
    FPR = 0
  } else{
    FPR = (FP/(FP + TN))
  }
}

# Set up number of simulations
nsim = 500

# Set up effect size sequence
seq = seq(0, 2, 0.2)

# Set up Parallel Computation
# Setup parallel backend to use many processors
cl <- makeCluster(10)
registerDoParallel(cl)

# This funnctio takes an arbitrary number of lists x all of which much have the same structure    
comb <- function(x, ...) {  
  mapply(rbind,x,...,SIMPLIFY=FALSE)
}

##########################################
##            1000 OBSERVATIONS         ##
##########################################

# Strong Instrument and Heterogeneity


system.time({
  matrix <- foreach(j=1:nsim,  .combine = 'comb', .multicombine = TRUE) %dopar% {
    
    # Load Packages and Functions
    library(MASS)
    library(bcf)
    library(inTrees)
    library(stabs)
    library(rpart)
    library(grf)
    library(causalTree)
    library(AER)
    
    # Initialize Matrices
    F_score_results <- matrix(NA, nrow = nsim, ncol = length(seq))
    FPR_results <- matrix(NA, nrow = nsim, ncol = length(seq))
    TPR_bcf_iv <- matrix(NA, nrow = nsim, ncol = length(seq))
    TPR_bcf_itt <- matrix(NA, nrow = nsim, ncol = length(seq))
    TPR_ctree <- matrix(NA, nrow = nsim, ncol = length(seq))
    
    for (i in seq)   {
    # Initialize number of rules
    nrules = 2  
    
    # Generate Dataset
    dataset <- generate_data(n = 2000, p = 10, rho = 0, effect_size = i, rules = nrules, compliance = 0.75)
    X <- cbind(dataset$x1, dataset$x2, dataset$x3, dataset$x4, dataset$x5,
               dataset$x6, dataset$x7, dataset$x8, dataset$x9, dataset$x10)
    colnames(X) <- c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10")
    y <- dataset$y
    z <- dataset$z
    w <- dataset$w
    
    # Implement HCT
    rule.sel.ctree <- hctree(X, y, z)
    
    # Implement BCF-ITT
    bcf_itt_results <- bcf_itt(X, y, z)
    rule.sel.bcfitt <- bcf_itt_results$node
    
    # Implement BCF-IV
    bcf_iv_results <- bcf_iv(y, w, z, X)
    rule.sel.bcfiv <- bcf_iv_results$node
    adj.pvalue <- bcf_iv_results$Adj_pvalue
    
    # Total Rules
    TR <- length(which(!is.na(adj.pvalue)))
    
    # BCF-IV Diagnostics
    # True Positives
    TP_bcfiv <- length(which(rule.sel.bcfiv[adj.pvalue <= 0.05]=="x1>=0.5 & x2>=0.5" |
                     rule.sel.bcfiv[adj.pvalue <= 0.05]=="x1< 0.5 & x2< 0.5" |
                     rule.sel.bcfiv[adj.pvalue <= 0.05]=="x2< 0.5 & x1< 0.5" |
                     rule.sel.bcfiv[adj.pvalue <= 0.05]=="x2>=0.5 & x1>=0.5" ))
  
    # False Negatives
    FN_bcfiv <- length(which(rule.sel.bcfiv[adj.pvalue > 0.05]=="x1>=0.5 & x2>=0.5" |
                       rule.sel.bcfiv[adj.pvalue > 0.05]=="x1< 0.5 & x2< 0.5" |
                       rule.sel.bcfiv[adj.pvalue > 0.05]=="x2< 0.5 & x1< 0.5" |
                       rule.sel.bcfiv[adj.pvalue > 0.05]=="x2>=0.5 & x1>=0.5" ))
  
    # False Positives
    FP_bcfiv = abs(TR - TP_bcfiv - length(which(!adj.pvalue <= 0.05)))

    # True Negatives
    TN_bcfiv = abs(TR - FN_bcfiv - length(which(!adj.pvalue > 0.05)))
    
    # True Positive Rate
    TPR_bcfiv = TP_bcfiv/nrules
  
    # F_score
    F_score <- f_score(TP = TP_bcfiv, FP = FP_bcfiv, FN = FN_bcfiv) 
  
    # False Positive Rate
    FPR <- fpr(FP = FP_bcfiv, TN = TN_bcfiv)
  
    # Store Results for Causal Rules
    F_score_results[j, which(seq==i)] <- F_score
  
    # Store Results False Positive Rate
    FPR_results[j, which(seq==i)] <- FPR
    
    # Store Results for True Positive Rate
    TPR_bcf_iv[j, which(seq==i)] <- TPR_bcfiv
    
    # Honest Causal Tree with IV Diagnostics
    TPctree <- length(which(rule.sel.ctree=="xx1>=0.5 & xx2>=0.5" |
                                          rule.sel.ctree=="xx1< 0.5 & xx2< 0.5" |
                                          rule.sel.ctree=="xx2< 0.5 & xx1< 0.5" |
                                          rule.sel.ctree=="xx2>=0.5 & xx1>=0.5" ))
    # True Positive Rate
    TPR_ctree[j, which(seq==i)] = TPctree/nrules
    
    # BCF-ITT Diagnostics
    TPbcfitt <- length(which(rule.sel.bcfitt[adj.pvalue <= 0.05]=="x1>=0.5 & x2>=0.5" |
                                        rule.sel.bcfitt[adj.pvalue <= 0.05]=="x1< 0.5 & x2< 0.5" |
                                        rule.sel.bcfitt[adj.pvalue <= 0.05]=="x2< 0.5 & x1< 0.5" |
                                        rule.sel.bcfitt[adj.pvalue <= 0.05]=="x2>=0.5 & x1>=0.5" ))
    # True Positive Rate
    TPR_bcf_itt[j, which(seq==i)] = TPbcfitt/nrules
    
  }
  # Return the values
  list(F_score_results, FPR_results, TPR_bcf_iv, TPR_ctree, TPR_bcf_itt)
  }
})

# Extract results
F_score_results <- na.omit(matrix[[1]])
FPR_results <- na.omit(matrix[[2]])
TPR_bcf_iv <- na.omit(matrix[[3]])
TPR_ctree <- na.omit(matrix[[4]])
TPR_bcf_itt <- na.omit(matrix[[5]])
  
# Get results for F-Score (Strong Instrument)
F_score_strong <- colMeans(F_score_results)

# Get results for FDR (Strong Instrument)
FPR_strong <- colMeans(FPR_results)

# Get results for TPR for BCF-IV (Strong Instrument)
TPR_bcf_iv_strong <- colMeans(TPR_bcf_iv)

# Get results for TPR for HCT-IV (Strong Instrument)
TPR_ctree_strong <- colMeans(TPR_ctree)

# Get results for TPR for BCF-ITT (Strong Instrument)
TPR_bcf_itt_strong <- colMeans(TPR_bcf_itt)


################################################################################
##########################################
##            4000 OBSERVATIONS         ##
##########################################

# Strong Instrument and Heterogeneity


system.time({
  matrix <- foreach(j=1:nsim,  .combine = 'comb', .multicombine = TRUE) %dopar% {
    
    # Load Packages and Functions
    library(MASS)
    library(bcf)
    library(inTrees)
    library(stabs)
    library(rpart)
    library(grf)
    library(causalTree)
    library(AER)
    
    # Initialize Matrices
    F_score_results <- matrix(NA, nrow = nsim, ncol = length(seq))
    FPR_results <- matrix(NA, nrow = nsim, ncol = length(seq))
    TPR_bcf_iv <- matrix(NA, nrow = nsim, ncol = length(seq))
    TPR_bcf_itt <- matrix(NA, nrow = nsim, ncol = length(seq))
    TPR_ctree <- matrix(NA, nrow = nsim, ncol = length(seq))
    
    for (i in seq)   {
      # Initialize number of rules
      nrules = 2  
      
      # Generate Dataset
      dataset <- generate_data(n = 8000, p = 10, rho = 0, effect_size = i, rules = nrules, compliance = 0.75)
      X <- cbind(dataset$x1, dataset$x2, dataset$x3, dataset$x4, dataset$x5,
                 dataset$x6, dataset$x7, dataset$x8, dataset$x9, dataset$x10)
      colnames(X) <- c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10")
      y <- dataset$y
      z <- dataset$z
      w <- dataset$w
      
      # Implement HCT
      rule.sel.ctree <- hctree(X, y, z)
      
      # Implement BCF-ITT
      bcf_itt_results <- bcf_itt(X, y, z)
      rule.sel.bcfitt <- bcf_itt_results$node
      
      # Implement BCF-IV
      bcf_iv_results <- bcf_iv(y, w, z, X)
      rule.sel.bcfiv <- bcf_iv_results$node
      adj.pvalue <- bcf_iv_results$Adj_pvalue
      
      # Total Rules
      TR <- length(which(!is.na(adj.pvalue)))
      
      # BCF-IV Diagnostics
      # True Positives
      TP_bcfiv <- length(which(rule.sel.bcfiv[adj.pvalue <= 0.05]=="x1>=0.5 & x2>=0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x1< 0.5 & x2< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x2< 0.5 & x1< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x2>=0.5 & x1>=0.5" ))
      
      # False Negatives
      FN_bcfiv <- length(which(rule.sel.bcfiv[adj.pvalue > 0.05]=="x1>=0.5 & x2>=0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x1< 0.5 & x2< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x2< 0.5 & x1< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x2>=0.5 & x1>=0.5" ))
      
      # False Positives
      FP_bcfiv = abs(TR - TP_bcfiv - length(which(!adj.pvalue <= 0.05)))
      
      # True Negatives
      TN_bcfiv = abs(TR - FN_bcfiv - length(which(!adj.pvalue > 0.05)))
      
      # True Positive Rate
      TPR_bcfiv = TP_bcfiv/nrules
      
      # F_score
      F_score <- f_score(TP = TP_bcfiv, FP = FP_bcfiv, FN = FN_bcfiv) 
      
      # False Positive Rate
      FPR <- fpr(FP = FP_bcfiv, TN = TN_bcfiv)
      
      # Store Results for Causal Rules
      F_score_results[j, which(seq==i)] <- F_score
      
      # Store Results False Positive Rate
      FPR_results[j, which(seq==i)] <- FPR
      
      # Store Results for True Positive Rate
      TPR_bcf_iv[j, which(seq==i)] <- TPR_bcfiv
      
      # Honest Causal Tree with IV Diagnostics
      TPctree <- length(which(rule.sel.ctree=="xx1>=0.5 & xx2>=0.5" |
                                rule.sel.ctree=="xx1< 0.5 & xx2< 0.5" |
                                rule.sel.ctree=="xx2< 0.5 & xx1< 0.5" |
                                rule.sel.ctree=="xx2>=0.5 & xx1>=0.5" ))
      # True Positive Rate
      TPR_ctree[j, which(seq==i)] = TPctree/nrules
      
      # BCF-ITT Diagnostics
      TPbcfitt <- length(which(rule.sel.bcfitt[adj.pvalue <= 0.05]=="x1>=0.5 & x2>=0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x1< 0.5 & x2< 0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x2< 0.5 & x1< 0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x2>=0.5 & x1>=0.5" ))
      # True Positive Rate
      TPR_bcf_itt[j, which(seq==i)] = TPbcfitt/nrules
      
    }
    # Return the values
    list(F_score_results, FPR_results, TPR_bcf_iv, TPR_ctree, TPR_bcf_itt)
  }
})

# Extract results
F_score_results <- na.omit(matrix[[1]])
FPR_results <- na.omit(matrix[[2]])
TPR_bcf_iv <- na.omit(matrix[[3]])
TPR_ctree <- na.omit(matrix[[4]])
TPR_bcf_itt <- na.omit(matrix[[5]])

# Get results for F-Score (Strong Instrument)
F_score_strong_4000 <- colMeans(F_score_results)

# Get results for FDR (Strong Instrument)
FPR_strong_4000 <- colMeans(FPR_results)

# Get results for TPR for BCF-IV (Strong Instrument)
TPR_bcf_iv_strong_4000 <- colMeans(TPR_bcf_iv)

# Get results for TPR for HCT-IV (Strong Instrument)
TPR_ctree_strong_4000 <- colMeans(TPR_ctree)

# Get results for TPR for BCF-ITT (Strong Instrument)
TPR_bcf_itt_strong_4000 <- colMeans(TPR_bcf_itt)

##########################################
##            4000 OBSERVATIONS         ##
##########################################

################################################################################
# Mild Instrument

system.time({
  matrix <- foreach(j=1:nsim,  .combine = 'comb', .multicombine = TRUE) %dopar% {
    
    # Load Packages and Functions
    library(MASS)
    library(bcf)
    library(inTrees)
    library(stabs)
    library(rpart)
    library(grf)
    library(causalTree)
    library(AER)
    
    # Initialize Matrices
    F_score_results <- matrix(NA, nrow = nsim, ncol = length(seq))
    FPR_results <- matrix(NA, nrow = nsim, ncol = length(seq))
    TPR_bcf_iv <- matrix(NA, nrow = nsim, ncol = length(seq))
    TPR_bcf_itt <- matrix(NA, nrow = nsim, ncol = length(seq))
    TPR_ctree <- matrix(NA, nrow = nsim, ncol = length(seq))
    
    for (i in seq)   {
      # Initialize number of rules
      nrules = 2  
      
      # Generate Dataset
      dataset <- generate_data(n = 2000, p = 10, rho = 0, effect_size = i, rules = nrules, compliance = 0.5)
      X <- cbind(dataset$x1, dataset$x2, dataset$x3, dataset$x4, dataset$x5,
                 dataset$x6, dataset$x7, dataset$x8, dataset$x9, dataset$x10)
      colnames(X) <- c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10")
      y <- dataset$y
      z <- dataset$z
      w <- dataset$w
      
      # Implement HCT
      rule.sel.ctree <- hctree(X, y, z)
      
      # Implement BCF-ITT
      bcf_itt_results <- bcf_itt(X, y, z)
      rule.sel.bcfitt <- bcf_itt_results$node
      
      # Implement BCF-IV
      bcf_iv_results <- bcf_iv(y, w, z, X)
      rule.sel.bcfiv <- bcf_iv_results$node
      adj.pvalue <- bcf_iv_results$Adj_pvalue
      
      # Total Rules
      TR <- length(which(!is.na(adj.pvalue)))
      
      # BCF-IV Diagnostics
      # True Positives
      TP_bcfiv <- length(which(rule.sel.bcfiv[adj.pvalue <= 0.05]=="x1>=0.5 & x2>=0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x1< 0.5 & x2< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x2< 0.5 & x1< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x2>=0.5 & x1>=0.5" ))
      
      # False Negatives
      FN_bcfiv <- length(which(rule.sel.bcfiv[adj.pvalue > 0.05]=="x1>=0.5 & x2>=0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x1< 0.5 & x2< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x2< 0.5 & x1< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x2>=0.5 & x1>=0.5" ))
      
      # False Positives
      FP_bcfiv = abs(TR - TP_bcfiv - length(which(!adj.pvalue <= 0.05)))
      
      # True Negatives
      TN_bcfiv = abs(TR - FN_bcfiv - length(which(!adj.pvalue > 0.05)))
      
      # True Positive Rate
      TPR_bcfiv = TP_bcfiv/nrules
      
      # F_score
      F_score <- f_score(TP = TP_bcfiv, FP = FP_bcfiv, FN = FN_bcfiv) 
      
      # False Positive Rate
      FPR <- fpr(FP = FP_bcfiv, TN = TN_bcfiv)
      
      # Store Results for Causal Rules
      F_score_results[j, which(seq==i)] <- F_score
      
      # Store Results False Positive Rate
      FPR_results[j, which(seq==i)] <- FPR
      
      # Store Results for True Positive Rate
      TPR_bcf_iv[j, which(seq==i)] <- TPR_bcfiv
      
      # Honest Causal Tree with IV Diagnostics
      TPctree <- length(which(rule.sel.ctree=="xx1>=0.5 & xx2>=0.5" |
                                rule.sel.ctree=="xx1< 0.5 & xx2< 0.5" |
                                rule.sel.ctree=="xx2< 0.5 & xx1< 0.5" |
                                rule.sel.ctree=="xx2>=0.5 & xx1>=0.5" ))
      # True Positive Rate
      TPR_ctree[j, which(seq==i)] = TPctree/nrules
      
      # BCF-ITT Diagnostics
      TPbcfitt <- length(which(rule.sel.bcfitt[adj.pvalue <= 0.05]=="x1>=0.5 & x2>=0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x1< 0.5 & x2< 0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x2< 0.5 & x1< 0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x2>=0.5 & x1>=0.5" ))
      # True Positive Rate
      TPR_bcf_itt[j, which(seq==i)] = TPbcfitt/nrules
      
    }
    # Return the values
    list(F_score_results, FPR_results, TPR_bcf_iv, TPR_ctree, TPR_bcf_itt)
  }
})

# Extract results
F_score_results <- na.omit(matrix[[1]])
FPR_results <- na.omit(matrix[[2]])
TPR_bcf_iv <- na.omit(matrix[[3]])
TPR_ctree <- na.omit(matrix[[4]])
TPR_bcf_itt <- na.omit(matrix[[5]])

# Get results for F-Score (mild Instrument)
F_score_mild <- colMeans(F_score_results)

# Get results for FDR (mild Instrument)
FPR_mild <- colMeans(FPR_results)

# Get results for TPR for BCF-IV (mild Instrument)
TPR_bcf_iv_mild <- colMeans(TPR_bcf_iv)

# Get results for TPR for HCT-IV (mild Instrument)
TPR_ctree_mild <- colMeans(TPR_ctree)

# Get results for TPR for BCF-ITT (mild Instrument)
TPR_bcf_itt_mild <- colMeans(TPR_bcf_itt)

################################################################################
# Weak Instrument

system.time({
  matrix <- foreach(j=1:nsim,  .combine = 'comb', .multicombine = TRUE) %dopar% {
    
    # Load Packages and Functions
    library(MASS)
    library(bcf)
    library(inTrees)
    library(stabs)
    library(rpart)
    library(grf)
    library(causalTree)
    library(AER)
    
    # Initialize Matrices
    F_score_results <- matrix(NA, nrow = nsim, ncol = length(seq))
    FPR_results <- matrix(NA, nrow = nsim, ncol = length(seq))
    TPR_bcf_iv <- matrix(NA, nrow = nsim, ncol = length(seq))
    TPR_bcf_itt <- matrix(NA, nrow = nsim, ncol = length(seq))
    TPR_ctree <- matrix(NA, nrow = nsim, ncol = length(seq))
    
    for (i in seq)   {
      # Initialize number of rules
      nrules = 2  
      
      # Generate Dataset
      dataset <- generate_data(n = 2000, p = 10, rho = 0, effect_size = i, rules = nrules, compliance = 0.25)
      X <- cbind(dataset$x1, dataset$x2, dataset$x3, dataset$x4, dataset$x5,
                 dataset$x6, dataset$x7, dataset$x8, dataset$x9, dataset$x10)
      colnames(X) <- c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10")
      y <- dataset$y
      z <- dataset$z
      w <- dataset$w
      
      # Implement HCT
      rule.sel.ctree <- hctree(X, y, z)
      
      # Implement BCF-ITT
      bcf_itt_results <- bcf_itt(X, y, z)
      rule.sel.bcfitt <- bcf_itt_results$node
      
      # Implement BCF-IV
      bcf_iv_results <- bcf_iv(y, w, z, X)
      rule.sel.bcfiv <- bcf_iv_results$node
      adj.pvalue <- bcf_iv_results$Adj_pvalue
      
      # Total Rules
      TR <- length(which(!is.na(adj.pvalue)))
      
      # BCF-IV Diagnostics
      # True Positives
      TP_bcfiv <- length(which(rule.sel.bcfiv[adj.pvalue <= 0.05]=="x1>=0.5 & x2>=0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x1< 0.5 & x2< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x2< 0.5 & x1< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x2>=0.5 & x1>=0.5" ))
      
      # False Negatives
      FN_bcfiv <- length(which(rule.sel.bcfiv[adj.pvalue > 0.05]=="x1>=0.5 & x2>=0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x1< 0.5 & x2< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x2< 0.5 & x1< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x2>=0.5 & x1>=0.5" ))
      
      # False Positives
      FP_bcfiv = abs(TR - TP_bcfiv - length(which(!adj.pvalue <= 0.05)))
      
      # True Negatives
      TN_bcfiv = abs(TR - FN_bcfiv - length(which(!adj.pvalue > 0.05)))
      
      # True Positive Rate
      TPR_bcfiv = TP_bcfiv/nrules
      
      # F_score
      F_score <- f_score(TP = TP_bcfiv, FP = FP_bcfiv, FN = FN_bcfiv) 
      
      # False Positive Rate
      FPR <- fpr(FP = FP_bcfiv, TN = TN_bcfiv)
      
      # Store Results for Causal Rules
      F_score_results[j, which(seq==i)] <- F_score
      
      # Store Results False Positive Rate
      FPR_results[j, which(seq==i)] <- FPR
      
      # Store Results for True Positive Rate
      TPR_bcf_iv[j, which(seq==i)] <- TPR_bcfiv
      
      # Honest Causal Tree with IV Diagnostics
      TPctree <- length(which(rule.sel.ctree=="xx1>=0.5 & xx2>=0.5" |
                                rule.sel.ctree=="xx1< 0.5 & xx2< 0.5" |
                                rule.sel.ctree=="xx2< 0.5 & xx1< 0.5" |
                                rule.sel.ctree=="xx2>=0.5 & xx1>=0.5" ))
      # True Positive Rate
      TPR_ctree[j, which(seq==i)] = TPctree/nrules
      
      # BCF-ITT Diagnostics
      TPbcfitt <- length(which(rule.sel.bcfitt[adj.pvalue <= 0.05]=="x1>=0.5 & x2>=0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x1< 0.5 & x2< 0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x2< 0.5 & x1< 0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x2>=0.5 & x1>=0.5" ))
      # True Positive Rate
      TPR_bcf_itt[j, which(seq==i)] = TPbcfitt/nrules
      
    }
    # Return the values
    list(F_score_results, FPR_results, TPR_bcf_iv, TPR_ctree, TPR_bcf_itt)
  }
})

# Extract results
F_score_results <- na.omit(matrix[[1]])
FPR_results <- na.omit(matrix[[2]])
TPR_bcf_iv <- na.omit(matrix[[3]])
TPR_ctree <- na.omit(matrix[[4]])
TPR_bcf_itt <- na.omit(matrix[[5]])

# Get results for F-Score (weak Instrument)
F_score_weak <- colMeans(F_score_results)

# Get results for FDR (weak Instrument)
FPR_weak <- colMeans(FPR_results)

# Get results for TPR for BCF-IV (weak Instrument)
TPR_bcf_iv_weak <- colMeans(TPR_bcf_iv)

# Get results for TPR for HCT-IV (weak Instrument)
TPR_ctree_weak <- colMeans(TPR_ctree)

# Get results for TPR for BCF-ITT (weak Instrument)
TPR_bcf_itt_weak <- colMeans(TPR_bcf_itt)

################################################################################
# Slight Heterogeneity


system.time({
  matrix <- foreach(j=1:nsim,  .combine = 'comb', .multicombine = TRUE) %dopar% {
    
    # Load Packages and Functions
    library(MASS)
    library(bcf)
    library(inTrees)
    library(stabs)
    library(rpart)
    library(grf)
    library(causalTree)
    library(AER)
    
    # Initialize Matrices
    F_score_results <- matrix(NA, nrow = nsim, ncol = length(seq))
    FPR_results <- matrix(NA, nrow = nsim, ncol = length(seq))
    TPR_bcf_iv <- matrix(NA, nrow = nsim, ncol = length(seq))
    TPR_bcf_itt <- matrix(NA, nrow = nsim, ncol = length(seq))
    TPR_ctree <- matrix(NA, nrow = nsim, ncol = length(seq))
    
    for (i in seq)   {
      # Initialize number of rules
      nrules = 4  
      
      # Generate Dataset
      dataset <- generate_data(n = 2000, p = 10, rho = 0, effect_size = i, rules = nrules, compliance = 0.75)
      X <- cbind(dataset$x1, dataset$x2, dataset$x3, dataset$x4, dataset$x5,
                 dataset$x6, dataset$x7, dataset$x8, dataset$x9, dataset$x10)
      colnames(X) <- c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10")
      y <- dataset$y
      z <- dataset$z
      w <- dataset$w
      
      # Implement HCT
      rule.sel.ctree <- hctree(X, y, z)
      
      # Implement BCF-ITT
      bcf_itt_results <- bcf_itt(X, y, z, inference_ratio = 0.2)
      rule.sel.bcfitt <- bcf_itt_results$node
      
      # Implement BCF-IV
      bcf_iv_results <- bcf_iv(y, w, z, X, inference_ratio = 0.2)
      rule.sel.bcfiv <- bcf_iv_results$node
      adj.pvalue <- bcf_iv_results$Adj_pvalue
      
      # Total Rules
      TR <- length(which(!is.na(adj.pvalue)))
      
      # BCF-IV Diagnostics
      
      # True Positives
      TP_bcfiv <- length(which(rule.sel.bcfiv[adj.pvalue <= 0.05]=="x1>=0.5 & x2>=0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x1< 0.5 & x2< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x2< 0.5 & x1< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x2>=0.5 & x1>=0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x1< 0.5 & x2>=0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x1>=0.5 & x2< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x2< 0.5 & x1>=0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x2>=0.5 & x1< 0.5"))
      
      
      # False Negatives
      FN_bcfiv <- length(which(rule.sel.bcfiv[adj.pvalue > 0.05]=="x1>=0.5 & x2>=0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x1< 0.5 & x2< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x2< 0.5 & x1< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x2>=0.5 & x1>=0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x1< 0.5 & x2>=0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x1>=0.5 & x2< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x2< 0.5 & x1>=0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x2>=0.5 & x1< 0.5"))
      
      # False Positives
      FP_bcfiv = abs(TR - TP_bcfiv - length(which(!adj.pvalue <= 0.05)))
      
      # True Negatives
      TN_bcfiv = abs(TR - FN_bcfiv - length(which(!adj.pvalue > 0.05)))
      
      # True Positive Rate
      TPR_bcfiv = TP_bcfiv/nrules
      
      # F_score
      F_score <- f_score(TP = TP_bcfiv, FP = FP_bcfiv, FN = FN_bcfiv) 
      
      # False Positive Rate
      FPR <- fpr(FP = FP_bcfiv, TN = TN_bcfiv)
      
      # Store Results for Causal Rules
      F_score_results[j, which(seq==i)] <- F_score
      
      # Store Results False Positive Rate
      FPR_results[j, which(seq==i)] <- FPR
      
      # Store Results for True Positive Rate
      TPR_bcf_iv[j, which(seq==i)] <- TPR_bcfiv
      
      # Honest Causal Tree with IV Diagnostics
      TPctree <- length(which(rule.sel.ctree=="xx1< 0.5 & xx2>=0.5" |
                                rule.sel.ctree=="xx1< 0.5 & xx2< 0.5" |
                                rule.sel.ctree=="xx2< 0.5 & xx1< 0.5" |
                                rule.sel.ctree=="xx2>=0.5 & xx1< 0.5" |
                                rule.sel.ctree=="xx1>=0.5 & xx2>=0.5" |
                                rule.sel.ctree=="xx1>=0.5 & xx2< 0.5" |
                                rule.sel.ctree=="xx2< 0.5 & xx1>=0.5" |
                                rule.sel.ctree=="xx2>=0.5 & xx1>=0.5" ))
      
      # True Positive Rate
      TPR_ctree[j, which(seq==i)] = TPctree/nrules
      
      # BCF-ITT Diagnostics
      TPbcfitt <- length(which(rule.sel.bcfitt[adj.pvalue <= 0.05]=="x1>=0.5 & x2>=0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x1< 0.5 & x2< 0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x2< 0.5 & x1< 0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x2>=0.5 & x1>=0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x1< 0.5 & x2>=0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x1>=0.5 & x2< 0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x2< 0.5 & x1>=0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x2>=0.5 & x1< 0.5"))
      
      
      # True Positive Rate
      TPR_bcf_itt[j, which(seq==i)] = TPbcfitt/nrules
      
    }
    # Return the values
    list(F_score_results, FPR_results, TPR_bcf_iv, TPR_ctree, TPR_bcf_itt)
  }
})

# Extract results
F_score_results <- na.omit(matrix[[1]])
FPR_results <- na.omit(matrix[[2]])
TPR_bcf_iv <- na.omit(matrix[[3]])
TPR_ctree <- na.omit(matrix[[4]])
TPR_bcf_itt <- na.omit(matrix[[5]])

# Get results for F-Score (slight Instrument)
F_score_slight <- colMeans(F_score_results)

# Get results for FDR (slight Instrument)
FPR_slight <- colMeans(FPR_results)

# Get results for TPR for BCF-IV (slight Instrument)
TPR_bcf_iv_slight <- colMeans(TPR_bcf_iv)

# Get results for TPR for HCT-IV (slight Instrument)
TPR_ctree_slight <- colMeans(TPR_ctree)

# Get results for TPR for BCF-ITT (slight Instrument)
TPR_bcf_itt_slight <- colMeans(TPR_bcf_itt)

##########################################
##            4000 OBSERVATIONS         ##
##########################################

################################################################################
# Slight Heterogeneity


system.time({
  matrix <- foreach(j=1:nsim,  .combine = 'comb', .multicombine = TRUE) %dopar% {
    
    # Load Packages and Functions
    library(MASS)
    library(bcf)
    library(inTrees)
    library(stabs)
    library(rpart)
    library(grf)
    library(causalTree)
    library(AER)
    
    # Initialize Matrices
    F_score_results <- matrix(NA, nrow = nsim, ncol = length(seq))
    FPR_results <- matrix(NA, nrow = nsim, ncol = length(seq))
    TPR_bcf_iv <- matrix(NA, nrow = nsim, ncol = length(seq))
    TPR_bcf_itt <- matrix(NA, nrow = nsim, ncol = length(seq))
    TPR_ctree <- matrix(NA, nrow = nsim, ncol = length(seq))
    
    for (i in seq)   {
      # Initialize number of rules
      nrules = 4  
      
      # Generate Dataset
      dataset <- generate_data(n = 8000, p = 10, rho = 0, effect_size = i, rules = nrules, compliance = 0.75)
      X <- cbind(dataset$x1, dataset$x2, dataset$x3, dataset$x4, dataset$x5,
                 dataset$x6, dataset$x7, dataset$x8, dataset$x9, dataset$x10)
      colnames(X) <- c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10")
      y <- dataset$y
      z <- dataset$z
      w <- dataset$w
      
      # Implement HCT
      rule.sel.ctree <- hctree(X, y, z)
      
      # Implement BCF-ITT
      bcf_itt_results <- bcf_itt(X, y, z)
      rule.sel.bcfitt <- bcf_itt_results$node
      
      # Implement BCF-IV
      bcf_iv_results <- bcf_iv(y, w, z, X)
      rule.sel.bcfiv <- bcf_iv_results$node
      adj.pvalue <- bcf_iv_results$Adj_pvalue
      
      # Total Rules
      TR <- length(which(!is.na(adj.pvalue)))
      
      # BCF-IV Diagnostics
      
      # True Positives
      TP_bcfiv <- length(which(rule.sel.bcfiv[adj.pvalue <= 0.05]=="x1>=0.5 & x2>=0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x1< 0.5 & x2< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x2< 0.5 & x1< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x2>=0.5 & x1>=0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x1< 0.5 & x2>=0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x1>=0.5 & x2< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x2< 0.5 & x1>=0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x2>=0.5 & x1< 0.5"))
      
      
      # False Negatives
      FN_bcfiv <- length(which(rule.sel.bcfiv[adj.pvalue > 0.05]=="x1>=0.5 & x2>=0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x1< 0.5 & x2< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x2< 0.5 & x1< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x2>=0.5 & x1>=0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x1< 0.5 & x2>=0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x1>=0.5 & x2< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x2< 0.5 & x1>=0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x2>=0.5 & x1< 0.5"))
      
      # False Positives
      FP_bcfiv = abs(TR - TP_bcfiv - length(which(!adj.pvalue <= 0.05)))
      
      # True Negatives
      TN_bcfiv = abs(TR - FN_bcfiv - length(which(!adj.pvalue > 0.05)))
      
      # True Positive Rate
      TPR_bcfiv = TP_bcfiv/nrules
      
      # F_score
      F_score <- f_score(TP = TP_bcfiv, FP = FP_bcfiv, FN = FN_bcfiv) 
      
      # False Positive Rate
      FPR <- fpr(FP = FP_bcfiv, TN = TN_bcfiv)
      
      # Store Results for Causal Rules
      F_score_results[j, which(seq==i)] <- F_score
      
      # Store Results False Positive Rate
      FPR_results[j, which(seq==i)] <- FPR
      
      # Store Results for True Positive Rate
      TPR_bcf_iv[j, which(seq==i)] <- TPR_bcfiv
      
      # Honest Causal Tree with IV Diagnostics
      TPctree <- length(which(rule.sel.ctree=="xx1< 0.5 & xx2>=0.5" |
                                rule.sel.ctree=="xx1< 0.5 & xx2< 0.5" |
                                rule.sel.ctree=="xx2< 0.5 & xx1< 0.5" |
                                rule.sel.ctree=="xx2>=0.5 & xx1< 0.5" |
                                rule.sel.ctree=="xx1>=0.5 & xx2>=0.5" |
                                rule.sel.ctree=="xx1>=0.5 & xx2< 0.5" |
                                rule.sel.ctree=="xx2< 0.5 & xx1>=0.5" |
                                rule.sel.ctree=="xx2>=0.5 & xx1>=0.5" ))
      
      # True Positive Rate
      TPR_ctree[j, which(seq==i)] = TPctree/nrules
      
      # BCF-ITT Diagnostics
      TPbcfitt <- length(which(rule.sel.bcfitt[adj.pvalue <= 0.05]=="x1>=0.5 & x2>=0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x1< 0.5 & x2< 0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x2< 0.5 & x1< 0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x2>=0.5 & x1>=0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x1< 0.5 & x2>=0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x1>=0.5 & x2< 0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x2< 0.5 & x1>=0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x2>=0.5 & x1< 0.5"))
      
      
      # True Positive Rate
      TPR_bcf_itt[j, which(seq==i)] = TPbcfitt/nrules
      
    }
    # Return the values
    list(F_score_results, FPR_results, TPR_bcf_iv, TPR_ctree, TPR_bcf_itt)
  }
})

# Extract results
F_score_results <- na.omit(matrix[[1]])
FPR_results <- na.omit(matrix[[2]])
TPR_bcf_iv <- na.omit(matrix[[3]])
TPR_ctree <- na.omit(matrix[[4]])
TPR_bcf_itt <- na.omit(matrix[[5]])

# Get results for F-Score (slight Instrument)
F_score_slight_4000 <- colMeans(F_score_results)

# Get results for FDR (slight Instrument)
FPR_slight_4000 <- colMeans(FPR_results)

# Get results for TPR for BCF-IV (slight Instrument)
TPR_bcf_iv_slight_4000 <- colMeans(TPR_bcf_iv)

# Get results for TPR for HCT-IV (slight Instrument)
TPR_ctree_slight_4000 <- colMeans(TPR_ctree)

# Get results for TPR for BCF-ITT (slight Instrument)
TPR_bcf_itt_slight_4000 <- colMeans(TPR_bcf_itt)

################################################################################
# Confounding


system.time({
  matrix <- foreach(j=1:nsim,  .combine = 'comb', .multicombine = TRUE) %dopar% {
    
    # Load Packages and Functions
    library(MASS)
    library(bcf)
    library(inTrees)
    library(stabs)
    library(rpart)
    library(grf)
    library(causalTree)
    library(AER)
    
    # Initialize Matrices
    F_score_results <- matrix(NA, nrow = nsim, ncol = length(seq))
    FPR_results <- matrix(NA, nrow = nsim, ncol = length(seq))
    TPR_bcf_iv <- matrix(NA, nrow = nsim, ncol = length(seq))
    TPR_bcf_itt <- matrix(NA, nrow = nsim, ncol = length(seq))
    TPR_ctree <- matrix(NA, nrow = nsim, ncol = length(seq))
    
    for (i in seq)   {
      # Initialize number of rules
      nrules = 2  
      
      # Generate Dataset
      dataset <- generate_data_confounding(n = 2000, p = 10, rho = 0, effect_size = i, compliance = 0.75)
      X <- cbind(dataset$x1, dataset$x2, dataset$x3, dataset$x4, dataset$x5,
                 dataset$x6, dataset$x7, dataset$x8, dataset$x9, dataset$x10)
      colnames(X) <- c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10")
      y <- dataset$y
      z <- dataset$z
      w <- dataset$w
      
      # Implement HCT
      rule.sel.ctree <- hctree(X, y, z)
      
      # Implement BCF-ITT
      bcf_itt_results <- bcf_itt_confounding(X, y, z)
      rule.sel.bcfitt <- bcf_itt_results$node
      
      # Implement BCF-IV
      bcf_iv_results <- bcf_iv_confounding(y, w, z, X)
      rule.sel.bcfiv <- bcf_iv_results$node
      adj.pvalue <- bcf_iv_results$Adj_pvalue
      
      # Total Rules
      TR <- length(which(!is.na(adj.pvalue)))
      
      # BCF-IV Diagnostics
      # True Positives
      TP_bcfiv <- length(which(rule.sel.bcfiv[adj.pvalue <= 0.05]=="x1>=0.5 & x2>=0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x1< 0.5 & x2< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x2< 0.5 & x1< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x2>=0.5 & x1>=0.5" ))
      
      # False Negatives
      FN_bcfiv <- length(which(rule.sel.bcfiv[adj.pvalue > 0.05]=="x1>=0.5 & x2>=0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x1< 0.5 & x2< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x2< 0.5 & x1< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x2>=0.5 & x1>=0.5" ))
      
      # False Positives
      FP_bcfiv = abs(TR - TP_bcfiv - length(which(!adj.pvalue <= 0.05)))
      
      # True Negatives
      TN_bcfiv = abs(TR - FN_bcfiv - length(which(!adj.pvalue > 0.05)))
      
      # True Positive Rate
      TPR_bcfiv = TP_bcfiv/nrules
      
      # F_score
      F_score <- f_score(TP = TP_bcfiv, FP = FP_bcfiv, FN = FN_bcfiv) 
      
      # False Positive Rate
      FPR <- fpr(FP = FP_bcfiv, TN = TN_bcfiv)
      
      # Store Results for Causal Rules
      F_score_results[j, which(seq==i)] <- F_score
      
      # Store Results False Positive Rate
      FPR_results[j, which(seq==i)] <- FPR
      
      # Store Results for True Positive Rate
      TPR_bcf_iv[j, which(seq==i)] <- TPR_bcfiv
      
      # Honest Causal Tree with IV Diagnostics
      TPctree <- length(which(rule.sel.ctree=="xx1>=0.5 & xx2>=0.5" |
                                rule.sel.ctree=="xx1< 0.5 & xx2< 0.5" |
                                rule.sel.ctree=="xx2< 0.5 & xx1< 0.5" |
                                rule.sel.ctree=="xx2>=0.5 & xx1>=0.5" ))
      # True Positive Rate
      TPR_ctree[j, which(seq==i)] = TPctree/nrules
      
      # BCF-ITT Diagnostics
      TPbcfitt <- length(which(rule.sel.bcfitt[adj.pvalue <= 0.05]=="x1>=0.5 & x2>=0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x1< 0.5 & x2< 0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x2< 0.5 & x1< 0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x2>=0.5 & x1>=0.5" ))
      # True Positive Rate
      TPR_bcf_itt[j, which(seq==i)] = TPbcfitt/nrules
      
    }
    # Return the values
    list(F_score_results, FPR_results, TPR_bcf_iv, TPR_ctree, TPR_bcf_itt)
  }
})

# Extract results
F_score_results <- na.omit(matrix[[1]])
FPR_results <- na.omit(matrix[[2]])
TPR_bcf_iv <- na.omit(matrix[[3]])
TPR_ctree <- na.omit(matrix[[4]])
TPR_bcf_itt <- na.omit(matrix[[5]])

# Get results for F-Score (confounding Instrument)
F_score_confounding <- colMeans(F_score_results)

# Get results for FDR (confounding Instrument)
FPR_confounding <- colMeans(FPR_results)

# Get results for TPR for BCF-IV (confounding Instrument)
TPR_bcf_iv_confounding <- colMeans(TPR_bcf_iv)

# Get results for TPR for HCT-IV (confounding Instrument)
TPR_ctree_confounding <- colMeans(TPR_ctree)

# Get results for TPR for BCF-ITT (confounding Instrument)
TPR_bcf_itt_confounding <- colMeans(TPR_bcf_itt)

################################################################################
# Covariance


system.time({
  matrix <- foreach(j=1:nsim,  .combine = 'comb', .multicombine = TRUE) %dopar% {
    
    # Load Packages and Functions
    library(MASS)
    library(bcf)
    library(inTrees)
    library(stabs)
    library(rpart)
    library(grf)
    library(causalTree)
    library(AER)
    
    # Initialize Matrices
    F_score_results <- matrix(NA, nrow = nsim, ncol = length(seq))
    FPR_results <- matrix(NA, nrow = nsim, ncol = length(seq))
    TPR_bcf_iv <- matrix(NA, nrow = nsim, ncol = length(seq))
    TPR_bcf_itt <- matrix(NA, nrow = nsim, ncol = length(seq))
    TPR_ctree <- matrix(NA, nrow = nsim, ncol = length(seq))
    
    for (i in seq)   {
      # Initialize number of rules
      nrules = 2  
      
      # Generate Dataset
      dataset <- generate_data(n = 2000, p = 10, rho = 0.25, effect_size = i, rules = nrules, compliance = 0.75)
      X <- cbind(dataset$x1, dataset$x2, dataset$x3, dataset$x4, dataset$x5,
                 dataset$x6, dataset$x7, dataset$x8, dataset$x9, dataset$x10)
      colnames(X) <- c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10")
      y <- dataset$y
      z <- dataset$z
      w <- dataset$w
      
      # Implement HCT
      rule.sel.ctree <- hctree(X, y, z)
      
      # Implement BCF-ITT
      bcf_itt_results <- bcf_itt(X, y, z)
      rule.sel.bcfitt <- bcf_itt_results$node
      
      # Implement BCF-IV
      bcf_iv_results <- bcf_iv(y, w, z, X)
      rule.sel.bcfiv <- bcf_iv_results$node
      adj.pvalue <- bcf_iv_results$Adj_pvalue
      
      # Total Rules
      TR <- length(which(!is.na(adj.pvalue)))
      
      # BCF-IV Diagnostics
      # True Positives
      TP_bcfiv <- length(which(rule.sel.bcfiv[adj.pvalue <= 0.05]=="x1>=0.5 & x2>=0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x1< 0.5 & x2< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x2< 0.5 & x1< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x2>=0.5 & x1>=0.5" ))
      
      # False Negatives
      FN_bcfiv <- length(which(rule.sel.bcfiv[adj.pvalue > 0.05]=="x1>=0.5 & x2>=0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x1< 0.5 & x2< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x2< 0.5 & x1< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x2>=0.5 & x1>=0.5" ))
      
      # False Positives
      FP_bcfiv = abs(TR - TP_bcfiv - length(which(!adj.pvalue <= 0.05)))
      
      # True Negatives
      TN_bcfiv = abs(TR - FN_bcfiv - length(which(!adj.pvalue > 0.05)))
      
      # True Positive Rate
      TPR_bcfiv = TP_bcfiv/nrules
      
      # F_score
      F_score <- f_score(TP = TP_bcfiv, FP = FP_bcfiv, FN = FN_bcfiv) 
      
      # False Positive Rate
      FPR <- fpr(FP = FP_bcfiv, TN = TN_bcfiv)
      
      # Store Results for Causal Rules
      F_score_results[j, which(seq==i)] <- F_score
      
      # Store Results False Positive Rate
      FPR_results[j, which(seq==i)] <- FPR
      
      # Store Results for True Positive Rate
      TPR_bcf_iv[j, which(seq==i)] <- TPR_bcfiv
      
      # Honest Causal Tree with IV Diagnostics
      TPctree <- length(which(rule.sel.ctree=="xx1>=0.5 & xx2>=0.5" |
                                rule.sel.ctree=="xx1< 0.5 & xx2< 0.5" |
                                rule.sel.ctree=="xx2< 0.5 & xx1< 0.5" |
                                rule.sel.ctree=="xx2>=0.5 & xx1>=0.5" ))
      # True Positive Rate
      TPR_ctree[j, which(seq==i)] = TPctree/nrules
      
      # BCF-ITT Diagnostics
      TPbcfitt <- length(which(rule.sel.bcfitt[adj.pvalue <= 0.05]=="x1>=0.5 & x2>=0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x1< 0.5 & x2< 0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x2< 0.5 & x1< 0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x2>=0.5 & x1>=0.5" ))
      # True Positive Rate
      TPR_bcf_itt[j, which(seq==i)] = TPbcfitt/nrules
      
    }
    # Return the values
    list(F_score_results, FPR_results, TPR_bcf_iv, TPR_ctree, TPR_bcf_itt)
  }
})

# Extract results
F_score_results <- na.omit(matrix[[1]])
FPR_results <- na.omit(matrix[[2]])
TPR_bcf_iv <- na.omit(matrix[[3]])
TPR_ctree <- na.omit(matrix[[4]])
TPR_bcf_itt <- na.omit(matrix[[5]])

# Get results for F-Score (covariance Instrument)
F_score_covariance <- colMeans(F_score_results)

# Get results for FDR (covariance Instrument)
FPR_covariance <- colMeans(FPR_results)

# Get results for TPR for BCF-IV (covariance Instrument)
TPR_bcf_iv_covariance <- colMeans(TPR_bcf_iv)

# Get results for TPR for HCT-IV (covariance Instrument)
TPR_ctree_covariance <- colMeans(TPR_ctree)

# Get results for TPR for BCF-ITT (covariance Instrument)
TPR_bcf_itt_covariance <- colMeans(TPR_bcf_itt)

################################################################################
# Misspecified Propensity Score


system.time({
  matrix <- foreach(j=1:nsim,  .combine = 'comb', .multicombine = TRUE) %dopar% {
    
    # Load Packages and Functions
    library(MASS)
    library(bcf)
    library(inTrees)
    library(stabs)
    library(rpart)
    library(grf)
    library(causalTree)
    library(AER)
    
    # Initialize Matrices
    F_score_results <- matrix(NA, nrow = nsim, ncol = length(seq))
    FPR_results <- matrix(NA, nrow = nsim, ncol = length(seq))
    TPR_bcf_iv <- matrix(NA, nrow = nsim, ncol = length(seq))
    TPR_bcf_itt <- matrix(NA, nrow = nsim, ncol = length(seq))
    TPR_ctree <- matrix(NA, nrow = nsim, ncol = length(seq))
    
    for (i in seq)   {
      # Initialize number of rules
      nrules = 2  
      
      # Generate Dataset
      dataset <- generate_data(n = 2000, p = 10, rho = 0, effect_size = i, rules = nrules, compliance = 0.75)
      X <- cbind(dataset$x1, dataset$x2, dataset$x3, dataset$x4, dataset$x5,
                 dataset$x6, dataset$x7, dataset$x8, dataset$x9, dataset$x10)
      colnames(X) <- c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10")
      y <- dataset$y
      z <- dataset$z
      w <- dataset$w
      
      # Implement HCT
      rule.sel.ctree <- hctree(X, y, z)
      
      # Implement BCF-ITT
      bcf_itt_results <- bcf_itt_miss_ps(X, y, z)
      rule.sel.bcfitt <- bcf_itt_results$node
      
      # Implement BCF-IV
      bcf_iv_results <- bcf_iv_miss_ps(y, w, z, X)
      rule.sel.bcfiv <- bcf_iv_results$node
      adj.pvalue <- bcf_iv_results$Adj_pvalue
      
      # Total Rules
      TR <- length(which(!is.na(adj.pvalue)))
      
      # BCF-IV Diagnostics
      # True Positives
      TP_bcfiv <- length(which(rule.sel.bcfiv[adj.pvalue <= 0.05]=="x1>=0.5 & x2>=0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x1< 0.5 & x2< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x2< 0.5 & x1< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x2>=0.5 & x1>=0.5" ))
      
      # False Negatives
      FN_bcfiv <- length(which(rule.sel.bcfiv[adj.pvalue > 0.05]=="x1>=0.5 & x2>=0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x1< 0.5 & x2< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x2< 0.5 & x1< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x2>=0.5 & x1>=0.5" ))
      
      # False Positives
      FP_bcfiv = abs(TR - TP_bcfiv - length(which(!adj.pvalue <= 0.05)))
      
      # True Negatives
      TN_bcfiv = abs(TR - FN_bcfiv - length(which(!adj.pvalue > 0.05)))
      
      # True Positive Rate
      TPR_bcfiv = TP_bcfiv/nrules
      
      # F_score
      F_score <- f_score(TP = TP_bcfiv, FP = FP_bcfiv, FN = FN_bcfiv) 
      
      # False Positive Rate
      FPR <- fpr(FP = FP_bcfiv, TN = TN_bcfiv)
      
      # Store Results for Causal Rules
      F_score_results[j, which(seq==i)] <- F_score
      
      # Store Results False Positive Rate
      FPR_results[j, which(seq==i)] <- FPR
      
      # Store Results for True Positive Rate
      TPR_bcf_iv[j, which(seq==i)] <- TPR_bcfiv
      
      # Honest Causal Tree with IV Diagnostics
      TPctree <- length(which(rule.sel.ctree=="xx1>=0.5 & xx2>=0.5" |
                                rule.sel.ctree=="xx1< 0.5 & xx2< 0.5" |
                                rule.sel.ctree=="xx2< 0.5 & xx1< 0.5" |
                                rule.sel.ctree=="xx2>=0.5 & xx1>=0.5" ))
      # True Positive Rate
      TPR_ctree[j, which(seq==i)] = TPctree/nrules
      
      # BCF-ITT Diagnostics
      TPbcfitt <- length(which(rule.sel.bcfitt[adj.pvalue <= 0.05]=="x1>=0.5 & x2>=0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x1< 0.5 & x2< 0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x2< 0.5 & x1< 0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x2>=0.5 & x1>=0.5" ))
      # True Positive Rate
      TPR_bcf_itt[j, which(seq==i)] = TPbcfitt/nrules
      
    }
    # Return the values
    list(F_score_results, FPR_results, TPR_bcf_iv, TPR_ctree, TPR_bcf_itt)
  }
})

# Extract results
F_score_results <- na.omit(matrix[[1]])
FPR_results <- na.omit(matrix[[2]])
TPR_bcf_iv <- na.omit(matrix[[3]])
TPR_ctree <- na.omit(matrix[[4]])
TPR_bcf_itt <- na.omit(matrix[[5]])

# Get results for F-Score (miss_ps Instrument)
F_score_miss_ps <- colMeans(F_score_results)

# Get results for FDR (miss_ps Instrument)
FPR_miss_ps <- colMeans(FPR_results)

# Get results for TPR for BCF-IV (miss_ps Instrument)
TPR_bcf_iv_miss_ps <- colMeans(TPR_bcf_iv)

# Get results for TPR for HCT-IV (miss_ps Instrument)
TPR_ctree_miss_ps <- colMeans(TPR_ctree)

# Get results for TPR for BCF-ITT (miss_ps Instrument)
TPR_bcf_itt_miss_ps <- colMeans(TPR_bcf_itt)


################################################################################
# Partially Weak Instrument (Scenario 1)
seq <- seq(0.5, 0.75, 0.025)

mu = rep(0, p)
rho = 0
p = 10
n = 2000
null = 0

system.time({
  matrix <- foreach(j=1:nsim,  .combine = 'comb', .multicombine = TRUE) %dopar% {
    
    # Load Packages and Functions
    library(MASS)
    library(bcf)
    library(inTrees)
    library(stabs)
    library(rpart)
    library(grf)
    library(causalTree)
    library(AER)
    
    # Initialize Matrices
    F_score_results <- matrix(NA, nrow = nsim, ncol = length(seq))
    FPR_results <- matrix(NA, nrow = nsim, ncol = length(seq))
    TPR_bcf_iv <- matrix(NA, nrow = nsim, ncol = length(seq))
    TPR_bcf_itt <- matrix(NA, nrow = nsim, ncol = length(seq))
    TPR_ctree <- matrix(NA, nrow = nsim, ncol = length(seq))
    
    for (i in seq)   {
      ###########################
      # Data Generating Process #
      ##########################
      
      nrules = 2
      
      # Generate Variables
      Sigma = matrix(rho, nrow = p, ncol = p) + diag(p)*(1-rho)
      rawvars = mvrnorm(n=n, mu=mu, Sigma=Sigma)
      pvars = pnorm(rawvars)
      binomvars = qbinom(pvars, 1, 0.5) 
      X = binomvars
      x1 = X[,1]
      x2 = X[,2]
      x3 = X[,3]
      x4 = X[,4]
      x5 = X[,5]
      x6 = X[,6]
      x7 = X[,7]
      x8 = X[,8]
      x9 = X[,9]
      x10 = X[,10]
      X <- cbind(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
      
      # Generate unit level observed exposure
      w0 <- numeric(n)
      w1 <- numeric(n)
      w1[x1==0 & x2==0] <- rbinom(n=sum(x1==0 & x2==0), 1, (0.5 -(i-0.5)))
      w1[x1==0 & x2==1] <- rbinom(n=sum(x1==0 & x2==1), 1, 0.5)
      w1[x1==1 & x2==0] <- rbinom(n=sum(x1==1 & x2==0), 1, 0.5)
      w1[x1==1 & x2==1] <- rbinom(n=sum(x1==1 & x2==1), 1, i)
      
      # Generate unit level potential outcome
      y0 <- rnorm(n, mean = 0, sd = 0.5)
      y1 <- numeric(n)  
      
      # Generate Heterogeneity
      y1[x1==0 & x2==0] <- y0[x1==0 & x2==0] + w1[x1==0 & x2==0] * 1
      y1[x1==0 & x2==1] <- y0[x1==0 & x2==1] + w1[x1==0 & x2==1] * null
      y1[x1==1 & x2==0] <- y0[x1==1 & x2==0] + w1[x1==1 & x2==0] * null
      y1[x1==1 & x2==1] <- y0[x1==1 & x2==1] + w1[x1==1 & x2==1] * -1
      
      # Generate Random Instrument
      z <- rbinom(n, 1, 0.5) 
      
      #  Unit level observed exposure and observed response
      w <- z * w1 + (1-z) * w0
      y <- z * y1 + (1-z) * y0
      
      # Observed data
      dataset <- as.data.frame(cbind(y, z, w, X))
      
      # Implement HCT
      rule.sel.ctree <- hctree(X, y, z)
      
      # Implement BCF-ITT
      bcf_itt_results <- bcf_itt(y, w, z, X)
      rule.sel.bcfitt <- bcf_itt_results$node
      
      # Implement BCF-IV
      bcf_iv_results <- bcf_iv(y, w, z, X)
      rule.sel.bcfiv <- bcf_iv_results$node
      adj.pvalue <- bcf_iv_results$Adj_pvalue
      
      # Total Rules
      TR <- length(which(!is.na(adj.pvalue)))
      
      # BCF-IV Diagnostics
      # True Positives
      TP_bcfiv <- length(which(rule.sel.bcfiv[adj.pvalue <= 0.05]=="x1>=0.5 & x2>=0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x1< 0.5 & x2< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x2< 0.5 & x1< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x2>=0.5 & x1>=0.5" ))
      
      # False Negatives
      FN_bcfiv <- length(which(rule.sel.bcfiv[adj.pvalue > 0.05]=="x1>=0.5 & x2>=0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x1< 0.5 & x2< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x2< 0.5 & x1< 0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x2>=0.5 & x1>=0.5" ))
      
      # False Positives
      FP_bcfiv = abs(TR - TP_bcfiv - length(which(!adj.pvalue <= 0.05)))
      
      # True Negatives
      TN_bcfiv = abs(TR - FN_bcfiv - length(which(!adj.pvalue > 0.05)))
      
      # True Positive Rate
      TPR_bcfiv = TP_bcfiv/nrules
      
      # F_score
      F_score <- f_score(TP = TP_bcfiv, FP = FP_bcfiv, FN = FN_bcfiv) 
      
      # False Positive Rate
      FPR <- fpr(FP = FP_bcfiv, TN = TN_bcfiv)
      
      # Store Results for Causal Rules
      F_score_results[j, which(seq==i)] <- F_score
      
      # Store Results False Positive Rate
      FPR_results[j, which(seq==i)] <- FPR
      
      # Store Results for True Positive Rate
      TPR_bcf_iv[j, which(seq==i)] <- TPR_bcfiv
      
      # Honest Causal Tree with IV Diagnostics
      TPctree <- length(which(rule.sel.ctree=="xx1>=0.5 & xx2>=0.5" |
                                rule.sel.ctree=="xx1< 0.5 & xx2< 0.5" |
                                rule.sel.ctree=="xx2< 0.5 & xx1< 0.5" |
                                rule.sel.ctree=="xx2>=0.5 & xx1>=0.5" ))
      # True Positive Rate
      TPR_ctree[j, which(seq==i)] = TPctree/nrules
      
      # BCF-ITT Diagnostics
      TPbcfitt <- length(which(rule.sel.bcfitt[adj.pvalue <= 0.05]=="x1>=0.5 & x2>=0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x1< 0.5 & x2< 0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x2< 0.5 & x1< 0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x2>=0.5 & x1>=0.5" ))
      # True Positive Rate
      TPR_bcf_itt[j, which(seq==i)] = TPbcfitt/nrules
      
    }
    # Return the values
    list(F_score_results, FPR_results, TPR_bcf_iv, TPR_ctree, TPR_bcf_itt)
  }
})

# Extract results
F_score_results <- na.omit(matrix[[1]])
FPR_results <- na.omit(matrix[[2]])
TPR_bcf_iv <- na.omit(matrix[[3]])
TPR_ctree <- na.omit(matrix[[4]])
TPR_bcf_itt <- na.omit(matrix[[5]])

# Get results for F-Score (pw1 Instrument)
F_score_pw1 <- colMeans(F_score_results)

# Get results for FDR (pw1 Instrument)
FPR_pw1 <- colMeans(FPR_results)

# Get results for TPR for BCF-IV (pw1 Instrument)
TPR_bcf_iv_pw1 <- colMeans(TPR_bcf_iv)

# Get results for TPR for HCT-IV (pw1 Instrument)
TPR_ctree_pw1 <- colMeans(TPR_ctree)

# Get results for TPR for BCF-ITT (pw1 Instrument)
TPR_bcf_itt_pw1 <- colMeans(TPR_bcf_itt)

################################################################################
# Partially Weak Instrument (Scenario 2)


system.time({
  matrix <- foreach(j=1:nsim,  .combine = 'comb', .multicombine = TRUE) %dopar% {
    
    # Load Packages and Functions
    library(MASS)
    library(bcf)
    library(inTrees)
    library(stabs)
    library(rpart)
    library(grf)
    library(causalTree)
    library(AER)
    
    # Initialize Matrices
    F_score_results <- matrix(NA, nrow = nsim, ncol = length(seq))
    FPR_results <- matrix(NA, nrow = nsim, ncol = length(seq))
    TPR_bcf_iv <- matrix(NA, nrow = nsim, ncol = length(seq))
    TPR_bcf_itt <- matrix(NA, nrow = nsim, ncol = length(seq))
    TPR_ctree <- matrix(NA, nrow = nsim, ncol = length(seq))
    
    for (i in seq)   {
      ###########################
      # Data Generating Process #
      ##########################
      nrules = 2
      
      # Generate Variables
      Sigma = matrix(rho, nrow = p, ncol = p) + diag(p)*(1-rho)
      rawvars = mvrnorm(n=n, mu=mu, Sigma=Sigma)
      pvars = pnorm(rawvars)
      binomvars = qbinom(pvars, 1, 0.5) 
      X = binomvars
      x1 = X[,1]
      x2 = X[,2]
      x3 = X[,3]
      x4 = X[,4]
      x5 = X[,5]
      X <- cbind(x1, x2, x3, x4, x5)
      
      # Generate unit level observed exposure
      w1 <- numeric(n)
      w1[x1==0] <- rbinom(n=sum(x1==0), 1, (0.5 -(i-0.5)))
      w1[x1==1] <- rbinom(n=sum(x1==1), 1, i)
      w0 <- numeric(n) 
      
      # Generate unit level potential outcomes
      y0 <- rnorm(n, 0, 1)
      
      # Generate Effect
      y1 <- y0 + 1

      # Generate Random Instrument
      z <- rbinom(n, 1, 0.5) 
      
      #  Unit level observed exposure and observed response
      w <- z * w1 + (1-z) * w0 
      y <- z * y1 + (1-z) * y0

      # Observed data
      dataset <- as.data.frame(cbind(y, z, w, X))
      
      # Implement HCT
      rule.sel.ctree <- hctree(X, y, z)
      
      # Implement BCF-ITT
      bcf_itt_results <- bcf_itt(X, y, z, max_depth = 1)
      rule.sel.bcfitt <- bcf_itt_results$node
      
      # Implement BCF-IV
      bcf_iv_results <- bcf_iv(y, w, z, X, max_depth = 1)
      rule.sel.bcfiv <- bcf_iv_results$node
      adj.pvalue <- bcf_iv_results$Adj_pvalue
      
      # Total Rules
      TR <- length(which(!is.na(adj.pvalue)))
      
      # BCF-IV Diagnostics
      # True Positives
      TP_bcfiv <- length(which(rule.sel.bcfiv[adj.pvalue <= 0.05]=="x1>=0.5" |
                                 rule.sel.bcfiv[adj.pvalue <= 0.05]=="x1< 0.5"))
      
      # False Negatives
      FN_bcfiv <- length(which(rule.sel.bcfiv[adj.pvalue > 0.05]=="x1>=0.5" |
                                 rule.sel.bcfiv[adj.pvalue > 0.05]=="x1< 0.5"))
      
      # False Positives
      FP_bcfiv = abs(TR - TP_bcfiv - length(which(!adj.pvalue <= 0.05)))
      
      # True Negatives
      TN_bcfiv = abs(TR - FN_bcfiv - length(which(!adj.pvalue > 0.05)))
      
      # True Positive Rate
      TPR_bcfiv = TP_bcfiv/nrules
      
      # F_score
      F_score <- f_score(TP = TP_bcfiv, FP = FP_bcfiv, FN = FN_bcfiv) 
      
      # False Positive Rate
      FPR <- fpr(FP = FP_bcfiv, TN = TN_bcfiv)
      
      # Store Results for Causal Rules
      F_score_results[j, which(seq==i)] <- F_score
      
      # Store Results False Positive Rate
      FPR_results[j, which(seq==i)] <- FPR
      
      # Store Results for True Positive Rate
      TPR_bcf_iv[j, which(seq==i)] <- TPR_bcfiv
      
      # Honest Causal Tree with IV Diagnostics
      TPctree <- length(which(rule.sel.ctree=="xx1< 0.5" |
                                rule.sel.ctree=="xx1>=0.5"))
      # True Positive Rate
      TPR_ctree[j, which(seq==i)] = TPctree/nrules
      
      # BCF-ITT Diagnostics
      TPbcfitt <- length(which(rule.sel.bcfitt[adj.pvalue <= 0.05]=="x1>=0.5" |
                                 rule.sel.bcfitt[adj.pvalue <= 0.05]=="x1< 0.5"))
      # True Positive Rate
      TPR_bcf_itt[j, which(seq==i)] = TPbcfitt/nrules
      
    }
    # Return the values
    list(F_score_results, FPR_results, TPR_bcf_iv, TPR_ctree, TPR_bcf_itt)
  }
})

# Extract results
F_score_results <- na.omit(matrix[[1]])
FPR_results <- na.omit(matrix[[2]])
TPR_bcf_iv <- na.omit(matrix[[3]])
TPR_ctree <- na.omit(matrix[[4]])
TPR_bcf_itt <- na.omit(matrix[[5]])

# Get results for F-Score (pw1 Instrument)
F_score_pw2 <- colMeans(F_score_results)

# Get results for FDR (pw1 Instrument)
FPR_pw2 <- colMeans(FPR_results)

# Get results for TPR for BCF-IV (pw1 Instrument)
TPR_bcf_iv_pw2 <- colMeans(TPR_bcf_iv)

# Get results for TPR for HCT-IV (pw1 Instrument)
TPR_ctree_pw2 <- colMeans(TPR_ctree)

# Get results for TPR for BCF-ITT (pw1 Instrument)
TPR_bcf_itt_pw2 <- colMeans(TPR_bcf_itt)
