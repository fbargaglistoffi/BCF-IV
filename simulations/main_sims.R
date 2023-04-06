##################
## Dependencies ##
##################

rm(list=ls())
library(BayesIV)
library(doParallel)
library(foreach)
library(dplyr)

###########################
## Initialize Parameters ##
###########################

n = 1000 # number of data points
p = 10 # number of covariates 
mu = rep(0, p)
rho = 0 # correlation within the covariates
null = 0
seq <- seq(0, 2, 0.2) # effect size
nsim = 1000 # number of datasets created
compliance = 0.75 # compliance rate

# Set up Parallel Computation
# Setup parallel backend to use many processors
cl <- makeCluster(20)
registerDoParallel(cl)

# This funnctio takes an arbitrary number of lists x all of which much have the same structure    
comb <- function(x, ...) {  
  mapply(rbind,x,...,SIMPLIFY=FALSE)
}

##########################################
##            1000 OBSERVATIONS         ##
##########################################

set.seed(2020)
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
    rules_double_bart <- data.frame()
    rules_bcf_iv <- data.frame()
    rules_hctree <- data.frame()
    tau_grf00 <- tau_grf11 <- data.frame()
    tau_bcf_iv00 <- tau_bcf_iv11 <- data.frame()
    se_grf00 <- se_grf11 <- data.frame()
    se_bcf_iv00 <- se_bcf_iv11 <- data.frame()
    
    for (i in seq)   {
      ###########################
      # Data Generating Process #
      ##########################
      
      
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
      w1 <- rbinom(n, 1, compliance)  
      w0 <- numeric(n) 
      
      # Generate unit level potential outcome
      y0 <- rnorm(n)  
      y1 <- numeric(n)  
      
      # Generate Heterogeneity
      y1[x1==0 & x2==0] <- y0[x1==0 & x2==0] + w1[x1==0 & x2==0] * i
      y1[x1==0 & x2==1] <- y0[x1==0 & x2==1] + w1[x1==0 & x2==1] * null
      y1[x1==1 & x2==0] <- y0[x1==1 & x2==0] + w1[x1==1 & x2==0] * null
      y1[x1==1 & x2==1] <- y0[x1==1 & x2==1] + w1[x1==1 & x2==1] * -i
      
      # Generate Random Instrument
      z <- rbinom(n, 1, 0.5) 
      
      #  Unit level observed exposure and observed response
      w <- z * w1 + (1-z) * w0  
      y <- z * y1 + (1-z) * y0 
      
      # Observed data
      dataset <- as.data.frame(cbind(y, z, w, X))
      
      # Implement check on the true discovery for HCT and BCF-IV
      rule.sel.double_bart <- bcf_iv(X, y, z, w)
      rule.sel.bcf.itt <- bcf_itt(X, y, z)
      rule.sel.ctree <- hctree(X, y, z)
      
      # Correct Rules
      rule.sel.double_bart <- length(which(rule.sel.double_bart=="xx1< 0.5 & xx2>=0.5" |
                                             rule.sel.double_bart=="xx1< 0.5 & xx2< 0.5" |
                                             rule.sel.double_bart=="xx2< 0.5 & xx1< 0.5" |
                                             rule.sel.double_bart=="xx2>=0.5 & xx1< 0.5" ))
      correct.rules.bcf.itt <- length(which(rule.sel.bcf.itt=="xx1< 0.5 & xx2>=0.5" |
                                             rule.sel.bcf.itt=="xx1< 0.5 & xx2< 0.5" |
                                             rule.sel.bcf.itt=="xx2< 0.5 & xx1< 0.5" |
                                             rule.sel.bcf.itt=="xx2>=0.5 & xx1< 0.5" ))
      correct.rules.ctree <- length(which(rule.sel.ctree=="xx1< 0.5 & xx2>=0.5" |
                                            rule.sel.ctree=="xx1< 0.5 & xx2< 0.5" |
                                            rule.sel.ctree=="xx2< 0.5 & xx1< 0.5" |
                                            rule.sel.ctree=="xx2>=0.5 & xx1< 0.5" ))
      
      # Store Results for Causal Rules
      rules_double_bart[j, which(seq==i)] <- rule.sel.double_bart 
      rules_bcf_iv[j, which(seq==i)] <- correct.rules.bcf.itt
      rules_hctree[j, which(seq==i)] <- correct.rules.ctree
      
      # Check on Precision
      # GRF
      iv.forest <- instrumental_forest(X, y, w, z)
      
      # Get posterior of treatment effects
      tau.hat <- predict(iv.forest, X, estimate.variance = TRUE)
      tau_grf00[j*2-1, which(seq==i)] <- mean(tau.hat$predictions[which(x1==0 & x2==0)])
      tau_grf11[j*2, which(seq==i)]  <- mean(tau.hat$predictions[which(x1==1 & x2==1)])
      se_grf00[j*2-1, which(seq==i)] <- mean(sqrt(tau.hat$variance[which(x1==0 & x2==0)]))
      se_grf11[j*2, which(seq==i)]  <- mean(sqrt(tau.hat$variance[which(x1==1 & x2==1)]))
      
      
      # Estimation for BCF-IV (2SLS-based)
      subset00 <- dataset[which(x1==0 & x2==0),]
      subset11 <- dataset[which(x1==1 & x2==1),]
      tau_bcf_iv00[j*2-1, which(seq==i)]  <- summary(ivreg(y ~  w  | z , data = subset00))$coef[1,1]
      tau_bcf_iv11[j*2, which(seq==i)]  <- summary(ivreg(y ~  w  | z , data = subset11))$coef[1,1]
      se_bcf_iv00[j*2-1, which(seq==i)]  <- summary(ivreg(y ~  w | z , data = subset00))$coef[1,2]
      se_bcf_iv11[j*2, which(seq==i)]  <- summary(ivreg(y ~  w | z , data = subset11))$coef[1,2]
    }
    # Return the values
    list(rules_bcf_iv, rules_hctree, tau_grf00, tau_grf11, tau_bcf_iv00, tau_bcf_iv11, se_grf00, se_grf11, se_bcf_iv00, se_bcf_iv11, rules_double_bart)
  }
})

# Extract results
rules_bcf <- na.omit(matrix[[1]])
rules_ctree <- na.omit(matrix[[2]])
rules_double_bart <- na.omit(matrix[[11]])
tau_grf_00 <- na.omit(matrix[[3]])
tau_grf_11 <- na.omit(matrix[[4]])
tau_bcf_iv_00 <- na.omit(matrix[[5]])
tau_bcf_iv_11 <- na.omit(matrix[[6]])
se_grf_00 <- na.omit(matrix[[7]])
se_grf_11 <- na.omit(matrix[[8]])
se_bcf_iv_00 <- na.omit(matrix[[9]])
se_bcf_iv_11 <- na.omit(matrix[[10]])

# Get results for Honest Causal Tree.IV
rules_ctree <- colMeans(rules_ctree)

# Get results for BCF-IV
rules_bcf_iv <- colMeans(rules_bcf)

# Get results for Double BART
rules_double_bart <- colMeans(rules_double_bart)


## Mean Squared Error & Bias
mse_tau_grf_00 <- c()
mse_tau_bcf_iv_00 <- c()
bias_tau_grf_00 <- c()
bias_tau_bcf_iv_00 <- c()

for(i in seq){
  mse_tau_grf_00[which(seq==i)] <- mean( (abs(tau_grf_00[,which(seq==i)]) - i)^2 , na.rm = TRUE )
  mse_tau_bcf_iv_00[which(seq==i)] <- mean( (abs(tau_bcf_iv_00[,which(seq==i)]) - i)^2 , na.rm = TRUE )
  bias_tau_grf_00[which(seq==i)] <- mean( (abs(tau_grf_00[,which(seq==i)]) - i) , na.rm = TRUE )
  bias_tau_bcf_iv_00[which(seq==i)] <- mean( (abs(tau_bcf_iv_00[,which(seq==i)]) - i) , na.rm = TRUE )
}

mse_tau_grf_11 <- c()
mse_tau_bcf_iv_11 <- c()
bias_tau_grf_11 <- c()
bias_tau_bcf_iv_11 <- c()

for(i in seq){
  mse_tau_grf_11[which(seq==i)] <- mean( (abs(tau_grf_11[,which(seq==i)]) - i)^2 , na.rm = TRUE )
  mse_tau_bcf_iv_11[which(seq==i)] <- mean( (abs(tau_bcf_iv_11[,which(seq==i)]) - i)^2 , na.rm = TRUE )
  bias_tau_grf_11[which(seq==i)] <- mean( (abs(tau_grf_11[,which(seq==i)]) - i) , na.rm = TRUE )
  bias_tau_bcf_iv_11[which(seq==i)] <- mean( (abs(tau_bcf_iv_11[,which(seq==i)]) - i) , na.rm = TRUE )
}

## Coverage
coverage_tau_grf_00 <- data.frame()
coverage_tau_bcf_iv_00 <- data.frame()

for(i in seq){
  for(j in 1:nrow(tau_grf_00)) {
    if (is.na(tau_grf_00[j,which(seq==i)])==FALSE) {
      coverage_tau_grf_00[j,which(seq==i)] <- between(i, abs(tau_grf_00[j,which(seq==i)]) - 1.96*se_grf_00[j,which(seq==i)],
                                                      abs(tau_grf_00[j,which(seq==i)]) + 1.96*se_grf_00[j,which(seq==i)])
    }
    else {
      coverage_tau_grf_00[j,which(seq==i)] <- NA
    }
  }  
}
coverage_tau_grf_00 <- colMeans(coverage_tau_grf_00, na.rm = TRUE) 

for(i in seq){
  for(j in 1:nrow(tau_bcf_iv_00)) {
    if (is.na(tau_bcf_iv_00[j,which(seq==i)])==FALSE) {
      coverage_tau_bcf_iv_00[j,which(seq==i)] <- between(i, abs(tau_bcf_iv_00[j,which(seq==i)]) - 1.96*se_bcf_iv_00[j,which(seq==i)],
                                                         abs(tau_bcf_iv_00[j,which(seq==i)]) + 1.96*se_bcf_iv_00[j,which(seq==i)])
    }
    else {
      coverage_tau_bcf_iv_00[j,which(seq==i)] <- NA
    }
  }  
}
coverage_tau_bcf_iv_00 <- colMeans(coverage_tau_bcf_iv_00, na.rm = TRUE) 

coverage_tau_grf_11 <- data.frame()
coverage_tau_bcf_iv_11 <- data.frame()

for(i in seq){
  for(j in 1:nrow(tau_grf_11)) {
    if (is.na(tau_grf_11[j,which(seq==i)])==FALSE) {
      coverage_tau_grf_11[j,which(seq==i)] <- between(i, abs(tau_grf_11[j,which(seq==i)]) - 1.96*se_grf_11[j,which(seq==i)],
                                                      abs(tau_grf_11[j,which(seq==i)]) + 1.96*se_grf_11[j,which(seq==i)])
    }
    else {
      coverage_tau_grf_11[j,which(seq==i)] <- NA
    }
  }  
}
coverage_tau_grf_11 <- colMeans(coverage_tau_grf_11, na.rm = TRUE) 

for(i in seq){
  for(j in 1:nrow(tau_bcf_iv_11)) {
    if (is.na(tau_bcf_iv_11[j,which(seq==i)])==FALSE) {
      coverage_tau_bcf_iv_11[j,which(seq==i)] <- between(i, abs(tau_bcf_iv_11[j,which(seq==i)]) - 1.96*se_bcf_iv_11[j,which(seq==i)],
                                                         abs(tau_bcf_iv_11[j,which(seq==i)]) + 1.96*se_bcf_iv_11[j,which(seq==i)])
    }
    else {
      coverage_tau_bcf_iv_11[j,which(seq==i)] <- NA
    }
  }  
}
coverage_tau_bcf_iv_11 <- colMeans(coverage_tau_bcf_iv_11, na.rm = TRUE) 

## Create a Matrix for the Results
results_sims <- cbind(rules_ctree, rules_bcf_iv, rules_double_bart,
                      mse_tau_grf_00, bias_tau_grf_00, coverage_tau_grf_00,
                      mse_tau_grf_11, bias_tau_grf_11, coverage_tau_grf_11,
                      mse_tau_bcf_iv_00, bias_tau_bcf_iv_00, coverage_tau_bcf_iv_00,
                      mse_tau_bcf_iv_11, bias_tau_bcf_iv_11, coverage_tau_bcf_iv_11)

colnames(results_sims) <- c("correct_rules hct-iv",
                            "correct rules bcf-iv", 
                            "correct rules double bart",
                            "MSE GRF 00",
                            "Bias GRF 00",
                            "Coverage GRF 00", 
                            "MSE GRF 11",
                            "Bias GRF 11",
                            "Coverage GRF 11",
                            "MSE BCF-IV 00",
                            "Bias BCF-IV 00",
                            "Coverage BCF-IV 00",
                            "MSE BCF-IV 11",
                            "Bias BCF-IV 11",
                            "Coverage BCF-IV 11")
write.csv(results_sims, file = "sim_0.75_1000.csv")

##########################
##   4000 Observations  ##
##########################

n =  4000

set.seed(2020)
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
    rules_double_bart <- data.frame()
    rules_bcf_iv <- data.frame()
    rules_hctree <- data.frame()
    tau_grf00 <- tau_grf11 <- data.frame()
    tau_bcf_iv00 <- tau_bcf_iv11 <- data.frame()
    se_grf00 <- se_grf11 <- data.frame()
    se_bcf_iv00 <- se_bcf_iv11 <- data.frame()
    
    for (i in seq)   {
      ###########################
      # Data Generating Process #
      ##########################
      
      
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
      w1 <- rbinom(n, 1, compliance)  
      w0 <- numeric(n) 
      
      # Generate unit level potential outcome
      y0 <- rnorm(n)  
      y1 <- numeric(n)  
      
      # Generate Heterogeneity
      y1[x1==0 & x2==0] <- y0[x1==0 & x2==0] + w1[x1==0 & x2==0] * i
      y1[x1==0 & x2==1] <- y0[x1==0 & x2==1] + w1[x1==0 & x2==1] * null
      y1[x1==1 & x2==0] <- y0[x1==1 & x2==0] + w1[x1==1 & x2==0] * null
      y1[x1==1 & x2==1] <- y0[x1==1 & x2==1] + w1[x1==1 & x2==1] * -i
      
      # Generate Random Instrument
      z <- rbinom(n, 1, 0.5) 
      
      #  Unit level observed exposure and observed response
      w <- z * w1 + (1-z) * w0
      y <- z * y1 + (1-z) * y0
      
      # Observed data
      dataset <- as.data.frame(cbind(y, z, w, X))
      
      # Implement check on the true discovery for HCT and BCF-IV
      rule.sel.double_bart <- bcf_iv(X, y, z, w)
      rule.sel.bcf.itt <- bcf_itt(X, y, z)
      rule.sel.ctree <- hctree(X, y, z)
      
      # Correct Rules
      rule.sel.double_bart <- length(which(rule.sel.double_bart=="xx1< 0.5 & xx2>=0.5" |
                                             rule.sel.double_bart=="xx1< 0.5 & xx2< 0.5" |
                                             rule.sel.double_bart=="xx2< 0.5 & xx1< 0.5" |
                                             rule.sel.double_bart=="xx2>=0.5 & xx1< 0.5" ))
      correct.rules.bcf.itt <- length(which(rule.sel.bcf.itt=="xx1< 0.5 & xx2>=0.5" |
                                             rule.sel.bcf.itt=="xx1< 0.5 & xx2< 0.5" |
                                             rule.sel.bcf.itt=="xx2< 0.5 & xx1< 0.5" |
                                             rule.sel.bcf.itt=="xx2>=0.5 & xx1< 0.5" ))
      correct.rules.ctree <- length(which(rule.sel.ctree=="xx1< 0.5 & xx2>=0.5" |
                                            rule.sel.ctree=="xx1< 0.5 & xx2< 0.5" |
                                            rule.sel.ctree=="xx2< 0.5 & xx1< 0.5" |
                                            rule.sel.ctree=="xx2>=0.5 & xx1< 0.5" ))
      
      # Store Results for Causal Rules
      rules_double_bart[j, which(seq==i)] <- rule.sel.double_bart 
      rules_bcf_iv[j, which(seq==i)] <- correct.rules.bcf.itt
      rules_hctree[j, which(seq==i)] <- correct.rules.ctree
      
      # Check on Precision
      # GRF
      iv.forest <- instrumental_forest(X, y, w, z)
      
      # Get posterior of treatment effects
      tau.hat <- predict(iv.forest, X, estimate.variance = TRUE)
      tau_grf00[j*2-1, which(seq==i)] <- mean(tau.hat$predictions[which(x1==0 & x2==0)])
      tau_grf11[j*2, which(seq==i)]  <- mean(tau.hat$predictions[which(x1==1 & x2==1)])
      se_grf00[j*2-1, which(seq==i)] <- mean(sqrt(tau.hat$variance[which(x1==0 & x2==0)]))
      se_grf11[j*2, which(seq==i)]  <- mean(sqrt(tau.hat$variance[which(x1==1 & x2==1)]))
      
      
      # Estimation for BCF-IV (2SLS-based)
      subset00 <- dataset[which(x1==0 & x2==0),]
      subset11 <- dataset[which(x1==1 & x2==1),]
      tau_bcf_iv00[j*2-1, which(seq==i)]  <- summary(ivreg(y ~  w | z , data = subset00))$coef[1,1]
      tau_bcf_iv11[j*2, which(seq==i)]  <- summary(ivreg(y ~  w | z , data = subset11))$coef[1,1]
      se_bcf_iv00[j*2-1, which(seq==i)]  <- summary(ivreg(y ~  w  | z , data = subset00))$coef[1,2]
      se_bcf_iv11[j*2, which(seq==i)]  <- summary(ivreg(y ~  w  | z , data = subset11))$coef[1,2]
    }
    # Return the values
    list(rules_bcf_iv, rules_hctree, tau_grf00, tau_grf11, tau_bcf_iv00, tau_bcf_iv11, se_grf00, se_grf11, se_bcf_iv00, se_bcf_iv11, rules_double_bart)
  }
})


# Extract results
rules_bcf <- na.omit(matrix[[1]])
rules_ctree <- na.omit(matrix[[2]])
rules_double_bart <- na.omit(matrix[[11]])
tau_grf_00 <- na.omit(matrix[[3]])
tau_grf_11 <- na.omit(matrix[[4]])
tau_bcf_iv_00 <- na.omit(matrix[[5]])
tau_bcf_iv_11 <- na.omit(matrix[[6]])
se_grf_00 <- na.omit(matrix[[7]])
se_grf_11 <- na.omit(matrix[[8]])
se_bcf_iv_00 <- na.omit(matrix[[9]])
se_bcf_iv_11 <- na.omit(matrix[[10]])

# Get results for Honest Causal Tree.IV
rules_ctree <- colMeans(rules_ctree)

# Get results for BCF-IV
rules_bcf_iv <- colMeans(rules_bcf)

# Get results for Double BART
rules_double_bart <- colMeans(rules_double_bart)


## Mean Squared Error & Bias
mse_tau_grf_00 <- c()
mse_tau_bcf_iv_00 <- c()
bias_tau_grf_00 <- c()
bias_tau_bcf_iv_00 <- c()

for(i in seq){
  mse_tau_grf_00[which(seq==i)] <- mean( (abs(tau_grf_00[,which(seq==i)]) - i)^2 , na.rm = TRUE )
  mse_tau_bcf_iv_00[which(seq==i)] <- mean( (abs(tau_bcf_iv_00[,which(seq==i)]) - i)^2 , na.rm = TRUE )
  bias_tau_grf_00[which(seq==i)] <- mean( (abs(tau_grf_00[,which(seq==i)]) - i) , na.rm = TRUE )
  bias_tau_bcf_iv_00[which(seq==i)] <- mean( (abs(tau_bcf_iv_00[,which(seq==i)]) - i) , na.rm = TRUE )
}

mse_tau_grf_11 <- c()
mse_tau_bcf_iv_11 <- c()
bias_tau_grf_11 <- c()
bias_tau_bcf_iv_11 <- c()

for(i in seq){
  mse_tau_grf_11[which(seq==i)] <- mean( (abs(tau_grf_11[,which(seq==i)]) - i)^2 , na.rm = TRUE )
  mse_tau_bcf_iv_11[which(seq==i)] <- mean( (abs(tau_bcf_iv_11[,which(seq==i)]) - i)^2 , na.rm = TRUE )
  bias_tau_grf_11[which(seq==i)] <- mean( (abs(tau_grf_11[,which(seq==i)]) - i) , na.rm = TRUE )
  bias_tau_bcf_iv_11[which(seq==i)] <- mean( (abs(tau_bcf_iv_11[,which(seq==i)]) - i) , na.rm = TRUE )
}

## Coverage
coverage_tau_grf_00 <- data.frame()
coverage_tau_bcf_iv_00 <- data.frame()

for(i in seq){
  for(j in 1:nrow(tau_grf_00)) {
    if (is.na(tau_grf_00[j,which(seq==i)])==FALSE) {
      coverage_tau_grf_00[j,which(seq==i)] <- between(i, abs(tau_grf_00[j,which(seq==i)]) - 1.96*se_grf_00[j,which(seq==i)],
                                                      abs(tau_grf_00[j,which(seq==i)]) + 1.96*se_grf_00[j,which(seq==i)])
    }
    else {
      coverage_tau_grf_00[j,which(seq==i)] <- NA
    }
  }  
}
coverage_tau_grf_00 <- colMeans(coverage_tau_grf_00, na.rm = TRUE) 

for(i in seq){
  for(j in 1:nrow(tau_bcf_iv_00)) {
    if (is.na(tau_bcf_iv_00[j,which(seq==i)])==FALSE) {
      coverage_tau_bcf_iv_00[j,which(seq==i)] <- between(i, abs(tau_bcf_iv_00[j,which(seq==i)]) - 1.96*se_bcf_iv_00[j,which(seq==i)],
                                                         abs(tau_bcf_iv_00[j,which(seq==i)]) + 1.96*se_bcf_iv_00[j,which(seq==i)])
    }
    else {
      coverage_tau_bcf_iv_00[j,which(seq==i)] <- NA
    }
  }  
}
coverage_tau_bcf_iv_00 <- colMeans(coverage_tau_bcf_iv_00, na.rm = TRUE) 

coverage_tau_grf_11 <- data.frame()
coverage_tau_bcf_iv_11 <- data.frame()

for(i in seq){
  for(j in 1:nrow(tau_grf_11)) {
    if (is.na(tau_grf_11[j,which(seq==i)])==FALSE) {
      coverage_tau_grf_11[j,which(seq==i)] <- between(i, abs(tau_grf_11[j,which(seq==i)]) - 1.96*se_grf_11[j,which(seq==i)],
                                                      abs(tau_grf_11[j,which(seq==i)]) + 1.96*se_grf_11[j,which(seq==i)])
    }
    else {
      coverage_tau_grf_11[j,which(seq==i)] <- NA
    }
  }  
}
coverage_tau_grf_11 <- colMeans(coverage_tau_grf_11, na.rm = TRUE) 

for(i in seq){
  for(j in 1:nrow(tau_bcf_iv_11)) {
    if (is.na(tau_bcf_iv_11[j,which(seq==i)])==FALSE) {
      coverage_tau_bcf_iv_11[j,which(seq==i)] <- between(i, abs(tau_bcf_iv_11[j,which(seq==i)]) - 1.96*se_bcf_iv_11[j,which(seq==i)],
                                                         abs(tau_bcf_iv_11[j,which(seq==i)]) + 1.96*se_bcf_iv_11[j,which(seq==i)])
    }
    else {
      coverage_tau_bcf_iv_11[j,which(seq==i)] <- NA
    }
  }  
}
coverage_tau_bcf_iv_11 <- colMeans(coverage_tau_bcf_iv_11, na.rm = TRUE) 

## Create a Matrix for the Results
results_sims <- cbind(rules_ctree, rules_bcf_iv, rules_double_bart,
                      mse_tau_grf_00, bias_tau_grf_00, coverage_tau_grf_00,
                      mse_tau_grf_11, bias_tau_grf_11, coverage_tau_grf_11,
                      mse_tau_bcf_iv_00, bias_tau_bcf_iv_00, coverage_tau_bcf_iv_00,
                      mse_tau_bcf_iv_11, bias_tau_bcf_iv_11, coverage_tau_bcf_iv_11)

colnames(results_sims) <- c("correct_rules hct-iv",
                            "correct rules bcf-iv", 
                            "correct rules double bart",
                            "MSE GRF 00",
                            "Bias GRF 00",
                            "Coverage GRF 00", 
                            "MSE GRF 11",
                            "Bias GRF 11",
                            "Coverage GRF 11",
                            "MSE BCF-IV 00",
                            "Bias BCF-IV 00",
                            "Coverage BCF-IV 00",
                            "MSE BCF-IV 11",
                            "Bias BCF-IV 11",
                            "Coverage BCF-IV 11")
write.csv(results_sims, file = "sim_0.75_4000.csv")