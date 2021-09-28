##### bcf_iv: function for discovery and estimation of heterogeneity in the Complier Average Causal Effect in Instrumental Variable settings
### INPUTs:y: outcome vector (n x 1)
###        w: vector of treatment received  (n x 1)
###        z: vector of treatment assigned (aka instrumental variable)  (n x 1)
###        x: covariates matrix (n x p) 
###        binary: Boolean to identify whether the outcome is binary (TRUE) or continuous/discrete (FALSE) 
###        n_burn: burn-in MCMC iterations for Bayesian Causal Forest
###        n_sim: iterations to save post burn-in for Bayesian Causal Forest
###        inference_ratio: % of observations to be assigned to the inference subsample (default is 0.5)
###        max_depth: maximal depth for the generated CART (default is 2)
###        cp: complexity parameter for the generated CART (default is 0.01)
###        minsplit: minimum observations needed to perform a binary split in the tree (default is 10)
###        adj_method: p-value adjustment method (default is "holm"), other options are "bonferroni", "hockberg", "hommel", "BH", "BY", "fdr", "none"
###        seed: random seed for reproducible results (default is 42)
### OUTPUTS: (1) a tree structure discovering the heterogeneity in the causal effects, (2) TSLS estimates of the Complier Average Causal Effects (CCACE) within its nodes, (3) p-values for each CCACE, (4) adjusted p-value for each CCACE, (5) Intention-to-Treat within all the nodes, (6) proportion of Compliers within all the nodes, (7) p-value for weak-iv test, (8) proportion of observations in the node 

bcf_iv <- function(y, w, z, x, binary = FALSE, n_burn = 500, n_sim = 500, inference_ratio = 0.5, 
                   max_depth = 2, cp = 0.01, minsplit = 10, adj_method = "holm", seed = 42) {
  
  # Upload the Packages
  require(bcf)
  require(rpart)
  require(lattice)
  require(rattle)
  require(AER)
  require(bartCause)
  
  ######################################################
  ####         Step 0: Initialize the Data          ####
  ######################################################
  
  # Split data into Discovery and Inference
  set.seed(seed)
  index <- sample(nrow(x), nrow(x)*inference_ratio, replace=FALSE)
  
  # Initialize total dataset
  iv.data <- as.data.frame(cbind(y, w, z, x))
  names(iv.data)[1:3] <- c("y", "w", "z")
  
  # Discovery and Inference Samples
  discovery <- iv.data[-index,]
  inference <- iv.data[index,]
  
  ######################################################
  ####  Step 1: Compute the Bayesian Causal Forest  ####
  ######################################################
  
  # Compute the Propensity Score though a Logistic Regression
  p.score <- glm(z ~ x[-index,],
                 family = binomial,
                 data = discovery)
  pihat <- predict(p.score, as.data.frame(x[-index,]))
  
  # Perform the Bayesian Causal Forest for the Proportion of Compliers (pic)
  pic_bcf <- quiet(bartc(w[-index], z[-index], x[-index,], n.samples = n_sim, n.burn = n_burn, n.chains = 2L))
  tau_pic <- extract(pic_bcf, type = "ite")
  pic <- apply(tau_pic, 2, mean)
  
  ######################################################
  ####     Continuous and Discrete Outcomes         ####
  ######################################################
  
  if (binary == FALSE){
    
    # Perform the Bayesian Causal Forest for the ITT
    itt_bcf <- quiet(bcf(y[-index], z[-index], x[-index,], x[-index,], pihat, nburn=n_burn, nsim=n_sim))
    tau_itt <- itt_bcf$tau
    itt <- colMeans(tau_itt)
    
    # Get posterior of treatment effects
    tauhat <- itt/pic
    exp <- as.data.frame(cbind(tauhat, x[-index,]))
    
    ######################################################
    ####  Step 2: Build a CART on the Unit Level CITT ####
    ######################################################
    
    fit.tree <- rpart(tauhat ~ .,
                      data = exp,
                      maxdepth = max_depth,
                      cp=cp,
                      minsplit=minsplit)
    
    ######################################################
    ####  Step 3: Extract the Causal Rules (Nodes)    ####
    ######################################################
    
    rules <- as.numeric(row.names(fit.tree$frame[fit.tree$numresp]))
    
    # Initialize Outputs
    bcfivMat <- as.data.frame(matrix(NA, nrow = length(rules), ncol=7))
    names(bcfivMat) <- c("node", "CCACE", "pvalue", "Weak_IV_test", "Pi_obs", "ITT", "Pi_compliers")
    
    # Generate Leaves Indicator
    lvs <- leaves <- numeric(length(rules)) 
    lvs[unique(fit.tree$where)] <- 1
    leaves[rules[lvs==1]] <- 1
    
    ######################################################
    ####  Step 4: Run an IV Regression on each Node   ####
    ######################################################
    
    # Run an IV Regression on the Root
    iv.root <- ivreg(y ~ w | z,  
                     data = inference)
    summary <- summary(iv.root, diagnostics = TRUE)
    iv.effect.root <-  summary$coef[2,1]
    p.value.root <- summary$coef[2,4]
    p.value.weak.iv.root <- summary$diagnostics[1,4]
    proportion.root <- 1
    compliers.root <- length(which(inference$z==inference$w))/nrow(inference)
    itt.root <- iv.effect.root*compliers.root
    
    # Store Results for Root
    bcfivMat[1,] <- c( NA , round(iv.effect.root, 4), round(p.value.root, 4), round(p.value.weak.iv.root, 4), round(proportion.root, 4), round(itt.root, 4), round(compliers.root, 4))
    
    # Initialize New Data
    names(inference) <- paste(names(inference), sep="")
    
    # Run a loop to get the rules (sub-populations)
    for (i in rules[-1]){
      # Create a Vector to Store all the Dimensions of a Rule
      sub <- as.data.frame(matrix(NA, nrow = 1,
                                  ncol = nrow(as.data.frame(path.rpart(fit.tree, node=i, print.it = FALSE)))-1))
      quiet(capture.output(for (j in 1:ncol(sub)){
        # Store each Rule as a Sub-population
        sub[,j] <- as.character(print(as.data.frame(path.rpart(fit.tree,node=i,print.it=FALSE))[j+1,1]))
        sub_pop <- noquote(paste(sub , collapse = " & "))
      }))
      
      subset <- with(inference, inference[which(eval(parse(text=sub_pop))),])
      
      # Run the IV Regression
      if (length(unique(subset$w))!= 1 | length(unique(subset$z))!= 1){
        iv.reg <- ivreg(y ~ w | z,  
                        data = subset)
        summary <- summary(iv.reg, diagnostics = TRUE)
        iv.effect <-  summary$coef[2,1]
        p.value <- summary$coef[2,4]
        p.value.weak.iv <- summary$diagnostics[1,4]
        compliers <- length(which(subset$z==subset$w))/nrow(subset)
        itt <- iv.effect*compliers
        
        # Proportion of observations in the node
        proportion.node <- nrow(subset)/nrow(inference)
        
        ######################################################
        ####   Step 5: Output the Values of each CCACE   ####
        ######################################################
        
        bcfivMat[i,] <- c(sub_pop, round(iv.effect, 4), round(p.value, 4), round(p.value.weak.iv, 4), round(proportion.node, 4), round(itt, 4), round(compliers, 4))
      }
      
      # Delete data
      rm(subset)
    }
    
    # Adjust P.values 
    bcfiv_correction <- cbind(as.data.frame(bcfivMat), leaves)
    adj <- round(p.adjust( as.numeric(bcfiv_correction$pvalue[which(bcfiv_correction$leaves==1)]) ,  paste(adj_method)), 5)
    Adj_pvalue <- rep(NA, length(rules)) 
    Adj_pvalue[which(bcfiv_correction$leaves==1)] <- adj
    
    # Store Results
    bcfivResults <- cbind(as.data.frame(bcfivMat), Adj_pvalue)
  }
  
  ######################################################
  ####           Binary Outcomes (Probit)           ####
  ######################################################
  
  if (binary == TRUE){
    
    # Perform the binary Bayesian Causal Forest for the ITT
    fit_itt <- quiet(bartc(y[-index], z[-index], x[-index,], n.samples = n_sim, n.burn = n_burn, n.chains = 2L))
    
    # Get posterior of treatment effects
    ites <- extract(fit_itt, type = "ite")
    itt <- apply(ites, 2, mean)
    
    # Get posterior of treatment effects
    tauhat <- itt/pic
    exp <- as.data.frame(cbind(tauhat, x[-index,]))
    
    ######################################################
    ####  Step 2: Build a CART on the Unit Level CITT ####
    ######################################################
    
    fit.tree <- rpart(tauhat ~ .,
                      data = exp,
                      maxdepth = max_depth,
                      cp=cp,
                      minsplit=minsplit)
    
    ######################################################
    ####  Step 3: Extract the Causal Rules (Nodes)    ####
    ######################################################
    
    rules <- as.numeric(row.names(fit.tree$frame[fit.tree$numresp]))
    
    # Initialize Outputs
    bcfivMat <- as.data.frame(matrix(NA, nrow = length(rules), ncol=7))
    names(bcfivMat) <- c("node", "CCACE", "pvalue", "Weak_IV_test", "Pi_obs", "ITT", "Pi_compliers")
    
    # Generate Leaves Indicator
    lvs <- leaves <- numeric(length(rules)) 
    lvs[unique(fit.tree$where)] <- 1
    leaves[rules[lvs==1]] <- 1
    
    ######################################################
    ####  Step 4: Run an IV Regression on each Node   ####
    ######################################################
    
    # Run an IV Regression on the Root
    iv.root <- ivreg(y ~ w | z,  
                     data = inference)
    summary <- summary(iv.root, diagnostics = TRUE)
    iv.effect.root <-  summary$coef[2,1]
    p.value.root <- summary$coef[2,4]
    p.value.weak.iv.root <- summary$diagnostics[1,4]
    proportion.root <- 1
    compliers.root <- length(which(inference$z==inference$w))/nrow(inference)
    itt.root <- iv.effect.root*compliers.root
    
    # Store Results for Root
    bcfivMat[1,] <- c( NA , round(iv.effect.root, 4), round(p.value.root, 4), round(p.value.weak.iv.root, 4), round(proportion.root, 4), round(itt.root, 4), round(compliers.root, 4))
    
    # Initialize New Data
    names(inference) <- paste(names(inference), sep="")
    
    # Run a loop to get the rules (sub-populations)
    for (i in rules[-1]){
      # Create a Vector to Store all the Dimensions of a Rule
      sub <- as.data.frame(matrix(NA, nrow = 1,
                                  ncol = nrow(as.data.frame(path.rpart(fit.tree, node=i, print.it = FALSE)))-1))
      quiet(capture.output(for (j in 1:ncol(sub)){
        # Store each Rule as a Sub-population
        sub[,j] <- as.character(print(as.data.frame(path.rpart(fit.tree,node=i,print.it=FALSE))[j+1,1]))
        sub_pop <- noquote(paste(sub , collapse = " & "))
      }))
      
      subset <- with(inference, inference[which(eval(parse(text=sub_pop))),])
      
      # Run the IV Regression
      if (length(unique(subset$w))!= 1 | length(unique(subset$z))!= 1){
        iv.reg <- ivreg(y ~ w | z,  
                        data = subset)
        summary <- summary(iv.reg, diagnostics = TRUE)
        iv.effect <-  summary$coef[2,1]
        p.value <- summary$coef[2,4]
        p.value.weak.iv <- summary$diagnostics[1,4]
        compliers <- length(which(subset$z==subset$w))/nrow(subset)
        itt <- iv.effect*compliers
        
        # Proportion of observations in the node
        proportion.node <- nrow(subset)/nrow(inference)
        
        ######################################################
        ####   Step 5: Output the Values of each CCACE   ####
        ######################################################
        
        bcfivMat[i,] <- c(sub_pop, round(iv.effect, 4), round(p.value, 4), round(p.value.weak.iv, 4), round(proportion.node, 4), round(itt, 4), round(compliers, 4))
      }
      
      # Delete data
      rm(subset)
    }
    
    # Adjust P.values 
    bcfiv_correction <- cbind(as.data.frame(bcfivMat), leaves)
    adj <- round(p.adjust( as.numeric(bcfiv_correction$pvalue[which(bcfiv_correction$leaves==1)]) ,  paste(adj_method)), 5)
    Adj_pvalue <- rep(NA, length(rules)) 
    Adj_pvalue[which(bcfiv_correction$leaves==1)] <- adj
    
    # Store Results
    bcfivResults <- cbind(as.data.frame(bcfivMat), Adj_pvalue)
  }
  
  # Return Results
  return(bcfivResults)
}

##### bcf_itt: function for discovery and estimation of heterogeneity in the Intention-To-Treat in Instrumental Variable settings
### INPUTs:y: outcome vector (n x 1)
###        z: vector of treatment assigned (aka instrumental variable)  (n x 1)
###        x: covariates matrix (n x p) 
###        n_burn: burn-in MCMC iterations for Bayesian Causal Forest
###        n_sim: iterations to save post burn-in for Bayesian Causal Forest
###        inference_ratio: % of observations to be assigned to the inference subsample (default is 0.5)
###        max_depth: maximal depth for the generated CART (default is 2)
###        adj_method: p-value adjustment method (default is "holm"), other options are "bonferroni", "hockberg", "hommel", "BH", "BY", "fdr", "none"
###        seed: random seed for reproducible results (default is 42)
### OUTPUTS: (1) a tree structure discovering the heterogeneity in the causal effects, (2) TSLS estimates of the Complier Average Causal Effects (CCACE) within its nodes, (3) p-values for each CCACE, (4) adjusted p-value for each CCACE, (5) p-value for weak-iv test

bcf_itt <- function(y, z, x, n_burn = 500, n_sim = 500, inference_ratio = 0.5,
                    max_depth = 2, adj_method = "holm", seed = 42) {
  
  ######################################################
  ####         Step 0: Initialize the Data          ####
  ######################################################
  
  # Split data into Discovery and Inference
  set.seed(seed)
  index <- sample(nrow(x), nrow(x)*inference_ratio, replace=FALSE)
  
  # Initialize total dataset
  iv.data <- as.data.frame(cbind(y, w, z, x))
  names(iv.data)[1:3] <- c("y", "w", "z")
  
  # Discovery and Inference Samples
  discovery <- iv.data[-index,]
  inference <- iv.data[index,]
  
  # Compute the Propensity Score though a Logistic Regression
  p.score <- glm(z ~ x[-index,],
                 family = binomial,
                 data = discovery)
  pihat <- predict(p.score, as.data.frame(x[-index,]))
  
  # BART for ITT
  itt_bcf <- quiet(bcf(y[-index], z[-index], x[-index,], x[-index,], pihat, nburn=n_burn, nsim=n_sim))
  tau_itt <- itt_bcf$tau
  itt <- colMeans(tau_itt)
  
  # Get posterior of treatment effects
  exp <- as.data.frame(cbind(itt, x[-index,]))
  
  ######################################################
  ####  Step 1: Build a CART on the Unit Level CITT ####
  ######################################################
  
  fit.tree <- rpart(itt ~ ., data = exp, maxdepth = max_depth)
  
  ######################################################
  ####  Step 2: Extract the Causal Rules (Nodes)    ####
  ######################################################
  
  rules <- as.numeric(row.names(fit.tree$frame[fit.tree$numresp]))
  
  # Initialize Outputs
  bcfittMat <- as.data.frame(matrix(NA, nrow = length(rules), ncol=4))
  names(bcfittMat) <- c("node", "CCACE", "pvalue", "Weak_IV_test")
  
  # Generate Leaves Indicator
  lvs <- leaves <- numeric(length(rules)) 
  lvs[unique(fit.tree$where)] <- 1
  leaves[rules[lvs==1]] <- 1
  
  ######################################################
  ####  Step 3: Run an IV Regression on each Node   ####
  ######################################################
  
  # Run an IV Regression on the Root
  iv.root <- ivreg(y ~ w | z,  
                   data = inference)
  summary <- summary(iv.root, diagnostics = TRUE)
  iv.effect.root <-  summary$coef[2,1]
  p.value.root <- summary$coef[2,4]
  p.value.weak.iv.root <- summary$diagnostics[1,4]
  proportion.root <- 1
  compliers.root <- length(which(inference$z==inference$w))/nrow(inference)
  itt.root <- iv.effect.root*compliers.root
  
  # Store Results for Root
  bcfittMat[1,] <- c( NA , round(iv.effect.root, 4), round(p.value.root, 4), round(p.value.weak.iv.root, 4), round(proportion.root, 4), round(itt.root, 4), round(compliers.root, 4))
  
  # Initialize New Data
  names(inference) <- paste(names(inference), sep="")
  
  # Run a loop to get the rules (sub-populations)
  for (i in rules[-1]){
    # Create a Vector to Store all the Dimensions of a Rule
    sub <- as.data.frame(matrix(NA, nrow = 1,
                                ncol = nrow(as.data.frame(path.rpart(fit.tree, node=i, print.it = FALSE)))-1))
    quiet(capture.output(for (j in 1:ncol(sub)){
      # Store each Rule as a Sub-population
      sub[,j] <- as.character(print(as.data.frame(path.rpart(fit.tree,node=i,print.it=FALSE))[j+1,1]))
      sub_pop <- noquote(paste(sub , collapse = " & "))
    }))
    
    subset <- with(inference, inference[which(eval(parse(text=sub_pop))),])
    
    # Run the IV Regression
    if (length(unique(subset$w))!= 1 | length(unique(subset$z))!= 1){
      iv.reg <- ivreg(y ~ w | z,  
                      data = subset)
      summary <- summary(iv.reg, diagnostics = TRUE)
      iv.effect <-  summary$coef[2,1]
      p.value <- summary$coef[2,4]
      p.value.weak.iv <- summary$diagnostics[1,4]
      
      ######################################################
      ####   Step 5: Output the Values of each CCACE   ####
      ######################################################
      
      bcfittMat[i,] <- c(sub_pop, round(iv.effect, 4), round(p.value, 4), round(p.value.weak.iv, 4))
    }
    # Delete data
    rm(subset)
  }
  
  # Adjust P.values 
  bcfitt_correction <- cbind(as.data.frame(bcfittMat), leaves)
  adj <- round(p.adjust( as.numeric(bcfitt_correction$pvalue[which(bcfitt_correction$leaves==1)]) ,  paste(adj_method)), 5)
  Adj_pvalue <- rep(NA, length(rules)) 
  Adj_pvalue[which(bcfitt_correction$leaves==1)] <- adj
  
  # Store Results
  bcfittResults <- cbind(as.data.frame(bcfittMat), Adj_pvalue)
  return(bcfittResults)
}


# Don't print in-function messages (source: https://r.789695.n4.nabble.com/Suppressing-output-e-g-from-cat-td859876.html)
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
}
