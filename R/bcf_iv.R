#' @title
#' Bayesian Causal Forest with Instrumental Variables
#'
#' @description
#' Function for discovery and estimation of heterogeneity in the Complier 
#' Average Causal Effect in Instrumental Variable settings.
#'
#' @param y outcome vector (n x 1).
#' @param w vector of treatment received  (n x 1).
#' @param z vector of treatment assigned (i.e., instrumental variable)  (n x 1).
#' @param x covariates matrix (n x p).
#' @param binary Boolean to identify whether the outcome is binary or not, i.e., 
#' continuous/discrete (default: FALSE).
#' @param n_burn burn-in MCMC iterations for Bayesian Causal Forest 
#' (default: 500).
#' @param n_sim iterations to save post burn-in for Bayesian Causal Forest 
#' (default: 500).
#' @param inference_ratio % of observations to be assigned to the inference 
#' subsample (default: 0.5).
#' @param max_depth maximal depth for the generated CART (default: 2).
#' @param cp complexity parameter for the generated CART (default: 0.01).
#' @param minsplit minimum observations needed to perform a binary split in the 
#' tree (default: 10).
#' @param adj_method p-value adjustment method (default: "holm"), other 
#' options are "bonferroni", "hockberg", "hommel", "BH", "BY", "fdr", "none".
#' @param seed random seed for reproducible results (default: 42).
#'
#'
#' @return
#' List with 8 elements:
#'   - a tree structure discovering the heterogeneity in the causal effects, 
#'   - TSLS estimates of the Complier Average Causal Effects (CCACE) within its nodes, 
#'   - p-values for each CCACE, 
#'   - adjusted p-value for each CCACE, 
#'   - Intention-to-Treat within all the nodes, 
#'   - proportion of Compliers within all the nodes, 
#'   - p-value for weak-iv test, 
#'   - proportion of observations in the node.
#'
bcf_iv <- function(y, w, z, x, binary = FALSE, n_burn = 500, n_sim = 500, 
                   inference_ratio = 0.5, max_depth = 2, cp = 0.01, 
                   minsplit = 10, adj_method = "holm", seed = 42) {
  
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


# Don't print in-function messages (source: https://r.789695.n4.nabble.com/Suppressing-output-e-g-from-cat-td859876.html)
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
}
