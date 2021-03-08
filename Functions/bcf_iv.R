bcf_iv <- function(y, w, z, x, max_depth, n_burn, n_sim, inference_ratio, binary = FALSE) {
  
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
                      maxdepth = max_depth)
    
    ######################################################
    ####  Step 3: Extract the Causal Rules (Nodes)    ####
    ######################################################
    
    rules <- as.numeric(row.names(fit.tree$frame[fit.tree$numresp]))
    
    # Exclude the Root
    rules <- rules[-1] 
    
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
    compliers.root <- length(which(z==w))/nrow(inference)
    itt.root <- iv.effect.root*compliers.root
    
    # Print Root Results
    cat(paste("The effect on the overall sample is", round(iv.effect.root, 4)),"\n")
    cat(paste("P-value", p.value.root),"\n")
    cat(paste("P-value Weak-Instrument Test", p.value.weak.iv.root),"\n")
    cat(paste("Proportion of observations in the node: ", "1.00"),"\n")
    cat(paste("Intention-to-treat effect: ", itt.root),"\n")
    cat(paste("Proportion of compliers in the node: ", compliers.root),"\n")
    
    # Initialize New Data
    names(inference) <- paste(names(inference), sep="")
    
    # Run a loop to get the rules (sub-populations)
    for (i in rules){
      # Create a Vector to Store all the Dimensions of a Rule
      sub <- as.data.frame(matrix(NA, nrow = 1,
                                  ncol = nrow(as.data.frame(path.rpart(fit.tree,node=i)))-1))
      quiet(capture.output(for (j in 1:ncol(sub)){
        # Store each Rule as a Sub-population
        sub[,j] <- as.character(print(as.data.frame(path.rpart(fit.tree,node=i))[j+1,1]))
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
        
        cat(paste("The conditional effect on the subpopulation is", round(iv.effect, 4)),"\n")
        cat(paste("P-value", p.value),"\n")
        cat(paste("P-value Weak-Instrument Test", p.value.weak.iv),"\n")
        cat(paste("Proportion of observations in the node: ", proportion.node),"\n")
        cat(paste("Intention-to-treat effect: ", itt),"\n")
        cat(paste("Proportion of compliers in the node: ", compliers),"\n")
        
      }
      
      # Delete data
      rm(subset)
      
    }
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
                      maxdepth = max_depth)
    
    ######################################################
    ####  Step 3: Extract the Causal Rules (Nodes)    ####
    ######################################################
    
    rules <- as.numeric(row.names(fit.tree$frame[fit.tree$numresp]))
    
    # Exclude the Root
    rules <- rules[-1] 
    
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
    compliers.root <- length(which(z==w))/nrow(inference)
    itt.root <- iv.effect.root*compliers.root
    
    # Print Root Results
    cat(paste("The effect on the overall sample is", round(iv.effect.root, 4)),"\n")
    cat(paste("P-value", p.value.root),"\n")
    cat(paste("P-value Weak-Instrument Test", p.value.weak.iv.root),"\n")
    cat(paste("Proportion of observations in the node: ", "1.00"),"\n")
    cat(paste("Intention-to-treat effect: ", itt.root),"\n")
    cat(paste("Proportion of compliers in the node: ", compliers.root),"\n")
    
    # Initialize New Data
    names(inference) <- paste(names(inference), sep="")
    
    # Run a loop to get the rules (sub-populations)
    for (i in rules){
      # Create a Vector to Store all the Dimensions of a Rule
      sub <- as.data.frame(matrix(NA, nrow = 1,
                                  ncol = nrow(as.data.frame(path.rpart(fit.tree,node=i)))-1))
      quiet(capture.output(for (j in 1:ncol(sub)){
        # Store each Rule as a Sub-population
        sub[,j] <- as.character(print(as.data.frame(path.rpart(fit.tree,node=i))[j+1,1]))
        sub_pop <- noquote(paste(sub , collapse = " & "))
      }))
      
      
      subset <- with(inference, inference[which(eval(parse(text=sub_pop))),])
      
      # Run IV Regression
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
        
        cat(paste("The conditional effect on the subpopulation is", round(iv.effect, 4)),"\n")
        cat(paste("P-value", p.value),"\n")
        cat(paste("P-value Weak-Instrument Test", p.value.weak.iv),"\n")
        cat(paste("Proportion of observations in the node: ", proportion.node),"\n")
        cat(paste("Intention-to-treat effect: ", itt),"\n")
        cat(paste("Proportion of compliers in the node: ", compliers),"\n")
        
      }
      
      # Delete data
      rm(subset)
      
    }
  }
  
}

bcf_itt <- function(y, w, z, x, max_depth, n_burn, n_sim, inference_ratio, binary = FALSE) {
  
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
  
  
  ######################################################
  ####     Continuous and Discrete Outcomes         ####
  ######################################################
  
  if (binary == FALSE){
    
    # Perform the Bayesian Causal Forest for the ITT
    bcf_fit <- quiet(bcf(y[-index], z[-index], x[-index,], x[-index,], pihat, nburn=n_burn, nsim=n_sim))
    
    # Compute Posterior
    tau_post <- bcf_fit$tau
    tauhat <- colMeans(tau_post)
    exp <- as.data.frame(cbind(tauhat, x[-index,]))
    
    ######################################################
    ####  Step 2: Build a CART on the Unit Level CITT ####
    ######################################################
    
    fit.tree <- rpart(tauhat ~ .,
                      data = exp,
                      maxdepth = max_depth)
    
    ######################################################
    ####  Step 3: Extract the Causal Rules (Nodes)    ####
    ######################################################
    
    rules <- as.numeric(row.names(fit.tree$frame[fit.tree$numresp]))
    
    # Exclude the Root
    rules <- rules[-1] 
    
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
    
    # Print Root Results
    cat(paste("The effect on the overall sample is", round(iv.effect.root, 4)),"\n")
    cat(paste("P-value", p.value.root),"\n")
    cat(paste("P-value Weak-Instrument Test", p.value.weak.iv.root),"\n")
    cat(paste("Proportion of observations in the node: ", "1.00"),"\n")
    
    # Initialize New Data
    names(inference) <- paste(names(inference), sep="")
    
    # Run a loop to get the rules (sub-populations)
    for (i in rules){
      # Create a Vector to Store all the Dimensions of a Rule
      sub <- as.data.frame(matrix(NA, nrow = 1,
                                  ncol = nrow(as.data.frame(path.rpart(fit.tree,node=i)))-1))
      quiet(capture.output(for (j in 1:ncol(sub)){
        # Store each Rule as a Sub-population
        sub[,j] <- as.character(print(as.data.frame(path.rpart(fit.tree,node=i))[j+1,1]))
        sub_pop <- noquote(paste(sub, collapse = " & "))
      }))
      
      
      subset <- with(inference, inference[which(eval(parse(text=sub_pop))),])
      
      # Run the IV Regression
      if (length(unique(subset$w))!= 1 | length(unique(subset$z))!= 1){
        iv.reg <- ivreg(y ~ w | z ,  
                        data = subset)
        summary <- summary(iv.reg, diagnostics = TRUE)
        iv.effect <-  summary$coef[2,1]
        p.value <- summary$coef[2,4]
        p.value.weak.iv <- summary$diagnostics[1,4]
        
        # Proportion of observations in the node
        proportion.node <- nrow(subset)/nrow(inference)
        
        ######################################################
        ####   Step 5: Output the Values of each CCACE   ####
        ######################################################
        
        cat(paste("The conditional effect on the subpopulation is", round(iv.effect, 4)),"\n")
        cat(paste("P-value", p.value),"\n")
        cat(paste("P-value Weak-Instrument Test", p.value.weak.iv),"\n")
        cat(paste("Proportion of observations in the node: ", proportion.node),"\n")
        
      }
      
      # Delete data
      rm(subset)
      
    }
  }
  
  ######################################################
  ####           Binary Outcomes (Probit)           ####
  ######################################################
  
  if (binary == TRUE){
    
    # Fit Binary BART (from BARTcause)
    fit_itt <- quiet(bartc(y[-index], z[-index], x[-index,], n.samples = n_sim, n.burn = n_burn, n.chains = 2L))
    
    # Get Posterior of Treatment Effects
    ites <- extract(fit_itt, type = "ite")
    tauhat <- apply(ites, 2, mean)
    exp <- as.data.frame(cbind(tauhat, x[-index,]))
    
    ######################################################
    ####  Step 2: Build a CART on the Unit Level CITT ####
    ######################################################
    
    fit.tree <- rpart(tauhat ~.,
                      data = exp,
                      maxdepth = max_depth)
    
    ######################################################
    ####  Step 3: Extract the Causal Rules (Nodes)    ####
    ######################################################
    
    rules <- as.numeric(row.names(fit.tree$frame[fit.tree$numresp]))
    
    # Exclude the Root
    rules <- rules[-1] 
    
    
    ######################################################
    ####  Step 4: Run an IV Regression on each Node   ####
    ######################################################
    
    # Run an IV Regression on the Root
    iv.root <- ivreg(y ~ w| z ,  
                     data = inference)
    summary <- summary(iv.root, diagnostics = TRUE)
    iv.effect.root <-  summary$coef[2,1]
    p.value.root <- summary$coef[2,4]
    p.value.weak.iv.root <- summary$diagnostics[1,4]
    
    # Print Root Results
    cat(paste("The effect on the overall sample is", round(iv.effect.root, 4)),"\n")
    cat(paste("P-value", p.value.root),"\n")
    cat(paste("P-value Weak-Instrument Test", p.value.weak.iv.root),"\n")
    cat(paste("Proportion of observations in the node: ", "1.00"),"\n")
    
    # Initialize New Data
    names(inference) <- paste(names(inference), sep="")
    
    # Run a loop to get the rules (sub-populations)
    for (i in rules){
      # Create a Vector to Store all the Dimensions of a Rule
      sub <- as.data.frame(matrix(NA, nrow = 1,
                                  ncol = nrow(as.data.frame(path.rpart(fit.tree,node=i)))-1))
       quiet(capture.output(for (j in 1:ncol(sub)){
        # Store each Rule as a Sub-population
        sub[,j] <- as.character(print(as.data.frame(path.rpart(fit.tree,node=i))[j+1,1]))
        sub_pop <- noquote(paste(sub , collapse = " & "))
      }))
      
      subset <- with(inference, inference[which(eval(parse(text=sub_pop))),])
      
      # Run the IV Regression
      if (length(unique(subset$w))!= 1 | length(unique(subset$z))!= 1){
        
        iv.reg <- ivreg(y ~ w | z ,  
                        data = subset)
        summary <- summary(iv.reg, diagnostics = TRUE)
        iv.effect <-  summary$coef[2,1]
        p.value <- summary$coef[2,4]
        p.value.weak.iv <- summary$diagnostics[1,4]
        
        # Proportion of observations in the node
        proportion.node <- nrow(subset)/nrow(inference)
        
        ######################################################
        ####   Step 5: Output the Values of each CCACE   ####
        ######################################################
        
        cat(paste("The conditional effect on the subpopulation is", round(iv.effect, 4)),"\n")
        cat(paste("P-value", p.value),"\n")
        cat(paste("P-value Weak-Instrument Test", p.value.weak.iv),"\n")
        cat(paste("Proportion of observations in the node: ", proportion.node),"\n")
        
      }
      
      # Delete data
      rm(subset)
    }
  }
  
}

# Don't print in-function messages (source: https://r.789695.n4.nabble.com/Suppressing-output-e-g-from-cat-td859876.html)
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
}
