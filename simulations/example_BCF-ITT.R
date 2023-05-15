# Generate synthetic dataset
dataset <- generate_dataset(n = 1000, 
                            p = 10, 
                            rho = 0, 
                            null = 0, 
                            effect_size = 2, 
                            compliance = 0.75)
y <- dataset[["y"]]
#w <- dataset[["w"]]
z <- dataset[["z"]]
X <- dataset[["X"]]

bcf_itt(y, z, X, 
        n_burn= 2000, 
        n_sim = 2000, 
        inference_ratio = 0.5, 
        binary = FALSE, 
        max_depth = 2)