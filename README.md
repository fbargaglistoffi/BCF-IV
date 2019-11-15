# Code for Bayesian Causal Forest with Instrumental Variable.

In this repository we provide the code for the BCF-IV function and for the application part of the paper "Heterogeneous causal effects with imperfect compliance: a novel Bayesian machine learning approach" by F.J. Bargagli-Stoffi, K. De Witte and G. Gnecco.
The paper can be found [here](https://arxiv.org/pdf/1905.12707.pdf).

# BCF-IV function

The function takes as inputs:

* <tt>`y`</tt>: the outcome variable;
* <tt>`w`</tt>: the reception of the treatment variable (binary);
* <tt>`z`</tt>: the assignment to the treatment variable (binary);
* <tt>`max_depth`</tt>: the maximal depth of the tree generated by the function;
* <tt>`n_burn`</tt>: the number of iterations discarded by the BCF-IV algorithm for the burn-in;
* <tt>`n_sim`</tt>: the number of iterations used by the BCF-IV algorithm  to get the posterior distribution of the estimands;
* <tt>`binary`</tt>: this option should be set to <tt>`TRUE`</tt> when the outcome variable is binary and to <tt>`FALSE`</tt> if the outcome variable is either discrete or continuous.

More details on the R code for the BCF-IV function can be found [here](https://github.com/barstoff/BCF-IV/blob/master/Functions/BCF-IV_in_detail.pdf).

Example usage:

```R
source("bcf-iv.R")

bcf_iv(y, w, z, x, max_depth = 2, n_burn= 2000, n_sim= 2000, binary = TRUE)

mm_bcf_iv(y, w, z, x, max_depth = 2, n_burn= 2000, n_sim= 2000, binary = TRUE)
```
