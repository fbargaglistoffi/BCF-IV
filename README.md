# Code for Bayesian Causal Forest with Instrumental Variable

**UPDATE: New Code out on December 14th, 2020**

In this repository we provide the code for the _BCF-IV_ function and for the application part of the paper _"Heterogeneous causal effects with imperfect compliance: a novel Bayesian machine learning approach"_ by F.J. Bargagli-Stoffi, K. De Witte and G. Gnecco. This function discovers the effect heterogeneity and provides two ways to estimate the heterogeneous causal effect (two-stage least squares and method-of-moments) in scenarios where the assignment mechanism is irregular.
The _BCF-IV_ function directly builds on the Bayesian Causal Forest algorithm by Hahn, Murray and Carvalho.


# BCF-IV function

The function takes as inputs:

* <tt>`y`</tt>: the outcome variable;
* <tt>`w`</tt>: the reception of the treatment variable (binary);
* <tt>`z`</tt>: the assignment to the treatment variable (binary);
* <tt>`max_depth`</tt>: the maximal depth of the tree generated by the function;
* <tt>`n_burn`</tt>: the number of iterations discarded by the BCF-IV algorithm for the burn-in;
* <tt>`n_sim`</tt>: the number of iterations used by the BCF-IV algorithm  to get the posterior distribution of the estimands;
* <tt>`binary`</tt>: this option should be set to <tt>`TRUE`</tt> when the outcome variable is binary and to <tt>`FALSE`</tt> if the outcome variable is either discrete or continuous.

The _mm_bcf_iv_ function returns the discovered sub-population, the conditional complier average treatment effect (CCACE), the p-value for this effect, the p-value for a weak-instrument test, the proportion of compliers, the conditional intention-to-treat effect (CITT) and the proportion of compliers in the node.

More details on the R code for the BCF-IV function can be found [here](https://github.com/barstoff/BCF-IV/blob/master/Functions/BCF-IV_in_detail.pdf).

# Example usage

```R
source("bcf-iv.R")

bcf_iv(y, w, z, x, max_depth = 2, n_burn= 2000, n_sim= 2000, binary = TRUE)

mm_bcf_iv(y, w, z, x, max_depth = 2, n_burn= 2000, n_sim= 2000, binary = TRUE)
```
#### References
* Falco J. Bargagli-Stoffi, Kristof De Witte, Giorgio Gnecco. <b>Heterogeneous causal effects with imperfect compliance: a novel Bayesian machine learning approach.</b> [<a href="https://arxiv.org/abs/1905.12707">link</a>]
* P. Richard Hahn, Jared S. Murray, Carlos Carvalho. <b>Bayesian regression tree models for causal inference: regularization, confounding, and heterogeneous effects.</b> [<a href="https://arxiv.org/abs/1706.09523">link</a>]

