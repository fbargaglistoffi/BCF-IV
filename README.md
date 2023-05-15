# Bayesian Causal Forest with Instrumental Variable [BayesIV]

In this repository we provide the code for the _BCF-IV_ and _BCF-ITT_ functions of the paper <a href="https://projecteuclid.org/journals/annals-of-applied-statistics/volume-16/issue-3/Heterogeneous-causal-effects-with-imperfect-compliance--A-Bayesian-machine/10.1214/21-AOAS1579.short"> _"Heterogeneous causal effects with imperfect compliance: a Bayesian machine learning approach"_ </a> by F.J. Bargagli-Stoffi, K. De Witte and G. Gnecco. 

## Getting Started

Installing the latest developing version: 

```r
library(devtools)
install_github("fbargaglistoffi/BCF-IV", ref="master")
```

Import:

```r
library("BayesIV")
```

## BCF-IV 

The _bcf-iv_ function discovers and estimates, in an interpretable manner, the effects heterogeneity in settings where the assignment mechanism is irregular (e.g., instrumental variable and fuzzy regression discontinuity scenarios). This function is directly built to discover and estimate the heterogeneity in the Complier Average Treatment Effects (CACE).
The function takes as inputs:

* <tt>`y`</tt>: the outcome variable;
* <tt>`w`</tt>: the reception of the treatment variable (binary);
* <tt>`z`</tt>: the assignment to the treatment variable (binary);
* <tt>`n_burn`</tt>: the number of iterations discarded by the BCF-IV algorithm for the burn-in;
* <tt>`n_sim`</tt>: the number of iterations used by the BCF-IV algorithm  to get the posterior distribution of the estimands;
* <tt>`inference_ratio`</tt>: the ratio of observations to be assigned to the interence subsample;
* <tt>`binary`</tt>: this option should be set to <tt>`TRUE`</tt> when the outcome variable is binary and to <tt>`FALSE`</tt> if the outcome variable is either discrete or continuous;
* <tt>`max_depth`</tt>: the maximal depth of the tree generated by the function;
* <tt>`cp`</tt>: complexity parameter for the generated CART (default is 0.01);
* <tt>`minsplit`</tt>: minimum observations needed to perform a binary split in the tree (default is 10);
* <tt>`adj_method`</tt>: p-value adjustment method (default is "holm"), other options are "bonferroni", "hockberg", "hommel", "BH", "BY", "fdr", "none";
* <tt>`seed`</tt> random seed for reproducible results (default is 42).

The _bcf_iv_ function returns the discovered sub-population, the conditional complier average treatment effect (CCACE), the p-value for this effect, the p-value for a weak-instrument test, the adjusted p-value, the proportion of compliers, the conditional intention-to-treat effect (CITT) and the proportion of compliers in the node.

## BCF-ITT 

The _bcf-itt_ function discovers the heterogeneity in the intention-to-treat (ITT) and then estimates the effect both for the conditional ITT and the conditional CACE for the discovered subgroups.
The function takes as inputs:

* <tt>`y`</tt>: the outcome variable;
* <tt>`w`</tt>: the reception of the treatment variable (binary);
* <tt>`z`</tt>: the assignment to the treatment variable (binary);
* <tt>`n_burn`</tt>: the number of iterations discarded by the BCF-IV algorithm for the burn-in;
* <tt>`n_sim`</tt>: the number of iterations used by the BCF-IV algorithm  to get the posterior distribution of the estimands;
* <tt>`inference_ratio`</tt>: the ratio of observations to be assigned to the interence subsample;
* <tt>`binary`</tt>: this option should be set to <tt>`TRUE`</tt> when the outcome variable is binary and to <tt>`FALSE`</tt> if the outcome variable is either discrete or continuous;
* <tt>`max_depth`</tt>: the maximal depth of the tree generated by the function.

The _bcf_itt_ function returns the discovered sub-population, the conditional complier average treatment effect (CCACE), the conditional intention-to-treat (CITT), the p-value for this effect, the p-value for a weak-instrument test,  the adjusted p-value, the proportion of compliers, the conditional intention-to-treat effect (CITT) and the proportion of compliers in the node.

## Examples

```R
# Generate the dataset
dataset <- generate_dataset(n = 1000, 
                            p = 10, 
                            rho = 0, 
                            null = 0, 
                            effect_size = 2, 
                            compliance = 0.75)
y <- dataset[["y"]]
w <- dataset[["w"]]
z <- dataset[["z"]]
X <- dataset[["X"]]

# BCF-IV
bcf_iv(y, w, z, X, 
       n_burn = 2000, 
       n_sim = 2000, 
       inference_ratio = 0.5, 
       binary = FALSE, 
       max_depth = 2, 
       adj_method = "holm")

# BCF-ITT
bcf_itt(y, z, X, 
        n_burn= 2000, 
        n_sim = 2000, 
        inference_ratio = 0.5, 
        binary = FALSE, 
        max_depth = 2)
```

For more exaustive example and synthetic simulations check the folder <a href="https://github.com/fbargaglistoffi/BCF-IV/tree/master/simulations">
`simulation/`</a>.

#### References
* Falco J. Bargagli-Stoffi, Kristof De Witte, Giorgio Gnecco. <b>Heterogeneous causal effects with imperfect compliance: a novel Bayesian machine learning approach.</b> [<a href="https://arxiv.org/abs/1905.12707">link</a>]
* P. Richard Hahn, Jared S. Murray, Carlos Carvalho. <b>Bayesian regression tree models for causal inference: regularization, confounding, and heterogeneous effects.</b> [<a href="https://arxiv.org/abs/1706.09523">link</a>]

