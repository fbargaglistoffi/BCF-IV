# Bayesian Causal Forest with Instrumental Variable [BayesIV]

In this repository we provide the code for the _BCF-IV_ and _BCF-ITT_ functions of the paper <a href="https://projecteuclid.org/journals/annals-of-applied-statistics/volume-16/issue-3/Heterogeneous-causal-effects-with-imperfect-compliance--A-Bayesian-machine/10.1214/21-AOAS1579.short"> _"Heterogeneous causal effects with imperfect compliance: a Bayesian machine learning approach"_ </a> by F.J. Bargagli-Stoffi, K. De Witte and G. Gnecco published in _The Annals of Applied Statistics_. 

The article has also been covered and summarized in two blog posts on <a href="https://www.r-bloggers.com/2021/12/heterogeneous-treatment-effects-with-instrumental-variables-a-causal-machine-learning-approach/"> _R-bloggers_ </a> and <a href="https://youngstats.github.io/post/2021/12/06/heterogeneous-treatment-effects-with-instrumental-variables-a-causal-machine-learning-approach/"> _YoungStatS_ </a>. Check them out for a coincise summary of the main novelties introduced in the paper.

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

**Attention**: `BayesIV` depends on `bcf` package, which, unfortunately, has just been removed from CRAN due to an un-addressed Issue. In order to run `BayesIV` package, manually install `bcf` package from GitHub following its [installation guideline](https://github.com/jaredsmurray/bcf).

## BCF-IV 

The _bcf-iv_ function discovers and estimates, in an interpretable manner, the effects heterogeneity in settings where the assignment mechanism is irregular (e.g., instrumental variable and fuzzy regression discontinuity scenarios). This function is directly built to discover and estimate the heterogeneity in the Complier Average Treatment Effects (CACE).
The function takes as inputs:

* `y`: the outcome variable;
* `w`: the reception of the treatment variable (binary);
* `z`: the assignment to the treatment variable (binary);
* `x`: the covariate matrix;
* `binary`: `TRUE` if the outcome variable is binary, `FALSE`otherwise (default: 
FALSE);
* `n_burn`: the number of iterations discarded by the BCF-IV algorithm for the 
burn-in (default: 500);
* `n_sim`: the number of iterations used by the BCF-IV algorithm  to get the 
posterior distribution of the estimands (default: 500);
* `inference_ratio`: the ratio of observations to be assigned to the interence 
subsample (default: 0.5);
* `max_depth`: the maximal depth of the tree generated by the function (default: 
2);
* '`cp`: complexity parameter for the generated CART (default: 0.01);
* `minsplit`: minimum observations needed to perform a binary split in the tree 
(default: 10);
* `adj_method`: p-value adjustment method, options are "holm", bonferroni", 
"hockberg", "hommel", "BH", "BY", "fdr", "none" (default: "holm");
* `seed`: random seed for reproducible results (default: 42).

The _bcf_iv_ function returns the discovered sub-population, the conditional complier average treatment effect (CCACE), the p-value for this effect, the p-value for a weak-instrument test, the adjusted p-value, the proportion of compliers, the conditional intention-to-treat effect (CITT) and the proportion of compliers in the node.

## BCF-ITT 

The _bcf-itt_ function discovers the heterogeneity in the intention-to-treat (ITT) and then estimates the effect both for the conditional ITT and the conditional CACE for the discovered subgroups.
The function takes as inputs:

* `y`: the outcome variable;
* `w`: the reception of the treatment variable (binary);
* `z`: the assignment to the treatment variable (binary);
* `x`: the covariate matrix;
* `binary`: `TRUE` if the outcome variable is binary, `FALSE`otherwise (default: 
FALSE);
* `max_depth`: the maximal depth of the generated CART (default: 
2);
* `n_burn`: the number of iterations discarded by the BCF-IV algorithm for the 
burn-in (default: 500);
* `n_sim`: the number of iterations used by the BCF-IV algorithm  to get the 
posterior distribution of the estimands (default: 500);
* `inference_ratio`: the ratio of observations to be assigned to the interence 
subsample (default: 0.5);
* `seed`: random seed for reproducible results (default: 42).

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
bcf_itt(y, w, z, X, 
        n_burn= 2000, 
        n_sim = 2000, 
        inference_ratio = 0.5, 
        binary = FALSE, 
        max_depth = 2)
```

For more exaustive synthetic examples check the folder <a href="https://github.com/fbargaglistoffi/BCF-IV/tree/master/examples">
`examples/`</a>.

## Reference
* Bargagli-Stoffi, F.J., De Witte, K. and Gnecco, G., 2022. Heterogeneous causal effects with imperfect compliance: a Bayesian machine learning approach. The Annals of Applied Statistics, 16(3), pp.1986-2009. </b> [<a href="https://projecteuclid.org/journals/annals-of-applied-statistics/volume-16/issue-3/Heterogeneous-causal-effects-with-imperfect-compliance--A-Bayesian-machine/10.1214/21-AOAS1579.short">paper</a>] [<a href="https://arxiv.org/abs/1905.12707">preprint</a>]

```
@article{bargagli2022heterogeneous,
  title={{Heterogeneous causal effects with imperfect compliance: a Bayesian machine learning approach}},
  author={Bargagli-Stoffi, Falco J and De Witte, Kristof and Gnecco, Giorgio},
  journal={The Annals of Applied Statistics},
  volume={16},
  number={3},
  pages={1986--2009},
  year={2022},
  publisher={Institute of Mathematical Statistics}
}
```
