
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Mitovolve

<!-- badges: start -->
<!-- badges: end -->

*Mitovolve* is developed to infer the evolutionary history of somatic
mitochondrial DNA (mtDNA) mutations of singe-cell genomic sequencing
data. The method models the evolutionary changes in the prevalence of a
mtDNA mutation over successive generation of cells. For each cell, we
assume that the mutant mtDNA duplications inherits into daughter cell
following a hypergeometric distribution. Then for a given number of
generations, the number of mutant mtDNA is a convolution of
hypergeometrics. In the observed data with multiple cells, a certain
number of mtDNA reads are obtained for each cell. For each number of
mtDNA reads observed per cell, the number of mutant mtDNA reads is
obtained by hypergeometric sampling. In this way, we obtain a model
distribution of the mutant allele fraction (MAF) observed in our data.

For the model with selective pressure, *Mitovolve* uses noncentral
hypergeometric distribution with the log odds ratio of selecting a
mutant mtDNA as a cubic polynomial function
$\theta(m) = \beta_{0}+\beta_{1}m+\beta_{2}m^{2}+\beta_{3}m^{3}$, where
m is the mutant allele fraction. To model the process without
evolutionary pressure, *Mitovolve* uses a central hypergeometric
distribution obtained by setting
$\beta_{0}$=$\beta_{1}$=$\beta_{2}$=$\beta_{3}$=0 in the above equation
so that there is no preference for a mutant mtDNA over a wildtype mtDNA
(the log odds ratio $\theta(m) = 0$). Positive values for $\theta(m)$
indicate a selective preference for the mutant mtDNAs; negative values
for $\theta(m)$ indicate a selective preference for the wildtype mtDNAs.

*Mitovolve* uses a likelihood framework to evaluate the fit of a series
of probabilistic models. The likelihood ratio test is used to evaluate
the hypothesis that there is no selective pressure by comparing
selective model with the best selection-free model. Finally, all best
fitted model and observed data can be visualized by histogram plot of
MAF.

## Installation

You can install the *Mitovolve* from github:

``` r
#install.packages("devtools")
devtools::install_github("yonghui-ni/Mitovolve")

```

## Load package

``` r
library(Mitovolve)
```
