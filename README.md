
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Mitovolve

<!-- badges: start -->
<!-- badges: end -->

*Mitovolve* is developed to infer the evolutionary history of somatic
mitochondrial DNA (mtDNA) mutations by singe-cell genomic sequencing
data. The method models the evolutionary changes in the prevalence of a
mtDNA mutation over successive generation of cells. For each cell
sequencing data, we assume that the mutant mtDNA duplications inherits
into daughter cell following a hypergeometric distribution. Then for a
given number of generations, the number of mutant mtDNA is a convolution
of hypergeometrics.In the observed data with multiple cells, a certain
number of mtDNA reads are obtained for each cell. For each number of
mtDNA reads observed per cell, the number of mutant mtDNA reads is
obtained by hypergeometric sampling. In this way, we obtain a model
distribution of the mutant allele fraction (MAF) observed in our data.

For the model with selective pressure, *Mitovolve* uses noncentral
hypergeometric distribution with the log odds ratio of selecting a
mutant mtDNA as a cubic polynomial function
$\theta(m) = \beta_{0}+\beta_{1}m+\beta_{2}m^{2}+\beta_{3}m^{3}$, where
m is the mutant allele fraction.To model the process without
evolutionary pressure, Mitovolve uses a central hypergeometric
distribution obtained by setting *{0}=*{1}=*{2}=*{3}=0 in the above
equation so that there is no preference for a mutant mtDNA over a
wildtype mtDNA (the log odds ratio $\theta(m) = 0$). Positive values for
(m) indicate a selective preference for the mutant mtDNAs; negative
values for (m) indicate a selective preference for the wildtype mtDNAs.

*Mitovolve* uses a likelihood framework to evaluate the fit of a series
of theoretical models. The likelihood ratio test is used to evaluate the
hypothesis that there is no selective pressure by comparing selective
model with the best selection-free model. Finally, all best fitted model
and observed data can be visualized by histogram plot of MAF.

and compares the to the observed distribution of wild-type and mutant
mtDNA reads covering a particular base-pair locus across individual
cells from the single cell sequencing data. Mitovolve computes the log
likelihood (logL) of a given model (i.e. based on starting MAF and
generation number with or without selection) by comparing the
theoretical distributions of mutant and wildtype mtDNA genomes per cell
(as shown in Figure 1G) to the observed distribution. For a
user-specified range of starting MAFs and generation numbers, Mitovolve
uses a likelihood ratio test to compare the best logL among models with
selection to the best logL among models without selection to determine
whether there is statistically significant evidence that the mtDNA
mutation was subject to selective pressure. Mitovolve also uses a series
of likelihood ratio tests to find models with fits that are not
significantly worse than that of the model with the best logL.

## Installation

You can install the *Mitovolve* from github:

``` r
# install.packages("devtools")
#devtools::install_github("yonghui-ni/Mitovolve")
```

## Load package

load the package

``` r
#library(Mitovolve)
```

## Example
