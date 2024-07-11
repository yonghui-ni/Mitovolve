
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
evolutionary pressure, Mitovolve uses a central hypergeometric
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

## Example

Following is the example of the package usage. We have the single-cell
RNA-seq data obtained from a pediatric subject’s leukemic cells
harboring a single tumor enriched mtDNA mutation at the time of
diagnosis. 1,888 leukemia cells are used for analysis, each cell with at
least 1 mutant read. First, we can use full.mle() in the package to do
model estimation by inputing observed reads data (we have example reads
data called “reads” in the package), the total copy number from 0 to
312, the generation from 1 to 45. Total 14398 model without selection
are estimated. Since the selective model estimation of 14398 models
needs long computation time, we specify the MAF and generation to
estimate less selection models in full.mle() step by setting
selection.mut.start = round(seq(0,312,312\*0.01) and selection.reps =
seq(33,43,1). By this step, we have result output from full.mle() and
save as an example called ‘res.tbl’ in our package.

``` r
head(reads)
#>      mut.read wt.read
#> [1,]        6      13
#> [2,]       75      35
#> [3,]        6      67
#> [4,]       85       5
#> [5,]        0       9
#> [6,]        4       2
head(res.tbl[which(res.tbl$generation %in% seq(33,43,1)),])
#>       mutant.start nmito.start generation null.nlogL alt.nlogL     beta3
#> 13588          128         312         43   6399.520  5602.600 0.8907272
#> 13587          127         312         43   6399.593        NA        NA
#> 13589          129         312         43   6400.153        NA        NA
#> 13586          126         312         43   6400.373        NA        NA
#> 13590          130         312         43   6401.489        NA        NA
#> 13585          125         312         43   6401.865  5602.598 0.8976096
#>           beta2    beta1      beta0
#> 13588 -1.573979 1.108664 -0.2496349
#> 13587        NA       NA         NA
#> 13589        NA       NA         NA
#> 13586        NA       NA         NA
#> 13590        NA       NA         NA
#> 13585 -1.575442 1.102004 -0.2436376
```

In the res.tbl, we have negative log-likelihood value of selection-free
model and selection model and coefficient estimates of selection model
cubic function. We can see the best selection-free model:

``` r
res.tbl[which(res.tbl$null.nlogL==min(res.tbl$null.nlogL)),]
#>       mutant.start nmito.start generation null.nlogL alt.nlogL beta3 beta2
#> 14213          127         312         45   6350.531        NA    NA    NA
#>       beta1 beta0
#> 14213    NA    NA
```

Then, we use likelihood ratio test to compares the fit of selection and
no selection models (the best no selection model). Here we specify the
generations we are interested in

``` r
res = P_value(res.tbl = res.tbl,
              baseline.model =  "best-non-selection",
              alt.model = "selection")
head(res$res.table.pval)
#>       mutant.start nmito.start generation null.nlogL alt.nlogL     beta3
#> 13588          128         312         43   6399.520  5602.600 0.8907272
#> 13585          125         312         43   6401.865  5602.598 0.8976096
#> 13591          131         312         43   6403.524  5602.602 0.8838910
#> 13582          122         312         43   6410.637  5602.595 0.9045556
#> 13594          134         312         43   6413.808  5602.604 0.8771039
#> 13275          128         312         42   6425.347  5602.612 0.9285524
#>           beta2    beta1      beta0 df pval
#> 13588 -1.573979 1.108664 -0.2496349  4    0
#> 13585 -1.575442 1.102004 -0.2436376  5    0
#> 13591 -1.572511 1.115270 -0.2556089  5    0
#> 13582 -1.576915 1.095291 -0.2376145  5    0
#> 13594 -1.571043 1.121828 -0.2615628  5    0
#> 13275 -1.635309 1.149901 -0.2588335  5    0
all(res$res.table.pval$pval<0.05)
#> [1] TRUE
```

Since we have 1089 selection models LR test p-values less than 0.05, we
could use clustering method on $\beta$ coefficients of these models to
get the best selection model. We can see that from 10 clusters, the best
selection models are both from low-starting and high-starting MAF. The
third row of res.tbl.best is the best model without selection.

``` r
res.tbl.best = get.best.model(res.tbl = res$res.table.pval,K = 10,show.best = TRUE)
#> [1] "The best models are from modeling with low-starting MAF and high-starting MAF"
res.tbl.best
#>   nmito.start mutant.start generation       beta0    beta1     beta2    beta3
#> 1         312          215         33 -0.61287045 1.882499 -2.267717 1.092491
#> 2         312           22         33  0.02113842 1.082474 -2.164404 1.611037
#> 3         312          128         43  0.00000000 0.000000  0.000000 0.000000
#>      nlogL
#> 1 5602.584
#> 2 5602.132
#> 3 6399.520
```

## Model Visulization

Plot of mtDNA mutation probability distribution over generations

``` r
res.tbl.best$clr = c("red","green","grey")        
plt.mtDNA.gens(res.tbl.best = res.tbl.best, ngen.show=5)
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

Histogram of mtDNA for the mutant allele fraction in the mtDNA mutation
model and beta-smoothed hisogram

``` r
plt.hist(read.data = reads,
           res.tbl.best = res.tbl.best,
           nbin = 1000,
           clr.scheme = "rainbow",
           show.logOR = TRUE,
           plot.title = "Histogram example")
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

``` r
plt.hist.smooth(read.data = reads,
                  res.tbl.best = res.tbl.best,
                  plot.title = "Beta-smoothed historam example")
```

<img src="man/figures/README-unnamed-chunk-9-2.png" width="100%" />

Plot the modeling of mtDNA mutation distribution respect to CDF or log
odds ratio

``` r
plt.mtDNA.model(res.tbl.best = res.tbl.best,
                read.data = reads,
                plot.type = "CDF")
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />
