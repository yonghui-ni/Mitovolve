#' @title Estimation of selection and non-selection models
#' @description Calculate the negative log-likelihood of the model and estimate coefficients for cubic function for log.psi
#' @param mut.read Initial mutant mitochondria read count for each single cell (vector)
#' @param wt.read Initial wildtype mitochondria read count for each single cell (vector)
#' @param nmito.start Mitochondrial copy number of parent cell (scalar)
#' @param nreps Number of generations (scalar)
#' @param selection.reps Number of generations for selection model estimation (vector between (0,nreps))
#' @param selection.mut.start  Number of mutant mitochondrial read count for selection model estimation (vector between (0,nmito.start))
#' @details The function calculates the negative log-likelihood of non-selection models of all numbers of mutant mitochondrial read count and all generations.
#' The function also calculates the negative log-likelihood of selection models with specified numbers of mutant mitochondrial read count and specified generation.
#' Coefficients for cubic function for log.psi (log odds ratio of mutant vs. wildtype read count) are estimated for selection model.
#' @seealso \code{\link{nlogL.null}} \code{\link{local.mle}}
#' @examples
#' set.seed(1234)
#' tot.read = 10                               ## total read count
#' mut.read = sample(tot.read,10,replace = TRUE)
#' wt.read = tot.read - mut.read
#' res = full.mle(mut.read = mut.read,
#'                wt.read = wt.read,
#'                nmito.start = 10,
#'                nreps = 5,
#'
#'               ## specify the generation in selection model estimation
#'               selection.reps = c(3,5),
#'
#'               ## specify the number of mutant mitochondrial read count in selection model estimation
#'               selection.mut.start = c(3,7))
#' @return A list of input (mut.read, wt.read, nmito.start, number of generations) and res.table
#' @export
full.mle=function(mut.read,
                  wt.read,
                  nmito.start,
                  nreps,
                  selection.reps,
                  selection.mut.start)
{
  # compute log-likelihood of no-selection model
  null.res=nlogL.null(mut.read,wt.read,
                      nmito.start,
                      nreps=nreps)
  res=null.res$res.table

  # define log-likelihood of selection model
  res$alt.nlogL=NA
  res$beta0=res$beta1=res$beta2=res$beta3=NA

  # select the selection models to be evaluated
  select= which(res$generation %in% selection.reps
                & res$mutant.start %in% selection.mut.start
                & res$null.nlogL !=Inf)



  k=0
  for (i in c(select))
  {
    k=k+1
    message(paste0("Checking selective model ",k," of ",
                   length(c(select)),": ",date()))
    mle.res=try(local.mle(mut.read=mut.read,
                          wt.read=wt.read,
                          mut.start=res$mutant.start[i],
                          wt.start=res$nmito.start[i]-res$mutant.start[i],
                          nreps=res$generation[i]))
    if (as.character(class(mle.res))!="try-error")
    {
      res$alt.nlogL[i]=mle.res$nlogL
      res$beta0[i]=mle.res$beta[1]
      res$beta1[i]=mle.res$beta[2]
      res$beta2[i]=mle.res$beta[3]
      res$beta3[i]=mle.res$beta[4]
    }

  }

  ## return output
  full.res=list(mut.read=mut.read,
                wt.read=wt.read,
                nmito.start=nmito.start,
                nreps=nreps,
                res.table=res)

  return(full.res)
}


#' @title Estimation for selection and non-selection models of specified number of mutant/wildtype mitochondrial read count
#' @description Calculate negative log-likelihood of the model and estimate cubic function coefficients for log odds ratio in non-central hypergeometric distribution
#' @param mut.read Initial mutant mitochondria read count for each single cell (vector)
#' @param wt.read Initial wildtype mitochondria read count for each single cell (vector)
#' @param mut.start Number of mutant mitochondrial read count (scalar)
#' @param wt.start Number of mutant mitochondrial read count (scalar)
#' @param nreps Number of generations (scalar)
#' @param beta0 Coefficients for cubic function for log.psi
#' @examples
#' set.seed(1234)
#' tot.read = 10                               ## total read count
#' mut.read = sample(tot.read,10,replace = TRUE)
#' wt.read = tot.read - mut.read   ## the sum of mut.start and wt.start should be total read count
#' mle = local.mle(mut.read = mut.read,
#'                 wt.read = wt.read,
#'                 mut.start = 3,
#'                 wt.start = 7,
#'                 nreps = 5,
#'                 beta0 = rep(0,4))
#' @seealso \code{\link{maf.pmf}},\code{\link{dcd.mtx}},\code{\link{read.pmf}},\code{\link{dnchg}},\code{\link{nlogL.maf}}
#' @return A list includes betas estimation for selection model, negative log-likelihood value of selection and non-selection model, model convergence information,and input values
#' @importFrom stats nlminb
#' @export
local.mle=function(mut.read,
                   wt.read,
                   mut.start,
                   wt.start,
                   nreps,
                   beta0=rep(0,4))

{

  nlm.res=nlminb(beta0,nlogL.maf,       # coefficient estimation start from 0 for selection model
                 mut.read=mut.read,
                 wt.read=wt.read,
                 mut.start=mut.start,
                 wt.start=wt.start,
                 nreps=nreps)
  beta=nlm.res$par

  null.nlogL=nlogL.maf(beta0,
                       mut.read=mut.read,
                       wt.read=wt.read,
                       mut.start=mut.start,
                       wt.start=wt.start,
                       nreps=nreps)

  res=list(beta=beta,
           nlogL=nlm.res$objective,
           null.nlogL=null.nlogL,
           wt.read=wt.read,
           mut.read=mut.read,
           wt.start=wt.start,
           mut.start=mut.start,
           nreps=nreps,
           convergence=nlm.res$convergence,
           message=nlm.res$message)

  return(res)
}


#' @title Calculate of negative log-likelihood of the model
#' @description Calculate of negative log-likelihood of the model
#' @param beta Coefficients for cubic function for log.psi
#' @param mut.read Initial mutant mitochondria read count for each single cell (vector)
#' @param wt.read Initial wildtype mitochondria read count for each single cell (vector)
#' @param mut.start Number of mutant mitochondrial read count (scalar)
#' @param wt.start Number of wildtype mitochondrial read count (scalar)
#' @param nreps Number of generations (scalar)
#' @return Negative log-likelihood value
#' @seealso \code{\link{maf.pmf}},\code{\link{dcd.mtx}},\code{\link{read.pmf}},\code{\link{dnchg}}
#' @examples
#' set.seed(1234)
#' tot.read = 10                       ## total read count
#' mut.read = sample(tot.read,10,replace = TRUE)
#' wt.read = tot.read - mut.read    ## the sum of mut.start and wt.start should be total read count
#' ## example for the model without selectio
#' neg_ll = nlogL.maf(beta = rep(0,4),
#'                    mut.read = mut.read,
#'                    wt.read = wt.read,
#'                    mut.start = 3,
#'                    wt.start = 7,
#'                    nreps = 5)
#' @export
nlogL.maf=function(beta,
                   mut.read,
                   wt.read,
                   mut.start,
                   wt.start,
                   nreps)
{
  MAF=maf.pmf(mut.start,wt.start,nreps,beta)
  res=-sum(log(read.pmf(mut.read,wt.read,MAF,nreps)))
  return(res)
}


#' @title Mitochondrial mutations prevalence modeling without selective pressure
#' @description Calculate negative log-likelihood for model of mutant mitochondria read count (from 1 to total copy number) and generation (from 1 to specified generation)
#' @param mut.read Initial mutant mitochondria read count for each single cell (vector)
#' @param wt.read Initial wildtype mitochondria read count for each single cell (vector)
#' @param nmito.start Mitochondrial copy number of parent cell (scalar)
#' @param nreps Number of generations (scalar)
#' @return A list includes all inputs (mut.read,wt.read,nmito.start,nreps) and a data frame with calculated negative log-likelihood value of all models .
#' @details The modeling without selective pressure of mutant mtDNA prevalence in subsequent generations is by iterative convolution of hypergeometric models, in which the mutant read is modeled by hypergeometric sampling for a given number of mutant/wildtype mitochondria read.
#' @seealso \code{\link{maf.pmf}},\code{\link{dcd.mtx}},\code{\link{read.pmf}},\code{\link{dnchg}}
#' @examples
#' tot.read = 20
#' mut.read = sample(tot.read,15,replace = TRUE)
#' wt.read = tot.read - mut.read
#' null.res = nlogL.null(mut.read,wt.read,nmito.start = 100,nreps = 50)
#' @export
nlogL.null=function(mut.read,wt.read,
                    nmito.start=100,
                    nreps=100)

{
  mtx=matrix(NA,nmito.start+1,nreps+1)
  colnames(mtx)=paste0("gen_",0:nreps)
  rownames(mtx)=paste0("mut_",0:nmito.start,"_wt_",nmito.start:0)
  for (i in 0:nmito.start)
  {
    message(paste0("Evaluating models for up to ",nreps," generations emerging ",
                   "from an initiating cell with ",
                   i," mutant mitochondria among ",nmito.start," mitochondria: ",date()))
    PMF=maf.pmf(i,nmito.start-i,nreps)
    for (j in 0:nreps)
      mtx[i+1,j+1]=sum(-log(read.pmf(mut.read,wt.read,PMF,j)))
  }
  vtr=as.vector(mtx)
  res=cbind.data.frame(mutant.start=rep(0:nmito.start,times=nreps+1),
                       nmito.start=nmito.start,
                       generation=rep(0:nreps,each=nmito.start+1),
                       null.nlogL=vtr)


  full.res=list(mut.read=mut.read,
                wt.read=wt.read,
                nmito.start=nmito.start,
                nreps=nreps,
                res.table=res)

  return(full.res)
}

#' @title Compute cumulative probability of mutant allele fraction (MAF)
#' @description Compute cumulative probability of MAF for a population of cells by a given MAF distribution matrix
#' @param maf A vector of values at which to evaluate the CDF
#' @param n.read Number of reads observed in the cells
#' @param maf.mtx Mutant allele fraction matrix from maf.pmf
#' @param nreps Number of generations, default is column of maf.mtx minus 1
#' @return The cumulative probability by given MAF
#' @seealso \code{\link{maf.pmf}}
#' @examples
#' set.seed(1234)
#' tot.read = 20
#' mut.read = sample(tot.read,15,replace = TRUE)
#' wt.read = tot.read - mut.read
#' maf.mtx = maf.pmf(mut.mito = 1,
#'           wt.mito = 19,
#'           nreps = 5,
#'           beta = rep(0,4)) ## example for the model without selection
#' maf = 0:tot.read/tot.read
#' cdf.maf(maf,n.read = tot.read,maf.mtx = maf.mtx)
#' @importFrom stats dhyper approx
#' @export
cdf.maf=function(maf,           # a vector of values at which to evaluate the cdf
                 n.read,        # a vector of the number of reads observed in the cells
                 maf.mtx,       # mutant allele fraction matrix from maf.pmf
                 nreps=NULL)    # number of generations (column of maf.mtx-1)

{
  if(is.null(nreps))
    nreps=ncol(maf.mtx)-1

  ncells=length(n.read)
  nmito=nrow(maf.mtx)-1  # number of mitochondria in the initiating cell
  y=rep(0,nrow(maf.mtx)) # initial value of the CDF
  x=0:(nmito)/nmito

  for (i in 1:ncells)
  {
    pr=rep(NA,n.read[i]+1)
    for (j in 0:n.read[i])
    {
      pr[j+1]=sum(exp(log(maf.mtx[,nreps+1])+dhyper(j,
                                                    0:nmito,
                                                    nmito:0,
                                                    n.read[i],
                                                    log=T)))
    }
    cum.pr=cumsum(pr)
    y=y+cum.pr[floor(x*n.read[i])+1]
  }
  y=y/ncells

  res=approx(x,y,xout=maf)$y
  return(res)

}




#' @title Probability distribution function for mutant read counts
#' @description Probability distribution function for mutant read counts
#' @param mut.read Mutant mitochondria read count for each single cell (vector)
#' @param wt.read Wildtype mitochondria read count for each single cell (vector)
#' @param maf.mtx Mutant allele fraction distribution matrix estimated from maf.pmf function
#' @param nreps Number of generations (scalar)
#' @return A vector of probability of mutant read counts
#' @examples
#' ## example for the model without selection
#' set.seed(1234)
#' maf.mtx = maf.pmf(mut.mito = 1,wt.mito = 19,nreps = 5,beta = rep(0,4))
#' tot.read = 10
#' mut.read = sample(tot.read,10,replace = TRUE)
#' wt.read = tot.read - mut.read
#' read.pmf(mut.read,wt.read,maf.mtx)
#' mut.all = 0:tot.read
#' wt.all = tot.read:0
#' sum(read.pmf(mut.all,wt.all,maf.mtx))
#' @seealso \code{\link{maf.pmf}}
#' @export
read.pmf=function(mut.read,    # mutant mitochondria read count for each single cell
                  wt.read,     # wildtype mitochondria read count for each single cell
                  maf.mtx,     # mutant allele fraction matrix from maf.pmf
                  nreps=NULL)  # number of replications

{
  if(is.null(nreps))
    nreps=ncol(maf.mtx)-1

  ncells=length(wt.read)
  tot.mito=nrow(maf.mtx)-1

  sum.mito=mut.read+wt.read
  if (any(sum.mito>tot.mito))
    stop("mut.read+wt.read cannot be greater than the total number of mtDNA (nrow(maf.mtx)-1) in the model.")

  res=rep(0,ncells)
  for (i in 0:tot.mito)
  {
           # mtDNAs per cell model  + mtDNAs sampled from cell (central HG distribution)
    log.pr=log(maf.mtx[i+1,nreps+1])+dhyper(mut.read,
                                            i,
                                            tot.mito-i,
                                            mut.read+wt.read,
                                            log=T)
    res=res+exp(log.pr)
  }
  return(res)

}




#' @title Compute distribution of mutant allele fraction for a population of cells
#' @description Compute distribution of mutant allele fraction for a population of cells descending over generations from an initiating cell
#' @param mut.mito Number of mutant mitochondrial read count
#' @param wt.mito Number of wildtype mitochondrial read count
#' @param nreps Number of generations (scalar)
#' @param beta Coefficients for cubic function for log.psi
#' @return Matrix of distribution of mutant allele fraction for a population of cells descending over generations
#' @examples
#' ## example for the model without selection
#' maf.mtx = maf.pmf(mut.mito = 1,
#'                   wt.mito = 19,
#'                   nreps = 5,
#'                   beta = rep(0,4))
#' @export
maf.pmf=function(mut.mito,
                 wt.mito,
                 nreps,
                 beta=rep(0,4))

{
  tot.mito=wt.mito+mut.mito                 # initial total number of mitochondrial genomes
  DCD=dcd.mtx(tot.mito,beta)                # daughter cell distribution matrix for this total MG copy number
  res.mtx=matrix(0,tot.mito+1,nreps+1)      # initialize the result matrix
  res.mtx[mut.mito+1,1]=1                   # initially, the cell has mut.mito copies (row mut.mito+1), so this entry has probability 1
  for (i in 2:(nreps+1))                    # loop over cellular replication cycles
  {
    #message(paste0("Working on generation ",i,": ",date()))
    res.mtx[,i]=t(DCD)%*%res.mtx[,i-1]      # generation i+1 is the product of the DCD matrix and generation i
  }
  return(res.mtx)                           # return the result matrix
}



#' @title Compute daughter cell mutant allele fraction distribution
#' @description Compute daughter cell mutant allele fraction distribution
#' @param n.mito Mitochondrial copy number of parent cell
#' @param beta Coefficients for cubic function for log.psi
#' @return Matrix of daughter cell mutant allele fraction distribution
#' @examples
#' maf.dcd = dcd.mtx(n.mito = 20,beta = rep(0,4)) ## example for the model without selection
#' @export
dcd.mtx=function(n.mito,
                 beta=rep(0,4))

{
  res.mtx=matrix(NA,n.mito+1,n.mito+1)         # initialize matrix result
  for (n.mut in 0:n.mito)                      # loop over possible number of mitochondrial genomes
  {
    mtx.row=n.mut+1                            # matrix row is number of mutant MGs +1
    n.wt=n.mito-n.mut                          # number of wildtype genomes
    maf=n.mut/n.mito
    log.psi=beta[1]+beta[2]*maf+beta[3]*maf^2+beta[4]*maf^3
    p1=rep(0,n.mito+1)
    for (j in 0:n.mito)
      p1[j+1]=dnchg(j,
                    2*n.mut,
                    2*n.wt,
                    n.mito,
                    exp(log.psi))
    res.mtx[mtx.row,]=p1
  }
  return(res.mtx)
}


#' @title Probability mass function of non-central hypergeometric distribution
#' @description Evaluate probability of a single point from non-central hypergeometric distribution
#' @param x  The location to evaluate the probability
#' @param n1 The size of group one
#' @param n2 The size of group two
#' @param m1 The size of both two groups
#' @param psi Odds ratio
#' @return The probability at point x
#' @examples
#' prob = dnchg(1,2,38,20,1)
#' @importFrom MCMCpack dnoncenhypergeom
#' @references J. G. Liao and Ori Rosen. 2001. “Fast and Stable Algorithms for Computing and Sampling From the Noncentral Hypergeometric Distribution." The American Statistician. 55: 366-369.
#' @export
dnchg=function(x,n1,n2,m1,psi)
{
  ll <- max(0, m1 - n2)
  uu <- min(n1, m1)
  if ((x>uu)||(x<ll)) return(0)
  res=MCMCpack::dnoncenhypergeom(x,n1,n2,m1,psi)
  return(res)
}




#' @title Compute P-value of likelihood ratio test between two models
#' @description Compute P-value of likelihood ratio (LR) test between two models
#' @param res.tbl res.table from {\link{full.mle}} function, containing negative log-likelihood values of each model
#' @param nreps Number of generations (scalar)
#' @param mutant.start Number of mutant mitochondrial read count (scalar)
#' @param baseline.model A character string specifying the null model, must be one of "best-non-selection"(default) or "best-selection"
#' @param alt.model A character string specifying the alternative model, must be one of "selection" (default) or "non-selection"
#' @return A list includes inputs (specified nreps and mutant.start for models in comparison),the best selection and non-selection model,
#' the baseline model and res.table.pval table with p-values from LR tests
#' @details The function calculates p-value of LR test of two models. E.g., for identifying selective pressure from mitochondrial DNA mutation, We assesses the goodness of fit of
#' selection model and the best non-selection by LR test.
#' @seealso \code{\link{full.mle}}
#' @examples
#' res = P_value(res.tbl = res.tbl,
#'                   nreps = NULL,
#'                   mutant.start = NULL,
#'                   baseline.model =  "best-non-selection",
#'                   alt.model = "selection")
#' @importFrom stats pchisq
#' @export
P_value=function(res.tbl,
                      nreps=NULL,
                      mutant.start=NULL,
                      baseline.model= "best-non-selection",
                      alt.model="selection")
{
  res.tbl<-res.tbl[!is.na(res.tbl$alt.nlogL),]

  suppressWarnings(
    if (is.null(nreps)==FALSE & max(nreps)>max(res.tbl$generation)){
      stop("nrep is out of range of generation in res.tbl")
    }else if (is.null(mutant.start)==FALSE & max(mutant.start)>max(res.tbl$mutant.start)){
      stop('mutant.start is out of range of mutant.start in res.tbl')
    }
  )

  req.clms=c("nmito.start","mutant.start",
             "generation","null.nlogL","alt.nlogL")

  if (any(!(req.clms%in%colnames(res.tbl))))
    stop("res.tbl missing some required columns.")

  # res.tbl subset
  keep.row.nrep<-which(res.tbl$generation %in% nreps)
  keep.row.mut<-which(res.tbl$mutant.start %in% mutant.start)

  if(is.null(nreps)==TRUE & is.null(mutant.start)==TRUE){
    res.tbl=res.tbl
  }else if (is.null(nreps)==TRUE & is.null(mutant.start)==FALSE){
    res.tbl=res.tbl[keep.row.mut,]
  }else if (is.null(nreps)==FALSE & is.null(mutant.start)==TRUE){
    res.tbl=res.tbl[keep.row.nrep,]
  }else{
    res.tbl=res.tbl[intersect(keep.row.nrep,keep.row.mut),]
  }


  # find the best models of non-selection and selection
  tmp<-res.tbl[order(res.tbl$null.nlogL),]
  best.null<-tmp[1,]

  tmp<-res.tbl[order(res.tbl$alt.nlogL),]
  best.alt<-tmp[1,]

  if(baseline.model=="best-non-selection"){
    baseline=best.null
  }else if(baseline.model=="best-selection"){
    baseline=best.alt
  }

  p1<-res.tbl$mutant.start==baseline$mutant.start
  p2<-res.tbl$nmito.start==baseline$nmito.start
  p3<-res.tbl$generation==baseline$generation
  tmp<-cbind(p1,p2,p3)

  # selection model vs. best non-selection
  if(baseline.model=="best-non-selection" & alt.model=="non-selection"){
    baseline=best.null
    df<-c()
    for (i in 1:nrow(res.tbl)){
      df[i]=3-table(tmp[i,])[["TRUE"]]
    }
    res.tbl$df<-df
    cs.stat=2*(res.tbl$null.nlogL-baseline$null.nlogL)
    res.tbl$pval=pchisq(cs.stat,df,lower.tail=F)


  }else if(baseline.model=="best-non-selection" & alt.model=="selection"){
    baseline=best.null
    df<-c()
    for (i in 1:nrow(res.tbl)){
      df[i]=7-table(tmp[i,])[["TRUE"]]
    }
    res.tbl$df<-df
    cs.stat=2*(baseline$null.nlogL-res.tbl$alt.nlogL)
    res.tbl$pval=pchisq(cs.stat,df,lower.tail=F)


  }else if(baseline.model=="best-selection" & alt.model=="selection"){
    baseline=best.alt
    df<-c()
    for (i in 1:nrow(res.tbl)){
      df[i]=7-table(tmp[i,])[["TRUE"]]
    }
    res.tbl$df<-df
    cs.stat=2*(res.tbl$alt.nlogL-baseline$alt.nlogL)
    res.tbl$pval=pchisq(cs.stat,df,lower.tail=F)


  }

  pval.res=list(nreps=nreps,
                mutant.start = mutant.start,
                baseline.model= baseline,
                best.non.selection = best.null,
                best.selection = best.alt,
                res.table.pval = res.tbl)

  return(pval.res)
}



#' @title Obtain the best selection model using clustering method
#' @description Clustering on the selection models and find the best one. The best selection model can be from low-starting mutant allele frequency (MAF)
#' model or high-starting MAF model or both
#' @param res.tbl res.table from {\link{full.mle}} function, containing negative log-likelihood values and estimated coefficients of each model
#' @param K Number of cluster
#' @param show.best A logical indicator. If TRUE then result shows the best model from all clusters. If FALSE then the result shows the best model of each cluster.
#' @return A data frame with the best selection model of low/high-starting MAF and the best no selection model at the last row.
#' @details Use hierarchical clustering method to cluster selection models and find the best one. The best model of each cluster could be with high-staring MAF (MAF>50%) or with low-starting MAF (MAF<=50%).
#' @seealso \code{\link{full.mle}}
#' @examples
#' res.tbl.best = get.best.model(res.tbl = res.tbl,K = 10,show.best = TRUE)
#' @importFrom stats dist hclust cutree
#' @importFrom dplyr %>% group_by slice
#' @importFrom rlang .data
#' @export
get.best.model=function(res.tbl,
                        K=10,
                        show.best=TRUE) # if FALSE,show all best models from each cluster
{

  req.clms=c("nmito.start","mutant.start",
             "generation","null.nlogL","alt.nlogL",
             "beta0","beta1","beta2","beta3")

  if (any(!(req.clms%in%colnames(res.tbl))))
    stop("res.tbl missing some required columns.")


  res.tbl<-res.tbl[,req.clms]

  # clustering on coefficients
  res.tbl<-res.tbl[!is.na(res.tbl$alt.nlogL),]
  betas=res.tbl[,grep("beta",colnames(res.tbl))]
  betas=as.matrix(betas)
  beta.dist=dist(betas)
  beta.hcl=hclust(beta.dist)


  beta.class=cutree(beta.hcl,K)
  res.tbl$beta.grp=beta.class


  # pick best selection model from each clusters
  res.tbl.best<-res.tbl %>%
                group_by(.data$beta.grp) %>%
                slice(which.min(.data$alt.nlogL))

  models=res.tbl.best[order(res.tbl.best$alt.nlogL),]
  models$start=ifelse(models$mutant.start>0.5*models$nmito.start,"high","low")
  res.tbl.best.2<-models %>%
                  group_by(.data$start) %>%
                  slice(which.min(.data$alt.nlogL))

  if (all(models$mutant.start<=0.5*models$nmito.start)==TRUE){
    print("The best model is from modeling with low-starting mutant allele frequency ")
    best.model=models[1,1:(ncol(res.tbl.best.2)-2)]
  } else if(all(models$mutant.start>0.5*models$nmito.start)==TRUE){
    print("The best model is from modeling with high-starting mutant allele frequency")
    best.model=models[1,1:(ncol(res.tbl.best.2)-2)]
  } else{
    print('The best models are from modeling with low-starting MAF and high-starting MAF')
    best.model=res.tbl.best.2[,1:(ncol(res.tbl.best.2)-2)]

  }

  best.null<-res.tbl[which(res.tbl$null.nlogL==min(res.tbl$null.nlogL)),-ncol(res.tbl)]


  if(show.best==FALSE){
    best.model.1<-data.frame(rbind(res.tbl.best[,1:(ncol(res.tbl.best)-1)],best.null))
    best.model.1[nrow(best.model.1),c("beta3","beta2","beta1","beta0")]=0
    best.model.1$nlogL=c(best.model.1$alt.nlogL[1:(nrow(best.model.1)-1)],best.model.1$null.nlogL[nrow(best.model.1)])
    best.model.1<-best.model.1[, !(names(best.model.1) %in% c("null.nlogL","alt.nlogL"))]
    return(best.model.1)
  } else if (show.best==TRUE){
    best.model<-data.frame(rbind(best.model,best.null))
    best.model[nrow(best.model),c("beta3","beta2","beta1","beta0")]=0
    best.model$nlogL=c(best.model$alt.nlogL[1:(nrow(best.model)-1)],best.model$null.nlogL[nrow(best.model)])
    best.model<-best.model[ , !(names(best.model) %in% c("null.nlogL","alt.nlogL"))]
    return(best.model)
  }

}
