#' @title Plot the modeling of mtDNA mutation distribution
#' @description Plot the modeling of mtDNA mutation distribution respect to CDF or log odds ratio
#' @param res.tbl.best A data frame from get.best.model() output with columns: nmito.start, mutant.start, generation, beta0-beta3, clr(optional), each row represents a model
#' @param read.data A matrix with columns: mut.read and wt.read, representing mutant and wildtype mitochondria read count for each single cell
#' @param plot.type Options must be one of "CDF" or "OR"
#' @param read.clr Color for observed data in CDF plot
#' @param brks Breaks for bar plot of observed reads
#' @param clr.scheme Default color scheme for models if no clr specified in res.tbl.best
#' @return Plot the selection and no selection modeling of mtDNA mutation distribution respect to CDF or log odds ratio
#' @seealso \code{\link{get.best.model}} to get the best selection model and best no selection model
#' @examples
#' res.tbl.best = get.best.model(res.tbl = res.tbl,K = 10,show.best = TRUE)
#' plt.mtDNA.model(res.tbl.best = res.tbl.best,
#'                  read.data = reads,
#'                  plot.type = "CDF")
#' @importFrom grDevices colorRamp rgb rainbow
#' @importFrom graphics layout screen par lines segments text
#' @export
plt.mtDNA.model=function(res.tbl.best,
                          read.data,
                          plot.type="CDF",
                          read.clr="black",
                          brks=0:25/25,
                          clr.scheme="rainbow")

{
  # check input
  req.clms=c("nmito.start","mutant.start",
             "generation",
             "beta0","beta1","beta2","beta3","nlogL")

  if (any(!(req.clms%in%colnames(res.tbl.best))))
    stop("res.tbl.best missing some required columns.")

  if (is.null(res.tbl.best$clr))
    res.tbl.best$clr=rainbow(nrow(res.tbl.best))


  # process observed data
  nread=read.data[,"mut.read"]+read.data[,"wt.read"]
  obs.maf=read.data[,"mut.read"]/nread
  obs.tbl=table(obs.maf)
  obs.cdf=cumsum(obs.tbl)/sum(obs.tbl)
  obs.pmf=obs.tbl/sum(obs.tbl)
  obs.maf=as.numeric(names(obs.cdf))


  # initialize objects to store model results
  mdl.maf=vector("list",nrow(res.tbl.best)) # store full MAF matrix
  plt.maf=vector("list",nrow(res.tbl.best)) # horizontal MAF axis values
  cdf=vector("list",nrow(res.tbl.best)) # CDF for plotting

  # compute plotting data for each model
  for (i in 1:nrow(res.tbl.best))
  {
    mdl.maf[[i]]=maf.pmf(mut.mito=res.tbl.best$mutant.start[i],
                         wt.mito=res.tbl.best$nmito.start[i]-res.tbl.best$mutant.start[i],
                         nreps=res.tbl.best$generation[i],
                         beta=unlist(res.tbl.best[i,c("beta0","beta1","beta2","beta3")]))
    maf=0:(res.tbl.best$nmito.start[i])
    plt.maf[[i]]=maf/res.tbl.best$nmito.start[i]
    cdf[[i]]=cdf.maf(plt.maf[[i]],nread,mdl.maf[[i]])
  }

  # Produce CDF plot
  if (any(plot.type%in%"CDF"))
  {
    layout(mat=matrix(1:2,2,1),
           heights=c(0.70,0.30))
    screen(1)
    par(mar=c(5,5,1,2))
    plot(c(0,1),c(0,1),
         ylab="Cumulative Distribution",
         xlab="Mutant Allele Fraction",
         type="n",las=1)
    for (i in 1:nrow(res.tbl.best))
      lines(plt.maf[[i]],cdf[[i]],
            col=res.tbl.best$clr[i],lwd=2)

    lines(obs.maf,obs.cdf,lwd=2,col=read.clr)

    screen(2)
    par(mar=c(1,5,1,2))
    plot(c(0,1),c(0,1),axes=F,
         xlab="",ylab="",type="n")
    segments(0.05,((nrow(res.tbl.best)+1):1)/(nrow(res.tbl.best)+3),
             0.10,col=c(res.tbl.best$clr,read.clr),lwd=2)
    text(0.20,((nrow(res.tbl.best)+2):1)/(nrow(res.tbl.best)+3),
         c("m0",paste0(res.tbl.best$mutant.start,"/",res.tbl.best$nmito.start),"data"),
         cex=0.75)

    text(0.30,((nrow(res.tbl.best)+2):1)/(nrow(res.tbl.best)+3),
         c("g",paste0(res.tbl.best$generation),"data"),
         cex=0.75)

    text(0.40,((nrow(res.tbl.best)+2):1)/(nrow(res.tbl.best)+3),
         c("b0",paste0(round(res.tbl.best$beta0,digits=2)),"data"),
         cex=0.75)

    text(0.50,((nrow(res.tbl.best)+2):1)/(nrow(res.tbl.best)+3),
         c("b1",paste0(round(res.tbl.best$beta1,digits=2)),"data"),
         cex=0.75)

    text(0.60,((nrow(res.tbl.best)+2):1)/(nrow(res.tbl.best)+3),
         c("b2",paste0(round(res.tbl.best$beta2,digits=2)),"data"),
         cex=0.75)

    text(0.70,((nrow(res.tbl.best)+2):1)/(nrow(res.tbl.best)+3),
         c("b3",paste0(round(res.tbl.best$beta3,digits=2)),"data"),
         cex=0.75)

    text(0.80,((nrow(res.tbl.best)+2):1)/(nrow(res.tbl.best)+3),
         c("nlogL",paste0(round(res.tbl.best$nlogL,digits=2)),"data"),
         cex=0.75)

  }

  if (any(plot.type%in%"OR"))
  {
    u=0:1000/1000
    y=matrix(NA,length(u),nrow(res.tbl.best))
    for (i in 1:nrow(res.tbl.best))
      y[,i]=res.tbl.best$beta0[i]+res.tbl.best$beta1[i]*u+
      res.tbl.best$beta2[i]*u^2+res.tbl.best$beta3[i]*u^3
    y.rng=range(as.vector(y))
    layout(mat=matrix(1:2,2,1),
           heights=c(0.70,0.30))

    screen(1)
    par(mar=c(5,5,1,2))
    plot(c(0,1),y.rng,
         ylab="Mutant Selection \n Log Odds Ratio",
         xlab="Mutant Allele Fraction",
         type="n",las=1)
    for (i in 1:nrow(res.tbl.best))
      lines(u,y[,i],col=res.tbl.best$clr[i],lwd=2)

    screen(2)
    par(mar=c(1,5,1,2))
    plot(c(0,1),c(0,1),axes=F,
         xlab="",ylab="",type="n")

    segments(0.05,((nrow(res.tbl.best)):1)/(nrow(res.tbl.best)+2),
             0.10,col=c(res.tbl.best$clr,read.clr),lwd=2)

    text(0.20,((nrow(res.tbl.best)+1):1)/(nrow(res.tbl.best)+2),
         c("m0",paste0(res.tbl.best$mutant.start,"/",res.tbl.best$nmito.start)),
         cex=0.75)

    text(0.30,((nrow(res.tbl.best)+1):1)/(nrow(res.tbl.best)+2),
         c("g",paste0(res.tbl.best$generation)),
         cex=0.75)

    text(0.40,((nrow(res.tbl.best)+1):1)/(nrow(res.tbl.best)+2),
         c("b0",paste0(round(res.tbl.best$beta0,digits=2))),
         cex=0.75)

    text(0.50,((nrow(res.tbl.best)+1):1)/(nrow(res.tbl.best)+2),
         c("b1",paste0(round(res.tbl.best$beta1,digits=2))),
         cex=0.75)

    text(0.60,((nrow(res.tbl.best)+1):1)/(nrow(res.tbl.best)+2),
         c("b2",paste0(round(res.tbl.best$beta2,digits=2))),
         cex=0.75)

    text(0.70,((nrow(res.tbl.best)+1):1)/(nrow(res.tbl.best)+2),
         c("b3",paste0(round(res.tbl.best$beta3,digits=2))),
         cex=0.75)

    text(0.80,((nrow(res.tbl.best)+1):1)/(nrow(res.tbl.best)+2),
         c("nlogL",paste0(round(res.tbl.best$nlogL,digits=2))),
         cex=0.75)

  }
}




#' @title Plot the mtDNA mutation probability distribution over generations
#' @description Plot the mtDNA mutation probability distribution over generations
#' @param res.tbl.best A data frame from get.best.model() output with columns: nmito.start, mutant.start, generation, beta0-beta3, clr(optional), each row represents a model
#' @param ngen.show The number of segment of generations to show
#' @param clr.scheme Default color scheme for models
#' @return Plot of mtDNA mutation probability distribution over generations
#' @seealso \code{\link{get.best.model}} to get the best selection model and best no selection model
#' @examples
#' res.tbl.best = get.best.model(res.tbl,K = 10,show.best = TRUE)
#' res.tbl.best$clr = c("red","green","grey")
#' plt.mtDNA.gens(res.tbl.best = res.tbl.best, ngen.show=5)
#' @importFrom grDevices colorRamp rgb rainbow
#' @importFrom graphics layout lines mtext par screen segments text
#' @export
plt.mtDNA.gens=function(res.tbl.best,
                         ngen.show=5,
                         clr.scheme="rainbow")
{

  # check input
  req.clms=c("nmito.start","mutant.start",
             "generation",
             "beta0","beta1","beta2","beta3","nlogL")

  if (any(!(req.clms%in%colnames(res.tbl.best))))
    stop("res.tbl.best missing some required columns.")

  if (is.null(res.tbl.best$clr))
    res.tbl.best$clr=rainbow(nrow(res.tbl.best))


  # Generate Null Model Progression
  best.null<-which(rowSums(res.tbl.best[,grep("beta",colnames(res.tbl.best))])==0)

  null.maf=maf.pmf(res.tbl.best$mutant.start[best.null],
                   res.tbl.best$nmito.start[best.null]-res.tbl.best$mutant.start[best.null],
                   res.tbl.best$generation[best.null])

  null.gens=round(seq(from=0,to=res.tbl.best$generation[best.null],length=ngen.show))
  null.mdl=res.tbl.best[best.null,]

  par(mfrow=c(2,2))
  # Cellular mtDNA copies in Selection-Free Model
  plot(c(0,nrow(null.maf)-1),-c(0,length(null.gens)),
       xlab=paste0("Cellular Mutant mtDNA (out of ",null.mdl$nmito.start,")"),
       ylab="",yaxt="n",type="n",
       main="Best selection-free Model",
       sub=paste0("m = ",null.mdl$mutant.start,"/",
                  null.mdl$nmito.start,"; g = ",
                  null.mdl$generation,"; -logL = ",
                  format(null.mdl$nlogL,digits=2)),
       cex.sub=0.75,cex.main=0.75,cex.lab=0.75,cex.axis=0.75)
  mtext(null.gens,
        at=-(1:length(null.gens))+0.5,
        side=2,las=1,adj=1.5,cex=0.75)
  mtext("Generation",
        at=0.25,side=2,las=1,cex=0.75)
  for (i in 1:length(null.gens))
  {
    h=null.maf[,null.gens[i]+1]/(1.1*max(null.maf[,null.gens[i]+1]))
    segments((1:nrow(null.maf))-1,
             (-i+h),
             (1:nrow(null.maf))-1,
             -i,col=res.tbl.best$clr[best.null])
  }


  # Generate Alternative Model Progression
  best.alt<-which(rowSums(res.tbl.best[,grep("beta",colnames(res.tbl.best))])!=0)
  for (i in best.alt){

    alt.beta=unlist(res.tbl.best[i,c("beta0","beta1","beta2","beta3")])

    alt.maf=maf.pmf(res.tbl.best$mutant.start[i],
                    res.tbl.best$nmito.start[i]-res.tbl.best$mutant.start[i],
                    res.tbl.best$generation[i],
                    alt.beta)

    alt.gens=round(seq(from=0,to=res.tbl.best$generation[i],length=ngen.show))
    alt.mdl=res.tbl.best[i,]


    plot(c(0,nrow(alt.maf)-1),-c(0,length(alt.gens)),
         xlab=paste0("Cellular Mutant mtDNA (out of ",alt.mdl$nmito.start,")"),
         ylab="",yaxt="n",type="n",
         main=paste0("Best selection Model"),
         sub=paste0("m = ",alt.mdl$mutant.start,"/",
                    alt.mdl$nmito.start,"; g = ",
                    alt.mdl$generation,"; -logL = ",
                    format(alt.mdl$nlogL,digits=2)),
         cex.sub=0.75,cex.main=0.75,cex.lab=0.75,cex.axis=0.75)
    mtext(alt.gens,
          at=-(1:length(alt.gens))+0.5,
          side=2,las=1,adj=1.5,cex=0.75)
    mtext("Generation",
          at=0.25,side=2,las=1,cex=0.75)
    for (i in 1:length(alt.gens))
    {
      h=alt.maf[,alt.gens[i]+1]/(1.1*max(alt.maf[,alt.gens[i]+1]))
      segments((1:nrow(alt.maf))-1,
               (-i+h),
               (1:nrow(alt.maf))-1,
               -i,col=alt.mdl$clr)
    }
  }

}


#' @title Plot histogram of mtDNA for the mutant allele fraction in the mtDNA mutation model
#' @description Plot the histogram of the selection model, selection-free model and observed reads data
#' @param read.data A matrix with columns: mut.read and wt.read, representing observed mutant and wildtype mitochondria read count for each single cell
#' @param res.tbl.best A data frame from get.best.model() output with columns: nmito.start, mutant.start, generation, beta0-beta3, clr(optional), each row represents a model
#' @param nbin Number of bin for histogram
#' @param clr.scheme Default color scheme for models
#' @param show.logOR A logical indicator. If TRUE, the plot showes log odds ratio of mutant mtDNA over wildtype mtDNA of the best selection model
#' @param plot.title A character string, default is NULL
#' @return Probability mass histogram for the mutant allele fraction in the observed read count of selection model and no selection model
#' @seealso \code{\link{get.best.model}}, \code{\link{sc.hist}}, \code{\link{hg.hist}}
#' @examples
#' res.tbl.best = get.best.model(res.tbl= res.tbl,K = 10,show.best = TRUE)
#' res.tbl.best$clr = c("red","green","grey")
#' plt.hist(read.data = reads,
#'           res.tbl.best = res.tbl.best,
#'           nbin = 1000,
#'           clr.scheme = "rainbow",
#'           show.logOR = TRUE,
#'           plot.title = NULL)
#' @importFrom patchwork plot_layout
#' @importFrom grDevices colorRamp rgb rainbow
#' @importFrom ggplot2 ggplot geom_line scale_color_manual labs theme geom_segment geom_text scale_color_identity theme_void aes xlim element_blank element_line element_rect element_text margin
#' @importFrom rlang .data
#' @export
plt.hist=function(read.data,res.tbl.best,nbin=1000,clr.scheme="rainbow",show.logOR=TRUE,plot.title=NULL)
{

  # check input
  req.clms=c("nmito.start","mutant.start",
             "generation",
             "beta0","beta1","beta2","beta3")

  if (any(!(req.clms%in%colnames(res.tbl.best))))
    stop("res.tbl.best missing some required columns.")

  if (nrow(res.tbl.best)<2)
    stop("res.tbl.best must have at least one selection model and one selection-free model")

  if (is.null(res.tbl.best$clr))
    res.tbl.best$clr=rainbow(nrow(res.tbl.best))


  obs.hist=sc.hist(read.data,nbin=nbin)  # read count matrix, one row per cell, first column with mutant counts, second with reference counts


  # Generate Null Model Progression
  best.null<-which(rowSums(res.tbl.best[,grep("beta",colnames(res.tbl.best))])==0)
  null.hist<-hg.hist(read.data,                # input data
                     res.tbl.best[best.null,], # list with mut0, wt0, gens, b0, b1, b2, b3
                     nbin=nbin)

  best.alt <- which(rowSums(res.tbl.best[, grep("beta", colnames(res.tbl.best))]) != 0)
  alt.hist<-list()
  for (i in best.alt){
    alt.hist[[i]]<-hg.hist(read.data,        # input data
                           res.tbl.best[i,], # list with mut0, wt0, gens, b0, b1, b2, b3
                           nbin=nbin)
  }

  # Combine histogram data into data frames for ggplot
  obs_df <- data.frame(x = obs.hist$x, y = obs.hist$y, Type = "Observed")
  null_df <- data.frame(x = null.hist$x, y = null.hist$y, Type = "No selection")
  alt_df <- lapply(1:nrow(res.tbl.best[best.alt,]), function(i) data.frame(x = alt.hist[[i]]$x, y = alt.hist[[i]]$y, Type = paste0("Selection with starting MAF ", res.tbl.best$mutant.start[i], "/", res.tbl.best$nmito.start[i])))
  alt_df <- do.call(rbind, alt_df)

  df<-rbind(obs_df,null_df,alt_df)


  p_main <- ggplot() +
    geom_line(data = df, aes(x = .data$x, y = .data$y, color = .data$Type), linewidth = 1)+
    geom_line(data = df, aes(x = .data$x, y = .data$y, color = .data$Type), linewidth = 1) +
    geom_line(data = df, aes(x = .data$x, y = .data$y, color = .data$Type), linewidth = 1) +
    scale_color_manual(values = c("black", res.tbl.best$clr[best.null], res.tbl.best$clr[-best.null]),
                       labels = c("Observed", "No Selection", paste0("Selection with starting MAF ", res.tbl.best$mutant.start[-best.null], "/", res.tbl.best$nmito.start[-best.null]))) +
    labs(x = "Mutant Allele Frequency", y = "Total-Read Adjusted Histogram",color=NULL,
         title = plot.title)+

    theme(
      panel.background = element_rect(fill="white", colour="white", linewidth=0.5,
                                      linetype="solid", color="white"),
      panel.border = element_blank(),

      plot.title = element_text(color="black", size=14, face="bold.italic", hjust=0),
      axis.title.x = element_text(color="black", size=14, face="plain"),
      axis.title.y = element_text(color="black", size=14, face="plain"),

      legend.key=element_rect(fill='white'),
      legend.position = c(0.8,0.9),

      axis.line = element_line(linewidth = 0.5, linetype = "solid",colour = "black"),
      axis.text.x = element_text(face="plain", color="black",
                                 size=12, angle=0),
      axis.text.y = element_text(face="plain", color="black",
                                 size=12, angle=0))+
    xlim(0, 1)


  ## calculate OR
  if(show.logOR==TRUE){
    y.max <- max(df$y)
    x_range <- range(df$x)

    best.model=res.tbl.best[-best.null,]

    u=(0:1000)/1000
    clr.func=colorRamp(c("cyan","skyblue","yellow","gold"))


    p_logOR<-list()
    for (i in best.alt){
      s<-best.model$beta0[i]+best.model$beta1[i]*u+best.model$beta2[i]*u^2+best.model$beta3[i]*u^3
      s01=exp(s)/(1+exp(s))
      clrs=rgb(clr.func(s01),maxColorValue=255)
      df_segments <- data.frame(u = u, y = rep(0.05 * y.max, length(u)), col = clrs)
      df_text <- data.frame(v = seq(0, 1, length.out = 11), y = rep(0.025 * y.max, 11),
                            sv = best.model$beta0[i] +
                              best.model$beta1[i] * seq(0, 1, length.out = 11) +
                              best.model$beta2[i] * seq(0, 1, length.out = 11)^2 +
                              best.model$beta3[i] * seq(0, 1, length.out = 11)^3)

      p_logOR[[i]] <- ggplot() +
        geom_segment(data = df_segments, aes(x = .data$u, xend = .data$u, y = 0, yend = .data$y, color = col)) +
        geom_text(data = df_text, aes(x = .data$v, y = .data$y, label = round(.data$sv, 1))) +
        scale_color_identity() +
        theme_void() +
        theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) + # Adjust margins
        labs(title = paste0("Odds ratio from selection model with starting MAF ", res.tbl.best$mutant.start[i], "/", res.tbl.best$nmito.start[i]))+
        xlim(x_range)
    }
  }
  # Arrange plots using patchwork
  plot_height<-c(9,rep(1/length(best.alt),length(best.alt)))
  if (show.logOR == TRUE) {
    p_main + p_logOR + plot_layout(ncol = 1, heights = plot_height)
  } else {
    p_main
  }

}


#' @title Calculate the histogram data for the mutant allele fraction in the mtDNA mutation model
#' @description Calculate the histogram data for the mutant allele fraction in the mtDNA mutation model
#' @param read.data Read count matrix, one row per cell, first column with mutant counts, second with reference counts
#' @param result.tbl Table with columns: nmito.start, mutant.start, generation, beta0-beta3, clr(optional)
#' @param nbin Number of bin for histogram
#' @return A list includes two numeric vectors: first with MAF point, second with adjusted mutant read frequency
#' @examples
#' res.tbl.best = get.best.model(res.tbl = res.tbl,K = 10,show.best = TRUE)
#' best.null = which(rowSums(res.tbl.best[,grep("beta",colnames(res.tbl.best))])==0)
#' null.hist = hg.hist(read.data = reads,                  ## example for the model without selection
#'                    result.tbl = res.tbl.best[best.null,],
#'                    nbin = 1000)
#' @export
hg.hist=function(read.data,
                 result.tbl,
                 nbin=10000)
{
  x=(0:nbin)/nbin
  y=rep(0,length(x))

  mdl.maf=maf.pmf(result.tbl$mutant.start,
                  result.tbl$nmito.start-result.tbl$mutant.start,
                  result.tbl$generation,
                  beta=c(result.tbl$beta0,
                         result.tbl$beta1,
                         result.tbl$beta2,
                         result.tbl$beta3))

  n=nrow(read.data)
  for (i in 1:n)
  {
    nr=sum(read.data[i,])
    rc.pmf=read.pmf(0:nr,nr:0,mdl.maf)
    rc.indx=ceiling(x*(nr+1))
    rc.indx[1]=1
    y=y+rc.pmf[rc.indx]*(nr+1)

  }
  y=y/n

  res=list(x=x,y=y)
  return(res)

}


#' @title Calculate the histogram data for the mutant allele fraction in the observed data
#' @description Calculate the histogram data for the mutant allele fraction in the observed data
#' @param read.data A matrix with columns: mut.read and wt.read, representing observed mutant and wildtype mitochondria read count for each single cell
#' @param nbin Number of bin for histogram
#' @return A data.frame: first column with MAF point, second with adjusted mutant read frequency
#' @examples
#' obs.hist = sc.hist(read.data = reads, nbin = 1000)
#' @export
sc.hist=function(read.data,nbin=10000)
{

  x=seq(from=0,to=1,length=nbin+1)  # set of MAF points for plotting
  y=rep(0,length(x))                # height of histogram at those points

  n=nrow(read.data)
  for (i in 1:n)
  {
    r=sum(read.data[i,1:2])
    m=read.data[i,1]
    indx=ceiling(nbin*m/(r+1)+2-(m==0)):floor(nbin*(m+1)/(r+1)+1)
    indx=unique(indx)
    hght=(r+1)
    y[indx]=y[indx]+hght
  }
  y=y/n

  res=cbind.data.frame(x=x,y=y)
  return(res)
}


#' @title Plot beta-smoothed histogram of mtDNA for the mutant allele fraction in the mtDNA mutation model
#' @description Plot beta-smoothed histogram of mtDNA for the mutant allele fraction in observed data and modeling
#' @param read.data A matrix with columns: mut.read and wt.read, representing observed mutant and wildtype mitochondria read count for each single cell
#' @param res.tbl.best A data frame from get.best.model() output with columns: nmito.start, mutant.start, generation, beta0-beta3, clr(optional), each row represents a model
#' @param plot.title A character string, default is NULL
#' @return Beta-smoothed histogram
#' @examples
#' res.tbl.best = get.best.model(res.tbl = res.tbl,K = 10,show.best = TRUE)
#' res.tbl.best$clr = c("red","green","grey")
#' plt.hist.smooth(read.data = reads,
#'                  res.tbl.best = res.tbl.best,
#'                  plot.title = NULL)
#' @seealso \code{\link{hg.hist.smooth}}
#' @importFrom grDevices colorRamp rgb rainbow
#' @importFrom stats dbeta
#' @importFrom ggplot2 ggplot geom_line scale_color_manual labs theme geom_segment geom_text scale_color_identity theme_void aes xlim element_blank element_line element_rect element_text margin
#' @importFrom rlang .data
#' @export
plt.hist.smooth=function(read.data,res.tbl.best,plot.title = NULL)
{


  # check input
  req.clms=c("nmito.start","mutant.start",
             "generation",
             "beta0","beta1","beta2","beta3")

  if (any(!(req.clms%in%colnames(res.tbl.best))))
    stop("res.tbl.best missing some required columns.")

  if (nrow(res.tbl.best)<2)
    stop("res.tbl.best must have one selection-free model and one selection model")

  if (is.null(res.tbl.best$clr))
    res.tbl.best$clr=rainbow(nrow(res.tbl.best))

  mut.read=read.data[,"mut.read"]
  wt.read=read.data[,"wt.read"]
  tot.read=mut.read+wt.read

  u=0:1000/1000
  obs.hist=rep(0,length(u))
  for (i in 1:length(tot.read)){
    obs.hist=obs.hist+dbeta(u,mut.read[i]+1,wt.read[i]+1)
  }
  obs.hist=obs.hist/length(tot.read)



  best.null<-which(rowSums(res.tbl.best[,grep("beta",colnames(res.tbl.best))])==0)
  null.maf=maf.pmf(res.tbl.best$mutant.start[best.null],
                   res.tbl.best$nmito.start[best.null]-res.tbl.best$mutant.start[best.null],
                   res.tbl.best$generation[best.null])
  null.mdl=res.tbl.best[best.null,]
  null.hist=hg.hist.smooth(null.maf,read.data,null.mdl)
  #null.hist=null.hist/length(tot.read)



  best.alt<-which(rowSums(res.tbl.best[,grep("beta",colnames(res.tbl.best))])!=0)
  alt.hist<-list()
  for (i in best.alt){
    alt.mdl<-res.tbl.best[i,]
    alt.beta=unlist(res.tbl.best[i,c("beta0","beta1","beta2","beta3")])

    alt.maf<-maf.pmf(res.tbl.best$mutant.start[i],
                     res.tbl.best$nmito.start[i]-res.tbl.best$mutant.start[i],
                     res.tbl.best$generation[i],
                     alt.beta)

    alt.hist[[i]]=hg.hist.smooth(alt.maf,read.data,alt.mdl)

    #alt.hist[[i]]=alt.hist[[i]]/length(tot.read)
  }



  obs_df <- data.frame(x = u, y = obs.hist, Type = "Observed")
  null_df <- data.frame(x = u, y = null.hist, Type = "No selection")
  alt_df <- lapply(1:nrow(res.tbl.best[best.alt,]), function(i) data.frame(x = u, y = alt.hist[[i]], Type = paste0("Selection with starting MAF ", res.tbl.best$mutant.start[i], "/", res.tbl.best$nmito.start[i])))
  alt_df <- do.call(rbind, alt_df)

  df<-rbind(obs_df,null_df,alt_df)

  ggplot() +
    geom_line(data = df, aes(x = .data$x, y = .data$y, color = .data$Type), linewidth = 1) +
    geom_line(data = df, aes(x = .data$x, y = .data$y, color = .data$Type), linewidth = 1) +
    geom_line(data = df, aes(x = .data$x, y = .data$y, color = .data$Type), linewidth = 1) +
    scale_color_manual(values = c("black", res.tbl.best$clr[best.null], res.tbl.best$clr[-best.null]),
                       labels = c("Observed", "No-selection", paste0("Selection with starting MAF ", res.tbl.best$mutant.start[-best.null], "/", res.tbl.best$nmito.start[-best.null]))) +
    labs(x = "Mutant Allele Frequency", y = "Beta-Smoothed Total-Read Adjusted Histogram",color=NULL,
         title = plot.title)+
    theme(
      panel.background = element_rect(fill="white", colour="white", linewidth=0.5,
                                      linetype="solid", color="white"),
      panel.border = element_blank(),

      plot.title = element_text(color="black", size=14, face="bold.italic", hjust=0),
      axis.title.x = element_text(color="black", size=14, face="plain"),
      axis.title.y = element_text(color="black", size=14, face="plain"),

      legend.key=element_rect(fill='white'),
      legend.position = c(0.8,0.9),

      axis.line = element_line(linewidth = 0.5, linetype = "solid",colour = "black"),
      axis.text.x = element_text(face="plain", color="black",
                                 size=12, angle=0),
      axis.text.y = element_text(face="plain", color="black",
                                 size=12, angle=0)
    )

}


#' @title Calculate the beta-smoothed histogram data for the mutant allele fraction in the mtDNA mutation model
#' @description Calculate the beta-smoothed histogram data for the mutant allele fraction in the mtDNA mutation model
#' @param maf Matrix of the distribution of mutant allele fraction for a population of cells
#' @param read.data A matrix with columns: mut.read and wt.read, representing observed mutant and wildtype mitochondria read count for each single cell
#' @param mdl Table with columns: nmito.start, mutant.start, generation, beta0-beta3, clr(optional)
#' @return A vector of beta-smoothed mutant read frequency
#' @examples
#' res.tbl.best = get.best.model(res.tbl = res.tbl,K = 10,show.best = TRUE)
#' best.null = which(rowSums(res.tbl.best[,grep("beta",colnames(res.tbl.best))])==0)
#' null.mdl = res.tbl[best.null,]
#' null.maf = maf.pmf(res.tbl.best$mutant.start[best.null],
#'                  res.tbl.best$nmito.start[best.null]-res.tbl.best$mutant.start[best.null],
#'                  res.tbl.best$generation[best.null])
#' smooth.hist = hg.hist.smooth(maf = null.maf,read.data = reads,mdl = null.mdl)
#' @seealso \code{\link{maf.pmf}}, \code{\link{read.pmf}}
#' @importFrom stats dbeta
#' @export
hg.hist.smooth=function(maf,read.data,mdl){


  tot.read=read.data[,"mut.read"]+read.data[,"wt.read"]


  u=0:1000/1000
  hist=rep(0,length(u))


  for (i in 1:length(tot.read))
  {

    pmf=read.pmf(0:tot.read[i],
                 tot.read[i]:0,
                 maf,
                 mdl$generation)

    for (j in 0:tot.read[i])
    {
      hist=hist+pmf[j+1]*dbeta(u,j+1,tot.read[i]-j+1)
    }

  }
  hist=hist/length(tot.read)
  return(hist)
}
