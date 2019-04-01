#' Spectrum: Fast Adaptive Spectral Clustering for Single and Multi-view Data
#'
#' Spectrum is a self-tuning spectral clustering method for single or multi-view data. Spectrum uses a new type of adaptive
#' density aware kernel that strengthens local connections in the graph. For integrating multi-view data and reducing noise
#' a tensor product graph data integration and diffusion procedure is used. Spectrum analyses eigenvector variance or distribution
#' to determine the number of clusters. Spectrum is well suited for a wide range of data, including both Gaussian and non-Gaussian
#' structures.
#'
#' @param data Data frame or list of data frames: contains the data with samples as columns and rows as features. For multi-view data a list of dataframes is to be supplied with the samples in the same order.
#' @param method Numerical value: 1 = default eigengap method (Gaussian clusters), 2 = multimodality gap method (Gaussian/ non-Gaussian clusters), 3 = no automatic method (see fixk param)
#' @param maxk Numerical value: the maximum number of expected clusters (default  = 10). This is data dependent, do not set excessively high.
#' @param fixk Numerical value: if we are just performing spectral clustering without automatic selection of K, set this parameter and method to 3
#' @param silent Logical flag: whether to turn off messages
#' @param showres Logical flag: whether to show the results on the screen
#' @param diffusion Logical flag: whether to perform graph diffusion to reduce noise, recommended for method 1 only
#' @param kerneltype Character string: 'density' (default) = adaptive density aware kernel, 'stsc' = Zelnik-Manor self-tuning kernel
#' @param NN Numerical value: kernel param, the number of nearest neighbours to use sigma parameters (default = 3)
#' @param NN2 Numerical value: kernel param, the number of nearest neighbours to use for the common nearest neigbours (default = 7)
#' @param showpca Logical flag: whether to show pca when running on one view
#' @param showheatmap Logical flag: whether to show heatmap of similarity matrix when running on one view
#' @param showdimred Logical flag: whether to show UMAP or t-SNE of final similarity matrix
#' @param visualisation Character string: what kind of dimensionality reduction to run on the similarity matrix (umap or tsne)
#' @param frac Numerical value: optk search param, fraction to find the last substantial drop (multimodality gap method param)
#' @param thresh Numerical value: optk search param, how many points ahead to keep searching (multimodality gap method param)
#' @param fontsize Numerical value: controls font size of the ggplot2 plots
#' @param dotsize Numerical value: controls the dot size of the ggplot2 plots
#' @param tunekernel Logical flag: whether to tune the kernel, only applies for method 2
#' @param clusteralg Character string: clustering algorithm for eigenvector matrix (GMM or km)
#' @param FASP Logical flag: whether to use Fast Approximate Spectral Clustering (for v. high sample numbers)
#' @param FASPk Numerical value: the number of centroids to compute when doing FASP
#'
#' @return A list, containing: 
#' 1) cluster assignments, in the same order as input data columns 
#' 2) eigenvector analysis results (either eigenvalues or dip test statistics)
#' 3) optimal K
#' 4) final similarity matrix
#' 5) eigenvectors and eigenvalues of graph Laplacian
#' @export
#'
#' @examples
#' res <- Spectrum(brain[[1]][,1:50])

Spectrum <- function(data,method=1,silent=FALSE,showres=TRUE,diffusion=TRUE,
                     kerneltype=c('density','stsc'),maxk=10,NN=3,NN2=7,
                     showpca=FALSE,showheatmap=FALSE,showdimred=FALSE,
                     visualisation=c('umap','tsne'),frac=2,thresh=7,
                     fontsize=18,dotsize=3,tunekernel=TRUE,clusteralg='GMM',
                     FASP=FALSE,FASPk=NULL,fixk=NULL){
  
  ###
  kerneltype <- match.arg(kerneltype)
  visualisation <- match.arg(visualisation)
  
  ### error handling
  if (class(data) != 'list'){
    datalist <- list(data) # just a single data source
  }else{
    datalist <- data
  }
  if ( length(datalist) > 1 & FASP == TRUE ) {
    stop("Error: FASP method works for only a single view")
  }
  if (is.null(FASPk) == TRUE & FASP == TRUE){
    stop("Error: FASP method requires a number of centroids to compute")
  }
  
  ### initial messages
  if (silent == FALSE){
    message('***Spectrum***')
    message(paste('detected views:',length(datalist)))
    message(paste('method:',method))
    message(paste('kernel:',kerneltype))
    if (FASP){
      message(paste('FASP method: TRUE'))
    }
  }
  
  ### if running FASP calculate the centroids
  if (FASP){
    cs <- kmeans(t(datalist[[1]]),centers=FASPk)
    csx <- cs$centers # these are the centroids to cluster
    cas <- cs$cluster # now samples are assigned to centroids
    datalist[[1]] <- data.frame(t(csx))
  }
  
  ### calculate individual kernels
  kernellist <- list()
  for (platform in seq(1,length(datalist))){
    ### calculate kernel
    if (silent == FALSE){
      message(paste('calculating kernel',platform))
    }
    if (kerneltype == 'stsc'){
      if (method == 2){
        if (tunekernel){
          NN <- kernfinder_local(datalist[[platform]],maxk=maxk,silent=silent,fontsize=fontsize,
                                 dotsize=dotsize,showres=showres)
        }else{
          NN <- 3
        }
      }
      kerneli <- rbfkernel_b(datalist[[platform]],K=NN,sigma=1)
    }else if (kerneltype == 'density'){ 
      if (method == 2){
        if (tunekernel){
          NN <- kernfinder_mine(datalist[[platform]],maxk=maxk,silent=silent,
                                showres=showres,fontsize=fontsize,dotsize=dotsize)
        }else{
          NN <- 3
        }
      }
      kerneli <- CNN_kernel_mine_b(datalist[[platform]],NN=NN,NN2=7)
    }
    if (silent == FALSE){
      message('done.')
    }
    ### save kernel in list
    colnames(kerneli) <- colnames(datalist[[1]])
    row.names(kerneli) <- colnames(datalist[[1]])
    kernellist[[platform]] <- kerneli
  }
  
  ### fuse and truncate/ make KNN graph (both)
  if (silent == FALSE){
    message('combining kernels if > 1 and making KNN graph...')
  }
  A <- Reduce('+', kernellist) # construct A by linear combination
  
  ### diffusion on KNN graph
  if (diffusion == TRUE){
    ## get KNN graph
    for (col in seq(1,ncol(A))){
      KNNs <- head(rev(sort(A[,col])),(10+1)) # find the KNNs (10 default)
      tokeep <- names(KNNs)
      A[!(names(A[,col])%in%tokeep),col] <- 0
    }
    A <- A/rowSums(A) # row normalise A
    if (silent == FALSE){
      message('done.')
    }
    ## diffusion iterations
    if (silent == FALSE){
      message('diffusing on tensor graph...')
    }
    Qt <- A
    im <- matrix(ncol=ncol(A),nrow=ncol(A))
    im[is.na(im)] <- 0
    diag(im) <- 1
    for (t in seq(2,5)){ # (5)
      Qt <- A%*%Qt%*%t(A)+im
    }
    A2 <- t(Qt)
    if (silent == FALSE){
      message('done.')
    }
    #diag(A2) <- 0 # removing self similarity makes no difference
  }else if (diffusion == FALSE){
    # if we are not doing TPG method use simple mean A
    A2 <- A/length(datalist)
  }
  
  ### calculate graph Laplacian of A2 (both)
  if (silent == FALSE){
    message('calculating graph laplacian...')
  }
  dv <- 1/sqrt(rowSums(A2))
  l <- dv * A2 %*% diag(dv)
  
  ### eigengap heuristic
  if (method == 1){
    if (silent == FALSE){
      message('getting eigendecomposition of L...')
    }
    decomp <- eigen(l)
    if (silent == FALSE){
      message('done.')
      message('examining eigenvalues to select K...')
    }
    evals <- as.numeric(decomp$values)
    diffs <- diff(evals)
    diffs <- diffs[-1]
    optk <- which.max(abs(diffs[1:maxk-1]))+1 # put a cap on the maximum number of clusters
    if (silent == FALSE){
      message(paste('optimal K:',optk))
    }
    nn <- maxk+1
    d <- data.frame(K=seq(1,maxk+1),evals=evals[1:nn])
    if (showres == TRUE){
      plot_egap(d,maxk=maxk,dotsize=dotsize,
                fontsize=fontsize)
    }
  }else if (method == 2){ # multimodality gap heuristic
    if (silent == FALSE){
      message('getting eigendecomposition of L...')
    }
    decomp <- eigen(l)
    if (silent == FALSE){
      message('done.')
      message('examining eigenvector distributions to select K...')
    }
    xi <- decomp$vectors[,1:(maxk+1)]
    res <- EM_finder(xi,silent=silent)
    d <- data.frame('K'=seq(1,maxk+1),'Z'=res[1:(maxk+1),2])
    if (showres == TRUE){
      plot_multigap(d,maxk=maxk,dotsize=dotsize,
                    fontsize=fontsize)
    }
    optk <- findk(res,maxk=maxk,frac=frac,thresh=thresh)
    if (silent == FALSE){
      message(paste('optimal K:',optk))
    }
  }else if (method == 3){ # no automatic K selection method
    decomp <- eigen(l)
    optk <- fixk
  } 
  
  ### select optimal eigenvectors
  xi <- decomp$vectors[,1:optk]
  # normalise rows
  yi <- xi/sqrt(rowSums(xi^2))
  # replace NA values (row = 0) with zeros
  yi[which(!is.finite(yi))] <- 0
  
  if (clusteralg == 'GMM'){
    ### GMM
    if (silent == FALSE){
      message('doing GMM clustering...')
    }
    gmm <- ClusterR::GMM(yi, optk, verbose = F, seed_mode = "random_spread") # use random spread          
    pr <- ClusterR::predict_GMM(yi, gmm$centroids, gmm$covariance_matrices, gmm$weights)
    names(pr)[3] <- 'cluster'
    if (0 %in% pr$cluster){
      pr$cluster <- pr$cluster+1
    }
    if (silent == FALSE){
      message('done.')
    }
  }else if (clusteralg == 'km'){
    ### k means
    if (silent == FALSE){
      message('doing k means clustering...')
    }
    pr <- kmeans(yi, optk)
    if (silent == FALSE){
      message('done.')
    }
  }
  
  ### display clusters using heatmap and tsne
  if (showres == TRUE){ # for any number of data sources
    if (showheatmap == TRUE){
      displayClusters(A2,group=pr$cluster,fsize=1)
    }
    if (showdimred == TRUE && visualisation == 'umap'){ # method 1
      message('running UMAP on similarity...')
      umap(A2,labels=as.factor(pr$cluster),axistextsize=fontsize,legendtextsize=fontsize,dotsize=dotsize)
      message('done.')
    }
    if (showdimred == TRUE && visualisation == 'tsne'){ # method 2
      message('running t-SNE on similarity...')
      tsne(A2,labels=as.factor(pr$cluster),axistextsize=fontsize,legendtextsize=fontsize,dotsize=dotsize)
      message('done.')
    }
  }
  if (length(datalist) == 1 && showres == TRUE){ # for one data source only
    if (showpca == TRUE){
      pca(datalist[[1]],labels=as.factor(pr$cluster),axistextsize=fontsize,legendtextsize=fontsize,dotsize=dotsize)
    }
  }
  
  if (method != 3){
    if (FASP){
      ### if running FASP assign original samples to centroid clusters
      casn <- cas
      casn <- casn[seq_along(casn)] <- pr$cluster[as.numeric(casn[seq_along(casn)])]
      names(casn) <- names(cas) 
      ### return results
      results <- list('allsample_assignments'=casn,'centroid_assignments'=pr$cluster,
                      'eigenvector_analysis'=d,'K'=optk,'similarity_matrix'=A2,
                      'eigendecomposition'=decomp)
    }else{
      ### return results
      results <- list('assignments'=pr$cluster,'eigenvector_analysis'=d,
                      'K'=optk,'similarity_matrix'=A2,'eigendecomposition'=decomp)
    }
  }else if (method == 3){
    if (FASP){
      ### if running FASP assign original samples to centroid clusters
      casn <- cas
      casn <- casn[seq_along(casn)] <- pr$cluster[as.numeric(casn[seq_along(casn)])]
      names(casn) <- names(cas) 
      ### return results
      results <- list('allsample_assignments'=casn,'centroid_assignments'=pr$cluster,
                      'K'=optk,'similarity_matrix'=A2,
                      'eigendecomposition'=decomp)
    }else{
      ### return results
      results <- list('assignments'=pr$cluster,
                      'K'=optk,'similarity_matrix'=A2,'eigendecomposition'=decomp)
    }
  }
  
  if (silent == FALSE){
    message('finished.')
  }
  
  return(results)
}
  

