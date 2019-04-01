### functions for Spectrum

if(getRversion() >= "2.15.1")  utils::globalVariables(c('K','PC1','PC2','X1','X2','Z','evals','x'))

#' CNN_kernel_mine_b: fast adaptive density aware kernel
#'
#' @param mat Matrix: matrix should have samples as columns and rows as features
#' @param NN Numerical value: the number of nearest neighbours to use when calculating local sigma
#' @param NN2 Numerical value: the number of nearest neighbours to use when calculating common nearest neighbours
#'
#' @return A kernel matrix
#' @export
#'
#' @examples
#' CNN_kern <- CNN_kernel_mine_b(blobs[,1:50])
CNN_kernel_mine_b <- function(mat, NN = 3, NN2 = 7) {
  n <- ncol(mat)
  ## need N nearest neighbour distance per sample (kn), 
  ## and names of NN2 nearest neigbours (nbs)
  nbs <- list()
  dm <- Rfast::Dist(t(mat))
  dimnames(dm) <- list(colnames(mat), colnames(mat))
  kn <- c()
  ## find kth nearest neighbour for each sample 1...N
  for (i in seq_len(n)) {
    # sort the vector to retrieve the N nearest neighbour and the names of the NN2 nearest neighbours
    sortedvec <- sort.int(dm[i, ])
    # append the NNth nearest neighbour distance
    kn <- c(kn, sortedvec[NN + 1])
    # append the names of the NN2 nearest neighbours
    nbs[[i]] <- names(sortedvec[2:(NN2+1)])
    names(nbs)[[i]] <- names(sortedvec)[1]
  }
  ## make the symmetrical matrix of kth nearest neighbours distances
  sigmamatrix <- kn %o% kn
  ## calculate the kernel using the local statistics (scale and density) of each sample pair
  out <- matrix(nrow = n, ncol = n)  # don't overwrite functions
  # calculate the numerator beforehand
  upper <- -dm^2
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      # shared nearest neighbours between ith and jth sample
      cnns <- length(intersect(nbs[[i]], nbs[[j]]))
      upperval <- upper[i, j]
      # retrieve sigma
      localsigma <- sigmamatrix[i, j]
      # calculate local affinity
      out[i, j] <- exp(upperval / (localsigma * (cnns + 1)))
    }
  }
  # reflect, diag = 1
  out <- pmax(out, t(out), na.rm = TRUE)
  diag(out) <- 1
  ## return kernel
  return(out)
}

#' rbfkernel_b: fast self-tuning kernel
#'
#' @param mat Matrix: matrix should have samples as columns and rows as features
#' @param K Numerical value: the number of nearest neighbours to use when calculating local sigma
#' @param sigma Numerical value: a global sigma, usually left to 1 which has no effect
#'
#' @return A kernel matrix
#' @export
#'
#' @examples
#' stsc_kern <- rbfkernel_b(blobs[,1:50])
rbfkernel_b <- function (mat, K = 3, sigma = 1) { # calculate gaussian kernel with local sigma
  n <- ncol(mat)
  NN <- K # nearest neighbours (2-3)
  dm <- Rfast::Dist(t(mat))
  kn <- c() # find kth nearest neighbour for each sample 1...N
  for (i in seq_len(n)) {
    sortedvec <- as.numeric(sort.int(dm[i, ]))
    sortedvec <- sortedvec[!sortedvec == 0]
    kn <- c(kn, sortedvec[NN])
  }
  sigmamatrix <- kn %o% kn # make the symmetrical matrix of kth nearest neighbours distances
  upper <- -dm^2 # calculate the numerator beforehand
  out <- matrix(nrow = n, ncol = n)
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      lowerval <- sigmamatrix[i, j] # retrieve sigma
      upperval <- upper[i, j]
      out[i, j] <- exp(upperval / (lowerval * sigma)) # calculate local affinity
    }
  }
  # reflect, diag = 1
  out <- pmax(out, t(out), na.rm = TRUE)
  diag(out) <- 1
  ## return kernel
  return(out)
}

### plot eigengap
plot_egap <- function(d,fontsize=22,width=26,maxk=maxk,dotsize=dotsize){
  py <- ggplot2::ggplot(data=d,aes(x=K,y=evals)) + ggplot2::geom_point(colour = 'grey', size = dotsize) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = fontsize, colour = 'black'),
          axis.text.x = ggplot2::element_text(size = fontsize, colour = 'black'),
          axis.title.x = ggplot2::element_text(size = fontsize),
          axis.title.y = ggplot2::element_text(size = fontsize),
          legend.text = ggplot2::element_text(size = fontsize),
          legend.title = ggplot2::element_text(size = fontsize),
          plot.title = ggplot2::element_text(size = fontsize, colour = 'black', hjust = 0.5),
          panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::ylab('Eigenvalue') +
    ggplot2::xlab('Eigenvector') +
    ggplot2::scale_x_continuous(limits=c(1,maxk+1),breaks=c(seq(1,maxk+1,by=1)))
  print(py)
}

#' tsne: A tsne function for similarity matrices or ordinary data
#'
#' @param mydata Data frame or matrix: kernel matrix or data frame with samples as columns, features as rows
#' @param labels Factor: to label the plot with colours
#' @param seed Numerical value: to repeat the results exactly, setting seed is required
#' @param perplex Numerical value: this is the perplexity parameter for tsne, it usually requires adjusting for each dataset
#' @param axistextsize Numerical value: axis text size
#' @param legendtextsize Numerical value: legend text size
#' @param dotsize Numerical value: dot size
#' @param similarity Logical flag: whether input is similarity matrix or not
#' 
#' @return A tsne plot object
#' @export
#'
#' @examples
#' ex_tsne <- tsne(blobs[,1:50],perplex=15,similarity=FALSE)
tsne <- function(mydata, labels=FALSE, perplex=15, seed=FALSE, axistextsize = 18,
                 legendtextsize = 18, dotsize = 3, similarity = TRUE){
  if (similarity == TRUE){
    # fix up similarity matrix and convert to distance
    mydata[lower.tri(mydata)] <- t(mydata)[lower.tri(mydata)]
    diag(mydata) <- 1
    mydata <- 1-mydata
  }
  if (seed != FALSE){
    set.seed(seed)
  }
  if (similarity == TRUE){
    tsne <- Rtsne::Rtsne(t(as.matrix(mydata)), dims = 2, perplexity=perplex, verbose=FALSE, max_iter = 500, is_distance = TRUE)
  }else{
    tsne <- Rtsne::Rtsne(t(as.matrix(mydata)), dims = 2, perplexity=perplex, verbose=FALSE, max_iter = 500)
  }
  scores <- data.frame(tsne$Y) # PC score matrix
  p <- ggplot2::ggplot(data = scores, aes(x = X1, y = X2) ) + ggplot2::geom_point(aes(colour = factor(labels)), size = dotsize) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_text(size = axistextsize, colour = 'black'),
          axis.text.x = ggplot2::element_text(size = axistextsize, colour = 'black'),
          axis.title.x = ggplot2::element_text(size = axistextsize),
          axis.title.y = ggplot2::element_text(size = axistextsize),
          legend.title = ggplot2::element_text(size = legendtextsize),
          legend.text = ggplot2::element_text(size = legendtextsize)) 
  print(p)
  return(p)
}

#' umap: A umap function for similarity matrices or ordinary data
#'
#' @param mydata Data frame or matrix: kernel matrix or data frame with samples as columns, features as rows
#' @param labels Factor: to label the plot with colours
#' @param axistextsize Numerical value: axis text size
#' @param legendtextsize Numerical value: legend text size
#' @param dotsize Numerical value: dot size
#' @param similarity Logical flag: whether input is similarity matrix or not
#' 
#' @return A umap plot object
#' @export
#'
#' @examples
#' ex_umap <- umap(blobs[,1:50],similarity=FALSE)
umap <- function(mydata, labels=FALSE, dotsize = 3, similarity=TRUE, axistextsize = 18, legendtextsize = 18){
  # fix up similarity matrix and convert to distance
  if (similarity == TRUE){
    mydata[lower.tri(mydata)] <- t(mydata)[lower.tri(mydata)]
    diag(mydata) <- 1
    mydata <- 1-mydata
  }
  # do umap
  if (similarity == TRUE){
    umap <- umap::umap(t(as.matrix(mydata)),input="dist")
  }else{
    umap <- umap::umap(t(as.matrix(mydata)))
  }
  # do plot
  scores <- data.frame(umap$layout) # PC score matrix
  p <- ggplot2::ggplot(data = scores, aes(x = X1, y = X2) ) + ggplot2::geom_point(aes(colour = factor(labels)), size = dotsize) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size = axistextsize, colour = 'black'),
                   axis.text.x = ggplot2::element_text(size = axistextsize, colour = 'black'),
                   axis.title.x = ggplot2::element_text(size = axistextsize),
                   axis.title.y = ggplot2::element_text(size = axistextsize),
                   legend.title = ggplot2::element_text(size = legendtextsize),
                   legend.text = ggplot2::element_text(size = legendtextsize)) 
  print(p)
  return(p)
}

#' pca: A pca function
#'
#' @param mydata Data frame or matrix: matrix or data frame with samples as columns, features as rows
#' @param labels Factor: to label the plot with colours
#' @param axistextsize Numerical value: axis text size
#' @param legendtextsize Numerical value: legend text size
#' @param dotsize Numerical value: dot size
#' 
#' @return A pca plot object
#' @export
#'
#' @examples
#' ex_pca <- pca(blobs[,1:50])
pca <- function(mydata, labels=FALSE, dotsize = 3, axistextsize = 18, legendtextsize = 18){
  pca1 = prcomp(t(mydata))
  scores <- data.frame(pca1$x) # PC score matrix
  p <- ggplot2::ggplot(data = scores, aes(x = PC1, y = PC2) ) + ggplot2::geom_point(aes(colour = factor(labels)), size = dotsize) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_text(size = axistextsize, colour = 'black'),
          axis.text.x = ggplot2::element_text(size = axistextsize, colour = 'black'),
          axis.title.x = ggplot2::element_text(size = axistextsize),
          axis.title.y = ggplot2::element_text(size = axistextsize),
          legend.title = ggplot2::element_text(size = legendtextsize),
          legend.text = ggplot2::element_text(size = legendtextsize)) 
  print(p)
}

### display similarity matrix in heatmap
displayClusters <- function(W, group, fsize = 1.5) {
  ind <- sort(as.vector(group),index.return=TRUE)
  ind <- ind$ix
  diag(W) <- 0
  W <- W / rowSums(W)
  W <- W + t(W)
  cols <- c('#FFFFFF',colorRampPalette(RColorBrewer::brewer.pal(9,'Blues'))(100))
  image(1:ncol(W),1:nrow(W),W[ind,ind],col=cols,xlab = 'Samples',ylab='Samples',
        cex.axis=fsize,cex.lab=fsize)
  title(main = paste('K =',max(group)), cex.main = fsize, font.main = 1)
}

### search multimodality diffs to find best k
findk <- function(res,maxk=maxk,thresh=4,frac=2){
  v <- diff(res[,2])
  v <- v[1:maxk] 
  # parameters for search used in paper
  # thresh <- 4
  # frac <- 2
  for (e in seq(2,length(v))){
    if (e == 2){ # store element
      saved <- v[e]
      index <- 2
    }else{
      currentindex <- e
      xx <- currentindex-index # how far are we ahead of minima
      if (xx >= thresh){
        # if we have a max followed by 4 negs, then stop search
        break
      }
      if (v[e] < saved | v[e] < saved/frac){
        # we have found a new max diff, so replace cache and index
        saved <- v[e]
        index <- e
      }
    }
  }
  optk <- index
  return(optk)
}

### tuner for locally adaptive density aware kernel
kernfinder_mine <- function(data,maxk=10,fontsize=fontsize,silent=silent,
                            showres=showres,dotsize=dotsize){ 
  if (silent == FALSE){
    message('finding optimal NN kernel parameter by examining eigenvector distributions')
  }
  rr <- c()
  for (param in seq(1,10)){
    if (silent == FALSE){
      message(paste('tuning kernel NN parameter:',param))
    }
    kern <- CNN_kernel_mine_b(data,NN=param,NN2=7)
    kern[which(!is.finite(kern))] <- 0 # deal with possible NaNs
    ## calculate difference between most multimodal eigenvectors and background
    dv <- 1/sqrt(rowSums(kern)) # D = diag(1/sqrt(rowSums(A)))
    l <- dv * kern %*% diag(dv) # L = D%*%A%*%D
    xi <- eigen(l)$vectors
    res <- matrix(nrow=ncol(xi),ncol=2)
    for (ii in seq(1,ncol(xi))){
      r <- diptest::dip.test(xi[,ii], simulate.p.value = FALSE, B = 2000)
      res[ii,1] <- r$p.value
      res[ii,2] <- r$statistic
    }
    ## chose kernel via maximum difference
    diffs <- diff(res[,2])
    diffs <- diffs[-1] # remove difference corresponding to first eigenvector
    tophit <- diffs[1:(maxk+1)][which.min((diffs[1:(maxk+1)]))]
    rr <- c(rr, tophit)
  }
  ## find index that yields lowest diff
  optimalparam <- which.min(rr)
  if (silent == FALSE){
    message(paste('optimal NN:'),optimalparam)
  }
  ## print tuning res
  d <- data.frame(x=seq(1,10),y=rr)
  py <- ggplot2::ggplot(data=d,aes(x=x,y=rr)) + ggplot2::geom_point(colour = 'black', size = dotsize) +
    ggplot2::theme_bw() +
    ggplot2::geom_line() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = fontsize, colour = 'black'),
          axis.text.x = ggplot2::element_text(size = fontsize, colour = 'black'),
          axis.title.x = ggplot2::element_text(size = fontsize),
          axis.title.y = ggplot2::element_text(size = fontsize),
          legend.text = ggplot2::element_text(size = fontsize),
          legend.title = ggplot2::element_text(size = fontsize),
          plot.title = ggplot2::element_text(size = fontsize, colour = 'black', hjust = 0.5),
          panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::ylab('D') +
    ggplot2::xlab('NN') +
    ggplot2::scale_x_continuous(limits=c(1,10),breaks=c(seq(1,10,by=1)))
  if (showres == TRUE){
    print(py)
  }
  ##
  return(optimalparam)
}

### tuner for self tuning kernel
kernfinder_local <- function(data,maxk=maxk,silent=silent,dotsize=dotsize,
                             fontsize=fontsize,showres=showres){ 
  if (silent == FALSE){
    message('finding optimal kernel NN parameter by examining eigenvectors')
  }
  rr <- c()
  for (param in seq(1,10)){
    if (silent == FALSE){
      message(paste('tuning NN parameter:',param))
    }
    kern <- rbfkernel_b(data,K=param,sigma=1) # sigma (0.75)
    ## calculate difference between most multimodal eigenvectors and background
    dv <- 1/sqrt(rowSums(kern)) # D = diag(1/sqrt(rowSums(A)))
    l <- dv * kern %*% diag(dv) # L = D%*%A%*%D
    xi <- eigen(l)$vectors
    res <- matrix(nrow=ncol(xi),ncol=2)
    for (ii in seq(1,ncol(xi))){
      r <- diptest::dip.test(xi[,ii], simulate.p.value = FALSE, B = 2000)
      res[ii,1] <- r$p.value
      res[ii,2] <- r$statistic
    }
    ## chose kernel via maximum difference
    diffs <- diff(res[,2])
    diffs <- diffs[-1] # remove difference corresponding to first eigenvector
    tophit <- diffs[1:(maxk+1)][which.min((diffs[1:(maxk+1)]))]
    rr <- c(rr, tophit)
  }
  ## find index that yields lowest diff
  optimalparam <- which.min(rr)
  if (silent == FALSE){
    message(paste('optimal NN:'),optimalparam)
  }
  ## print tuning res
  d <- data.frame(x=seq(1,10),y=rr)
  py <- ggplot2::ggplot(data=d,aes(x=x,y=rr)) + ggplot2::geom_point(colour = 'black', size = dotsize) +
    ggplot2::theme_bw() +
    ggplot2::geom_line() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = fontsize, colour = 'black'),
                   axis.text.x = ggplot2::element_text(size = fontsize, colour = 'black'),
                   axis.title.x = ggplot2::element_text(size = fontsize),
                   axis.title.y = ggplot2::element_text(size = fontsize),
                   legend.text = ggplot2::element_text(size = fontsize),
                   legend.title = ggplot2::element_text(size = fontsize),
                   plot.title = ggplot2::element_text(size = fontsize, colour = 'black', hjust = 0.5),
                   panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::ylab('D') +
    ggplot2::xlab('NN') +
    ggplot2::scale_x_continuous(limits=c(1,10),breaks=c(seq(1,10,by=1)))
  if (showres == TRUE){
    print(py)
  }
  ## this is the maximum method, better to look at the distribution
  optimalparam <- which.min(rr)
  if (silent == FALSE){
    message(paste('optimal NN:'),optimalparam)
  }
  return(optimalparam)
}

### evaluate multimodality of eigenvectors to get max gap
EM_finder <- function(xi,silent=silent){ # accepts the eigenvector decomposition of L
  if (silent == FALSE){
    message('finding informative eigenvectors...')
  }
  res <- matrix(nrow=ncol(xi),ncol=2)
  for (ii in seq(1,ncol(xi))){
    r <- diptest::dip.test(xi[,ii], simulate.p.value = FALSE, B = 2000)
    res[ii,1] <- r$p.value
    res[ii,2] <- r$statistic
  }
  if (silent == FALSE){
    message('done.')
  }
  return(res)
}

### plot eigengap
plot_multigap <- function(d,fontsize=fontsize,width=26,maxk=maxk,
                          dotsize=dotsize){
  py <- ggplot2::ggplot(data=d,aes(x=K,y=Z)) + ggplot2::geom_point(colour = 'grey', size=dotsize) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = fontsize, colour = 'black'),
          axis.text.x = ggplot2::element_text(size = fontsize, colour = 'black'),
          axis.title.x = ggplot2::element_text(size = fontsize),
          axis.title.y = ggplot2::element_text(size = fontsize),
          legend.text = ggplot2::element_text(size = fontsize),
          legend.title = ggplot2::element_text(size = fontsize),
          plot.title = ggplot2::element_text(size = fontsize, colour = 'black', hjust = 0.5),
          panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::ylab('Z') +
    ggplot2::xlab('Eigenvector') +
    ggplot2::scale_x_continuous(limits=c(1,maxk+1),breaks=c(seq(1,maxk+1,by=1)))
  print(py)
}

#' kernel_pca: A kernel pca function
#'
#' @param datam Dataframe or matrix: a data frame with samples as columns, rows as features, or a kernel matrix
#' @param labels Factor: to label the plot with colours
#' @param axistextsize Numerical value: axis text size
#' @param legendtextsize Numerical value: legend text size
#' @param dotsize Numerical value: dot size
#' @param kernel Logical flag: whether the input is a kernel or not
#' 
#' @return A kernel PCA plot
#' @export
#'
#' @examples
#' ex_kernel_pca <- kernel_pca(blobs[,1:50], kernel=FALSE)
kernel_pca <- function(datam, labels = FALSE, axistextsize = 18, legendtextsize = 18, dotsize = 3,
                       kernel = TRUE){
  if (kernel == FALSE){
    km <- CNN_kernel_mine_b(datam)
  }else{
    km <- datam
  }
  # input a NXN kernel
  m <- nrow(km)
  # center kernel matrix
  kc <- t(t(km - colSums(km)/m) - rowSums(km)/m) + sum(km)/m^2
  # compute eigenvectors
  res <- eigen(kc/m,symmetric=TRUE)
  # multiply each eigenvector by the square root of its eigenvalue (get PCs)
  features <- m
  ret <- suppressWarnings(t(t(res$vectors[,1:features])/sqrt(res$values[1:features])))
  # do plot
  scores <- data.frame(ret)
  colnames(scores)[1:2] <- c('PC1','PC2')
  p1 <- ggplot(data = scores, aes(x = PC1, y = PC2) ) + geom_point(aes(colour = factor(labels)),size=3) + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = axistextsize, colour = 'black'),
          axis.text.x = element_text(size = axistextsize, colour = 'black'),
          axis.title.x = element_text(size = axistextsize),
          axis.title.y = element_text(size = axistextsize),
          legend.title = element_text(size = legendtextsize),
          legend.text = element_text(size = legendtextsize)) 
  print(p1)
  #
  return(p1)
}
