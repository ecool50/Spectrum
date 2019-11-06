## plotting functions for spectrum

if(getRversion() >= "2.15.1")  utils::globalVariables(c('K','PC1','PC2','X1','X2','Z','evals','x'))

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
#' @param similarity Logical flag: whether the input is a similarity matrix or not
#' 
#' @return A kernel PCA plot
#' @export
#'
#' @examples
#' ex_kernel_pca <- kernel_pca(blobs[,1:50], similarity=FALSE)
kernel_pca <- function(datam, labels = FALSE, axistextsize = 18, legendtextsize = 18, dotsize = 3,
                       similarity = TRUE){
  if (similarity == FALSE){
    km <- CNN_kernel(datam)
  }else{
    km <- datam
  }
  # input a NXN kernel
  m <- nrow(km)
  # center kernel matrix
  kc <- t(t(km - colSums(km)/m) - rowSums(km)/m) + sum(km)/m^2
  # compute eigenvectors
  res <- eigen(kc/m,symmetric=TRUE)
  # multiply each eigenvector by the square root of its eigenvalue
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
