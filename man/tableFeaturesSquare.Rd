% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tableFeaturesSquare.R
\name{tableFeaturesSquare}
\alias{tableFeaturesSquare}
\title{Title
Show selected features as Table (for Squared one)}
\usage{
tableFeaturesSquare(Z, index, id)
}
\arguments{
\item{Z}{Tensor of features}

\item{index}{List that includes selected features and P-values}

\item{id}{feature to be shown}
}
\value{
Table list of selected features
}
\description{
Title
Show selected features as Table (for Squared one)
}
\examples{
omics1 <- matrix(runif(100000),ncol=10)
dimnames(omics1) <- list(seq_len(10000),seq_len(10))
omics2 <- matrix(runif(100000),ncol=10)
dimnames(omics2) <- dimnames(omics1)
Multi <- list(omics1,omics2)
Z <- PrepareSummarizedExperimentTensorSquare(
sample=matrix(colnames(omics1),1),
feature=list(omics1=rownames(omics1),
omics2=rownames(omics2)),
value=convertSquare(Multi),
sampleData=list(NA))
HOSVD <- computeHosvdSqure(Z)
cond <- list(0,rep(seq_len(2),each=5),c("A","B"))
input_all <- selectSingularValueVectorLarge(HOSVD,cond,input_all=c(1,1))
index <- selectFeatureSquare(HOSVD,input_all,Multi,de=c(0.1,0.1),
interact=FALSE)
index[[1]]$index[1:100]<-TRUE
index[[1]]$p.value[1:100] <- 1e-3
tableFeaturesSquare(Z,index,1)
}
