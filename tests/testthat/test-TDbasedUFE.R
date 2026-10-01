test_that("multiplication works", {
require(GenomicRanges)
require(rTensor)
f <- system.file("extdata", "quickstart_counts.rds",
    package = "TDbasedUFE", mustWork = TRUE)
example <- readRDS(f)
counts <- example$counts
samples <- example$samples
stopifnot(identical(dim(counts), c(10000L, 6L)))
stopifnot(identical(colnames(counts), samples$run))

Z <- PrepareSummarizedExperimentTensor(matrix(samples$sample,c(3,2)),
    rownames(counts),array(counts,c(dim(counts)[1],3,2)))
HOSVD <- computeHosvd(Z)
input_all <- selectSingularValueVectorSmall(HOSVD,input_all=c(1,2))
index <- selectFeature(HOSVD,input_all)


Z <- PrepareSummarizedExperimentTensor(matrix(samples$sample,c(6,1)),
rownames(counts),array(counts,c(dim(counts)[1],6)))
HOSVD <- computeHosvd(Z)
cond <- list(0,rep(c("A","B"),each=3))
input_all <- selectSingularValueVectorLarge(HOSVD,cond,input_all=2)

#MultiOmics
require(MOFAdata)
data("CLL_data")
data("CLL_covariates")
Z <- PrepareSummarizedExperimentTensorSquare(
    sample=matrix(colnames(CLL_data$Drugs),1),
    feature=list(Drugs=rownames(CLL_data$Drugs),
    Methylation=rownames(CLL_data$Methylation),
    mRNA=rownames(CLL_data$mRNA),Mutations=rownames(CLL_data$Mutations)),
    value=convertSquare(CLL_data),sampleData=list(CLL_covariates[,1]))
HOSVD <- computeHosvdSqure(Z)
cond <- list(attr(Z,"sampleData")[[1]],attr(Z,"sampleData")[[1]],seq_len(4))
input_all <- selectSingularValueVectorLarge(HOSVD,cond,input_all=c(8,1))
index <- selectFeatureSquare(HOSVD,input_all,CLL_data,
    de=c(0.3,0.03,0.1,0.1),interact=FALSE)
})
