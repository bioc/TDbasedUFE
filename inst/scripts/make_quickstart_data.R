# Run from the package root.
# Requires the original six-sample tximportData example.

stopifnot(
    requireNamespace("tximportData", quietly = TRUE),
    requireNamespace("tximport", quietly = TRUE)
)

d <- system.file("extdata", package = "tximportData")
samples <- read.table(
    file.path(d, "samples.txt"),
    header = TRUE
)
stopifnot(nrow(samples) == 6L)

samples$condition <- factor(rep(c("A", "B"), each = 3))
rownames(samples) <- samples$run

files <- file.path(d, "salmon", samples$run, "quant.sf.gz")
names(files) <- samples$run
stopifnot(all(file.exists(files)))

tx2gene <- read.csv(
    file.path(d, "tx2gene.gencode.v27.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
)

txi <- tximport::tximport(
    files,
    type = "salmon",
    tx2gene = tx2gene,
    importer = function(x) {
        read.delim(
            x,
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
    }
)

stopifnot(
    nrow(txi$counts) >= 10000L,
    ncol(txi$counts) == 6L,
    identical(rownames(txi$counts), rownames(txi$abundance))
)

counts <- txi$counts[seq_len(10000), , drop = FALSE]
stopifnot(identical(colnames(counts), samples$run))

example <- list(
    counts = counts,
    samples = samples,
    source = list(
        package = "tximportData",
        version = as.character(packageVersion("tximportData")),
        tximport_version = as.character(packageVersion("tximport")),
        original_dimension = dim(txi$counts),
        selected_rows = c(1L, 10000L),
        countsFromAbundance = txi$countsFromAbundance
    )
)

output <- "inst/extdata/quickstart_counts.rds"
saveRDS(example, output, compress = "xz")

reloaded <- readRDS(output)
stopifnot(
    identical(example, reloaded),
    identical(dim(reloaded$counts), c(10000L, 6L))
)

cat("Original matrix:", dim(txi$counts), "\n")
cat("Saved matrix:", dim(reloaded$counts), "\n")
cat("Saved file:", output, "\n")
cat("File size (bytes):", file.info(output)$size, "\n")
cat("Exact reload verification: OK\n")
cat("tximportData license:",
    packageDescription("tximportData")$License, "\n")
