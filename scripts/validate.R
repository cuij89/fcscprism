# Check the web tool's pure-R pipeline against the original CytoML / FlowSOM
# pipeline on the same 1,888 cells.
#
#   Rscript webtool/scripts/validate.R
#
# Needs webtool/testdata (make_testdata.R) and Out/cd8_dp_clusters.csv.

suppressPackageStartupMessages({ library(flowCore) })

root <- normalizePath(file.path(dirname(sub("--file=", "",
  grep("--file=", commandArgs(FALSE), value = TRUE))), "..", ".."))
for (f in list.files(file.path(root, "webtool", "app", "R"), full.names = TRUE))
  source(f)
td <- file.path(root, "webtool", "testdata")

hr <- function(s) cat("\n==== ", s, " ", strrep("=", 60 - nchar(s)), "\n", sep = "")

# ---- 1. FCS reader vs flowCore ----------------------------------------------
hr("1. FCS reader vs flowCore::read.FCS")
files <- c(list.files(td, "\\.fcs$", full.names = TRUE),
           file.path(root, "Data", c("export_T1_Lymphocytes.fcs",
                                     "export_T7_Lymphocytes.fcs")))
for (f in files) {
  t0 <- proc.time()[["elapsed"]]
  a <- fcs_read(f)
  t1 <- proc.time()[["elapsed"]]
  b <- read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)
  eb <- exprs(b)
  same_names <- identical(unname(a$params$name), unname(colnames(eb)))
  d <- max(abs(a$exprs - eb))
  cat(sprintf("  %-32s %8d x %2d  names %s  max|diff| %.3g  (%.2fs)\n",
              basename(f), nrow(a$exprs), ncol(a$exprs),
              if (same_names) "ok " else "DIFF", d, t1 - t0))
}

# ---- 2. load test set the way the app does ----------------------------------
s_tab <- read.csv(file.path(td, "samples.csv"), stringsAsFactors = FALSE)
panel <- read.csv(file.path(td, "panel.csv"), stringsAsFactors = FALSE)
fcs <- lapply(file.path(td, s_tab$file), fcs_read)
names(fcs) <- s_tab$sample
ch <- apply_panel(channel_table(fcs), panel)
sel <- ch[ch$use, ]
cat(sprintf("\n  %d samples, %d cells, %d channels selected\n",
            length(fcs), sum(sapply(fcs, function(f) nrow(f$exprs))), nrow(sel)))

orig <- read.csv(file.path(root, "Out", "cd8_dp_clusters.csv"), check.names = FALSE)
stopifnot(nrow(orig) == sum(sapply(fcs, function(f) nrow(f$exprs))))

# ---- 3. exported raw values -> FlowJo biex reproduces the original matrix ----
hr("2. raw export + FlowJo biex == original pipeline input")
tr <- readRDS(file.path(td, "flowjo_transforms.rds"))
biex <- function(raw) {
  out <- raw
  for (j in seq_len(ncol(raw))) {
    chn <- sel$channel[j]
    f <- tr[[chn]]
    if (is.null(f)) stop("no FlowJo transform for ", chn)
    out[, j] <- f(raw[, j])
  }
  out
}
dat_b <- build_matrix(fcs, sel$channel, sel$marker, cofactor = biex)
mk <- intersect(sel$marker, names(orig))
dev <- max(abs(dat_b$expr[, mk] - as.matrix(orig[, mk])))
cat(sprintf("  max |biex(raw export) - original| over %d markers: %.3g\n", length(mk), dev))
cat(sprintf("  sample order matches: %s\n", identical(dat_b$sample, orig$sample)))

# ---- 4. own SOM + consensus vs FlowSOM + ConsensusClusterPlus -----------------
hr("3. clustering on the identical biex input, k = 6")
Zb <- zscore(dat_b$expr[, mk])
som <- run_som(Zb, dat_b$sample, seed = 42)
cons <- consensus_hc(som$codes, max_k = 6, seed = 42)
cl_b <- cons[["6"]][som$node]
cat(sprintf("  ARI(web tool, original) = %.3f\n", ari(cl_b, orig$cluster)))
print(table(web = cl_b, original = orig$cluster))

# seed sensitivity of the original method itself, for scale
suppressPackageStartupMessages(library(FlowSOM))
fs_run <- function(Z, sample, seed) {
  set.seed(seed)
  per <- min(table(sample))
  train <- unlist(lapply(split(seq_len(nrow(Z)), sample), function(i) sample(i, per)))
  f <- FlowSOM::ReadInput(Z[train, , drop = FALSE], scale = FALSE)
  f <- FlowSOM::BuildSOM(f, colsToUse = seq_len(ncol(Z)), xdim = 8, ydim = 8, silent = TRUE)
  codes <- f$map$codes
  node <- nearest_node(Z, codes)
  cc <- suppressMessages(ConsensusClusterPlus::ConsensusClusterPlus(
    t(codes), maxK = 6, reps = 100, pItem = 0.9, pFeature = 1, title = tempdir(),
    plot = NULL, verbose = FALSE, clusterAlg = "hc", distance = "euclidean", seed = seed))
  cc[[6]]$consensusClass[node]
}
cat("\n  original method, different seeds, ARI vs the seed-42 original:\n")
for (s in c(42, 1, 2, 3)) {
  cat(sprintf("    FlowSOM seed %-3d ARI = %.3f\n", s, ari(fs_run(Zb, dat_b$sample, s), orig$cluster)))
}
cat("  web tool, different seeds:\n")
for (s in c(1, 2, 3)) {
  sm <- run_som(Zb, dat_b$sample, seed = s)
  cl <- consensus_hc(sm$codes, max_k = 6, seed = s)[["6"]][sm$node]
  cat(sprintf("    web seed %-3d ARI = %.3f\n", s, ari(cl, orig$cluster)))
}

# ---- 5. arcsinh instead of biex ---------------------------------------------
hr("4. arcsinh cofactor instead of FlowJo biex, k = 6")
for (cf in c(150, 500, 1000, 3000, 6000)) {
  res <- run_pipeline(fcs, ch, cofactor = cf, max_k = 6)
  cl <- clusters_at(res, 6)
  cat(sprintf("  cofactor %5d  ARI vs original = %.3f   (umap %.1fs, clustering %.1fs)\n",
              cf, ari(cl, orig$cluster), res$timing[["umap"]], res$timing[["cluster"]]))
}
