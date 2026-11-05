# Analysis core for the web tool. Only uses packages that have WebAssembly
# builds in the webR repository (uwot, kohonen) plus base R, so the same code
# runs in R on a desktop and in the browser.

# ---- input -------------------------------------------------------------------

# Channel table for a set of FCS files: name, FCS description ($PnS), and a
# guess at whether it is a fluorescence parameter worth analysing.
channel_table <- function(fcs_list) {
  # the channel set most files share is the reference, so one odd file (a
  # concatenated file, a control) is named instead of everything else
  keys <- vapply(fcs_list, function(f) paste(f$params$name, collapse = "\r"), "")
  ref <- names(sort(table(keys), decreasing = TRUE))[1]
  p <- fcs_list[[match(ref, keys)]]$params
  same <- keys == ref
  if (!all(same)) {
    msg <- vapply(names(fcs_list)[!same], function(s) {
      nm <- fcs_list[[s]]$params$name
      extra <- setdiff(nm, p$name); missing <- setdiff(p$name, nm)
      show <- function(x) if (length(x) > 3) paste0(paste(x[1:3], collapse = ", "), ", ... (", length(x), ")") else paste(x, collapse = ", ")
      sprintf("%s - %d channels vs %d%s%s", s, length(nm), length(p$name),
              if (length(extra)) paste0("; extra: ", show(extra)) else "",
              if (length(missing)) paste0("; missing: ", show(missing)) else "")
    }, "")
    stop(sprintf("%d of %d files do not match the channels of the other %d: %s. ",
                 sum(!same), length(same), sum(same), paste(msg, collapse = " | ")),
         "Remove them (concatenated files and controls usually differ), ",
         "or re-export every sample with the same FlowJo settings.")
  }
  skip <- grepl("^(FSC|SSC|Time|Event|AF|Original|Sample|Width|Height)",
                p$name, ignore.case = TRUE) |
          grepl("-(H|W)$", p$name)
  data.frame(channel = p$name,
             marker = ifelse(nzchar(p$desc), p$desc, p$name),
             use = !skip,
             stringsAsFactors = FALSE)
}

# FlowJo names exports "export_<sample>_<population>.fcs". Split that into
# sample and population; anything else is one sample of population "all".
parse_export_name <- function(file) {
  stem <- sub("^export_", "", sub("[.]fcs$", "", file, ignore.case = TRUE))
  has <- grepl("_", stem)
  data.frame(sample = ifelse(has, sub("_[^_]*$", "", stem), stem),
             population = ifelse(has, sub("^.*_", "", stem), "all"),
             stringsAsFactors = FALSE)
}

# Merge an uploaded panel table (channel, marker[, use]) into the channel table.
apply_panel <- function(ch, panel) {
  if (is.null(panel)) return(ch)
  names(panel) <- tolower(names(panel))
  i <- match(ch$channel, panel$channel)
  hit <- !is.na(i)
  if ("marker" %in% names(panel)) ch$marker[hit] <- panel$marker[i[hit]]
  if ("use" %in% names(panel))
    ch$use[hit] <- as.logical(panel$use[i[hit]])
  ch
}

asinh_tf <- function(x, cofactor) asinh(x / cofactor)

# Stack the selected channels of every sample, transformed, optionally
# downsampled per sample.
build_matrix <- function(fcs_list, channels, markers, cofactor = 6000,
                         max_per_sample = Inf, seed = 42) {
  set.seed(seed)
  parts <- lapply(names(fcs_list), function(s) {
    m <- fcs_list[[s]]$exprs[, channels, drop = FALSE]
    if (nrow(m) > max_per_sample) m <- m[sort(sample.int(nrow(m), max_per_sample)), , drop = FALSE]
    list(sample = rep(s, nrow(m)), m = m)
  })
  raw <- do.call(rbind, lapply(parts, `[[`, "m"))
  colnames(raw) <- markers
  expr <- if (is.function(cofactor)) cofactor(raw) else asinh_tf(raw, cofactor)
  list(sample = unlist(lapply(parts, `[[`, "sample")), raw = raw, expr = expr)
}

zscore <- function(x) {
  z <- scale(x)
  z[is.na(z)] <- 0
  z
}

# ---- embedding ---------------------------------------------------------------

run_umap <- function(Z, n_neighbors = 15, min_dist = 0.1, seed = 42) {
  set.seed(seed)
  uwot::umap(Z, n_neighbors = min(n_neighbors, nrow(Z) - 1), min_dist = min_dist,
             metric = "euclidean", n_threads = 1, n_sgd_threads = 1)
}

# ---- clustering --------------------------------------------------------------

# SOM trained on an equal number of cells per sample (so one large sample
# cannot claim private nodes), then every cell mapped to its nearest node.
run_som <- function(Z, sample, xdim = 8, ydim = 8, balanced = TRUE, rlen = 10,
                    seed = 42) {
  set.seed(seed)
  idx <- seq_len(nrow(Z))
  train <- if (balanced) {
    per <- min(table(sample))
    unlist(lapply(split(idx, sample), function(i) i[sample.int(length(i), per)]))
  } else idx
  grid <- kohonen::somgrid(xdim, ydim, topo = "rectangular")
  som <- kohonen::som(Z[train, , drop = FALSE], grid = grid, rlen = rlen,
                      alpha = c(0.05, 0.01))
  codes <- som$codes[[1]]
  list(codes = codes, node = nearest_node(Z, codes), n_train = length(train))
}

nearest_node <- function(Z, codes) {
  # squared euclidean distance via the expansion |a|^2 - 2ab + |b|^2
  d <- -2 * Z %*% t(codes)
  d <- sweep(d, 2, rowSums(codes^2), "+")
  max.col(-d, ties.method = "first")
}

# Consensus hierarchical clustering of the SOM nodes, as ConsensusClusterPlus
# does with clusterAlg = "hc": resample pItem of the nodes `reps` times,
# average-linkage cluster each subsample, and record how often each pair of
# nodes lands together. The final partition clusters 1 - consensus.
# Returns one node -> metacluster vector per k in 2..max_k.
consensus_hc <- function(codes, max_k = 15, reps = 100, p_item = 0.9, seed = 42) {
  set.seed(seed)
  n <- nrow(codes)
  max_k <- min(max_k, n - 1)
  ks <- 2:max_k
  M <- lapply(ks, function(k) matrix(0, n, n))
  I <- matrix(0, n, n)
  D <- as.matrix(dist(codes))
  for (r in seq_len(reps)) {
    s <- sort(sample.int(n, floor(n * p_item)))
    hc <- hclust(as.dist(D[s, s]), method = "average")
    I[s, s] <- I[s, s] + 1
    for (j in seq_along(ks)) {
      cl <- cutree(hc, ks[j])
      M[[j]][s, s] <- M[[j]][s, s] + outer(cl, cl, "==")
    }
  }
  out <- lapply(seq_along(ks), function(j) {
    C <- M[[j]] / pmax(I, 1)
    diag(C) <- 1
    hc <- hclust(as.dist(1 - C), method = "average")
    relabel_by_size(cutree(hc, ks[j]))
  })
  names(out) <- ks
  out
}

# Number clusters 1..k from largest to smallest, so labels are stable to read.
relabel_by_size <- function(cl) {
  o <- order(-tabulate(cl))
  match(cl, o)
}

# ---- description -------------------------------------------------------------

# Per-cluster median of each marker (transformed scale), scored against ALL
# cells: (cluster median - overall median) / overall SD, plus AUC of the
# cluster vs every other cell.
cluster_markers <- function(expr, cluster) {
  centre <- apply(expr, 2, median)
  spread <- apply(expr, 2, sd)
  spread[spread == 0] <- 1
  cl <- sort(unique(cluster))
  med <- t(sapply(cl, function(k) apply(expr[cluster == k, , drop = FALSE], 2, median)))
  if (length(cl) == 1) med <- matrix(med, nrow = 1, dimnames = list(NULL, colnames(expr)))
  rownames(med) <- cl
  z <- sweep(sweep(med, 2, centre), 2, spread, "/")
  auc <- t(sapply(cl, function(k) {
    inn <- cluster == k
    apply(expr, 2, function(v) {
      r <- rank(v)
      n1 <- sum(inn); n0 <- sum(!inn)
      if (n0 == 0) return(NA_real_)
      (sum(r[inn]) - n1 * (n1 + 1) / 2) / (n1 * n0)
    })
  }))
  if (length(cl) == 1) auc <- matrix(auc, nrow = 1, dimnames = list(NULL, colnames(expr)))
  rownames(auc) <- cl
  list(median = med, z = z, auc = auc)
}

# Top markers per cluster in words: highest and lowest AUC.
cluster_summary <- function(cm, cluster, n = 4) {
  cl <- rownames(cm$auc)
  sizes <- table(factor(cluster, levels = cl))
  data.frame(
    cluster = cl,
    cells = as.integer(sizes),
    percent = round(100 * as.numeric(sizes) / length(cluster), 1),
    high = vapply(cl, function(k) {
      a <- sort(cm$auc[k, ], decreasing = TRUE)[seq_len(n)]
      paste(sprintf("%s (%.2f)", names(a), a), collapse = ", ")
    }, ""),
    low = vapply(cl, function(k) {
      a <- sort(cm$auc[k, ])[seq_len(n)]
      paste(sprintf("%s (%.2f)", names(a), a), collapse = ", ")
    }, ""),
    stringsAsFactors = FALSE, row.names = NULL
  )
}

# Cluster share within each sample, then mean and SEM across samples per group.
# The per-sample mean is the number to report: pooling lets one large sample
# set the answer.
cluster_proportions <- function(sample, cluster, groups = NULL, exclude = NULL) {
  keep <- !cluster %in% exclude
  tab <- table(sample = sample[keep], cluster = droplevels(as.factor(cluster[keep])))
  pct <- 100 * prop.table(tab, 1)
  per_sample <- as.data.frame.matrix(round(pct, 2))
  per_sample <- cbind(sample = rownames(per_sample), n = as.integer(rowSums(tab)),
                      per_sample)
  if (is.null(groups)) groups <- setNames(rep("all", nrow(pct)), rownames(pct))
  g <- groups[rownames(pct)]
  g[is.na(g)] <- "unassigned"
  summ <- do.call(rbind, lapply(split(seq_len(nrow(pct)), g), function(i) {
    x <- pct[i, , drop = FALSE]
    data.frame(group = g[i[1]], cluster = colnames(pct), n_samples = length(i),
               mean = colMeans(x),
               sem = if (length(i) > 1) apply(x, 2, sd) / sqrt(length(i)) else NA_real_,
               pooled = as.numeric(100 * colSums(tab[i, , drop = FALSE]) / sum(tab[i, ])),
               row.names = NULL)
  }))
  list(per_sample = per_sample, summary = summ)
}

# Adjusted Rand index between two partitions.
ari <- function(a, b) {
  t <- table(a, b)
  c2 <- function(x) x * (x - 1) / 2
  s <- sum(c2(t)); sa <- sum(c2(rowSums(t))); sb <- sum(c2(colSums(t)))
  e <- sa * sb / c2(sum(t))
  (s - e) / ((sa + sb) / 2 - e)
}

# ---- one call for the whole run ---------------------------------------------

run_pipeline <- function(fcs_list, ch, cofactor = 6000, max_per_sample = Inf,
                         n_neighbors = 15, min_dist = 0.1, xdim = 8, ydim = 8,
                         balanced = TRUE, max_k = 15, k = 6, seed = 42,
                         progress = function(msg) NULL) {
  sel <- ch[ch$use, ]
  progress("Building matrix")
  dat <- build_matrix(fcs_list, sel$channel, sel$marker, cofactor, max_per_sample, seed)
  Z <- zscore(dat$expr)
  t0 <- proc.time()[["elapsed"]]
  progress("UMAP")
  emb <- run_umap(Z, n_neighbors, min_dist, seed)
  t1 <- proc.time()[["elapsed"]]
  progress("SOM")
  som <- run_som(Z, dat$sample, xdim, ydim, balanced, seed = seed)
  progress("Consensus clustering")
  cons <- consensus_hc(som$codes, max_k = max_k, seed = seed)
  t2 <- proc.time()[["elapsed"]]
  list(sample = dat$sample, raw = dat$raw, expr = dat$expr, Z = Z,
       umap = emb, som = som, consensus = cons,
       params = list(cofactor = cofactor, max_per_sample = max_per_sample,
                     n_neighbors = n_neighbors, min_dist = min_dist,
                     xdim = xdim, ydim = ydim, balanced = balanced,
                     max_k = max_k, seed = seed, channels = sel$channel,
                     markers = sel$marker),
       timing = c(umap = t1 - t0, cluster = t2 - t1))
}

clusters_at <- function(res, k) res$consensus[[as.character(k)]][res$som$node]
