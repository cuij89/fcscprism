# Minimal FCS 3.0 / 3.1 reader in base R.
#
# flowCore has no WebAssembly build, so the web tool cannot use it. This reads
# the HEADER, the TEXT keywords and the DATA segment directly. It covers what
# FlowJo and the Cytek Aurora write: list mode, float / double / fixed-width
# integer data, either byte order.

fcs_read <- function(path, max_events = Inf) {
  con <- file(path, "rb")
  on.exit(close(con))
  size <- file.info(path)$size

  header <- rawToChar(readBin(con, "raw", 58))
  version <- substr(header, 1, 6)
  if (!grepl("^FCS[23]\\.", version))
    stop(basename(path), ": not an FCS file (header '", version, "')")
  off <- function(a, b) suppressWarnings(as.numeric(trimws(substr(header, a, b))))
  text_begin <- off(11, 18); text_end <- off(19, 26)
  data_begin <- off(27, 34); data_end <- off(35, 42)

  seek(con, text_begin)
  txt_raw <- readBin(con, "raw", text_end - text_begin + 1)
  kw <- fcs_parse_text(txt_raw)

  # offsets past 99,999,999 bytes do not fit the header and live in TEXT
  if (is.na(data_begin) || data_begin == 0 || data_end == 0) {
    data_begin <- as.numeric(kw[["$BEGINDATA"]])
    data_end   <- as.numeric(kw[["$ENDDATA"]])
  }

  if (!is.null(kw[["$MODE"]]) && toupper(kw[["$MODE"]]) != "L")
    stop(basename(path), ": only list-mode FCS is supported")

  npar <- as.integer(kw[["$PAR"]])
  ntot <- as.integer(kw[["$TOT"]])
  dtype <- toupper(kw[["$DATATYPE"]])
  byteord <- gsub("\\s", "", kw[["$BYTEORD"]])
  endian <- if (byteord %in% c("1,2,3,4", "1,2")) "little" else "big"

  pkey <- function(i, s) { v <- kw[[sprintf("$P%d%s", i, s)]]; if (is.null(v)) "" else v }
  pn <- vapply(seq_len(npar), pkey, "", s = "N")
  ps <- vapply(seq_len(npar), pkey, "", s = "S")
  pb <- as.integer(vapply(seq_len(npar), pkey, "", s = "B"))

  n_read <- min(ntot, max_events)
  seek(con, data_begin)
  vals <- switch(dtype,
    "F" = readBin(con, "numeric", n = npar * n_read, size = 4, endian = endian),
    "D" = readBin(con, "numeric", n = npar * n_read, size = 8, endian = endian),
    "I" = {
      if (length(unique(pb)) != 1 || !pb[1] %in% c(8, 16, 32))
        stop(basename(path), ": mixed or unusual integer widths are not supported")
      v <- readBin(con, "integer", n = npar * n_read, size = pb[1] / 8,
                   signed = pb[1] == 32, endian = endian)
      if (pb[1] == 32) v[v < 0] <- v[v < 0] + 2^32
      # honour $PnR bit masks for 16-bit data
      as.numeric(v)
    },
    stop(basename(path), ": $DATATYPE '", dtype, "' not supported")
  )
  if (length(vals) != npar * n_read)
    stop(basename(path), ": DATA segment shorter than $TOT x $PAR")

  m <- matrix(vals, ncol = npar, byrow = TRUE)
  colnames(m) <- pn
  list(
    exprs = m,
    params = data.frame(name = pn, desc = ps, stringsAsFactors = FALSE),
    keywords = kw,
    n_total = ntot,
    file = basename(path)
  )
}

# TEXT segment: the first byte is the delimiter; a doubled delimiter inside a
# value is an escaped literal (FCS forbids empty values, so a doubled
# delimiter is never two separators).
fcs_parse_text <- function(raw) {
  raw <- raw[raw != as.raw(0)]
  delim <- rawToChar(raw[1])
  txt <- rawToChar(raw[-1])
  Encoding(txt) <- "latin1"
  txt <- enc2utf8(txt)
  ph <- "\u0001"
  txt <- gsub(paste0(delim, delim), ph, txt, fixed = TRUE)
  tok <- strsplit(txt, delim, fixed = TRUE)[[1]]
  tok <- gsub(ph, delim, tok, fixed = TRUE)
  if (length(tok) %% 2 == 1) tok <- tok[-length(tok)]
  keys <- toupper(trimws(tok[c(TRUE, FALSE)]))
  vals <- tok[c(FALSE, TRUE)]
  stats::setNames(as.list(vals), keys)
}
