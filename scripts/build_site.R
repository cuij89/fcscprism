# Bundle webtool/app into a static site (webtool/site) that runs entirely in
# the browser via Shinylive / webR. Upload the contents of site/ to any static
# host; nothing in it contains data.
#
#   Rscript webtool/scripts/build_site.R
root <- normalizePath(file.path(dirname(sub("--file=", "",
  grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
unlink(file.path(root, "site"), recursive = TRUE)
shinylive::export(file.path(root, "app"), file.path(root, "site"))
