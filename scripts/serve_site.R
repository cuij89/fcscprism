# Preview the built site locally. Shinylive needs http://, not file://.
# The port is picked automatically: Windows reserves some port ranges
# (check with: netsh interface ipv4 show excludedportrange protocol=tcp).
#
#   Rscript webtool/scripts/serve_site.R      then open the URL it prints
root <- normalizePath(file.path(dirname(sub("--file=", "",
  grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
port <- httpuv::randomPort(host = "127.0.0.1")
cat(sprintf("Serving webtool/site at http://127.0.0.1:%d  (Ctrl+C to stop)\n", port))
httpuv::runStaticServer(file.path(root, "site"), host = "127.0.0.1", port = port,
                        browse = interactive())
