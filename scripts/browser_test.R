# End-to-end test of the built site in headless Chrome: loads it, uploads
# testdata/, runs the analysis and saves screenshots to test_shots/.
#
#   Rscript webtool/scripts/browser_test.R     (needs chromote, Chrome)

library(chromote)
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
site <- file.path(root, "site"); td <- file.path(root, "testdata")
port <- httpuv::randomPort(host = "127.0.0.1")   # Windows reserves some port ranges
srv <- httpuv::startServer("127.0.0.1", port, list(staticPaths = list("/" = httpuv::staticPath(site, indexhtml = TRUE))))
shot_dir <- file.path(root, "test_shots"); dir.create(shot_dir, showWarnings = FALSE)

b <- ChromoteSession$new(width = 1400, height = 1000)
t0 <- Sys.time()
b$Page$navigate(sprintf("http://127.0.0.1:%d/", port))
js <- function(code) b$Runtime$evaluate(code, returnByValue = TRUE)$result$value
doc <- "(document.querySelector('iframe') ? document.querySelector('iframe').contentDocument : document)"
wait_for <- function(cond, timeout = 300, what = cond) {
  t <- Sys.time()
  repeat {
    v <- tryCatch(js(cond), error = function(e) FALSE)
    if (isTRUE(v)) return(as.numeric(difftime(Sys.time(), t, units = "secs")))
    if (difftime(Sys.time(), t, units = "secs") > timeout) stop("timeout waiting for: ", what)
    later::run_now(0.5)
  }
}
cat(sprintf("wasm mime: %s\n", js("fetch('shinylive/webr/R.wasm',{method:'HEAD'}).then(r=>r.headers.get('content-type'))") %||% "?"))
s <- wait_for(sprintf("!!%s && !!%s.querySelector('#fcs')", doc, doc), what = "app UI")
cat(sprintf("app ready after %.1f s (first load, nothing cached)\n", s))
b$screenshot(file.path(shot_dir, "01_loaded.png"))

# upload files into the fileInputs inside the iframe
dom <- b$DOM$getDocument(depth = -1, pierce = TRUE)
frame_doc <- function() {
  d <- b$DOM$getDocument(depth = -1, pierce = TRUE)$root
  find <- function(n) { if (!is.null(n$contentDocument)) return(n$contentDocument$nodeId)
    for (c in n$children) { r <- find(c); if (!is.null(r)) return(r) }; NULL }
  find(d) %||% d$nodeId
}
upload <- function(sel, files) {
  fd <- frame_doc()
  nid <- b$DOM$querySelector(fd, sel)$nodeId
  b$DOM$setFileInputFiles(files = as.list(normalizePath(files, winslash = "/")), nodeId = nid)
}
fcs <- list.files(td, "[.]fcs$", full.names = TRUE)
upload("#fcs", fcs)
s <- wait_for(sprintf("%s.body.innerText.includes('Read 8 files')", doc), what = "FCS read")
cat(sprintf("8 FCS files read in %.1f s\n", s))
upload("#panel_csv", file.path(td, "panel.csv"))
upload("#samples_csv", file.path(td, "samples.csv"))
wait_for(sprintf("!!%s.querySelector('#sname_1') && %s.querySelector('#sname_1').value=='T1'", doc, doc), what = "sample table")
later::run_now(2)
b$screenshot(file.path(shot_dir, "02_data.png"))
js(sprintf("%s.querySelector('a[data-value=\"2 · Channels\"]').click()", doc)); later::run_now(2)
wait_for(sprintf("%s.querySelectorAll('input[id^=use_]').length > 0", doc), what = "channel table")
later::run_now(1)
n_used <- js(sprintf("[...%s.querySelectorAll('input[id^=use_]')].filter(x=>x.checked).length", doc))
cat(sprintf("channels ticked after panel.csv: %s\n", n_used))
b$screenshot(file.path(shot_dir, "03_channels.png"))

js(sprintf("%s.querySelector('a[data-value=\"3 · Run\"]').click()", doc)); later::run_now(1)
js(sprintf("%s.querySelector('#run').click()", doc))
s <- wait_for(sprintf("%s.querySelector('#run_log') && %s.querySelector('#run_log').innerText.startsWith('Done')", doc, doc), timeout = 600, what = "run")
cat(sprintf("analysis finished after %.1f s (wall clock)\n", s))
cat(js(sprintf("%s.querySelector('#run_log').innerText", doc)), "\n")
later::run_now(3)
b$screenshot(file.path(shot_dir, "04_umap.png"))
for (tab in c("Markers", "Proportions")) {
  js(sprintf("[...%s.querySelectorAll('a.nav-link')].find(a=>a.innerText.trim()=='%s').click()", doc, tab))
  later::run_now(4)
  b$screenshot(file.path(shot_dir, sprintf("05_%s.png", tolower(tab))))
}
cat(sprintf("total %.1f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
err <- js(sprintf("[...%s.querySelectorAll('.shiny-output-error')].map(e=>e.innerText).join(' | ')", doc))
cat("output errors:", if (nzchar(err)) err else "none", "\n")
b$close()
httpuv::stopServer(srv)
