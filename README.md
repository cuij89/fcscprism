# fcscprism — browser-based UMAP and clustering for flow cytometry (prototype)

A Shiny app that runs entirely in the browser through Shinylive / webR. Users
pick their own FCS files, and the files are read and analysed on their own
computer. Nothing is uploaded, and nothing is stored on the server.

**Gating is not done here.** Each user gates in FlowJo and exports the
population to analyse as one FCS file per sample. The tool starts from those
files.

## Layout

```
app/app.R            Shiny UI and server
app/R/fcs.R          FCS 3.0/3.1 reader in base R (flowCore has no WebAssembly build)
app/R/pipeline.R     arcsinh → z-score → UMAP (uwot) → SOM (kohonen) → consensus
                     hierarchical metaclustering (own implementation, same algorithm
                     as ConsensusClusterPlus clusterAlg = "hc") → markers, proportions
scripts/make_testdata.R   build testdata/ from Data/ with CytoML (Q2 gate, raw values)
scripts/validate.R        compare against the original CytoML/FlowSOM pipeline
scripts/build_site.R      shinylive::export → site/
scripts/serve_site.R      preview site/ locally (prints the URL)
scripts/browser_test.R    headless-Chrome end-to-end test of site/
testdata/  site/  test_shots/   generated, git-ignored
```

## Use

```bash
Rscript webtool/scripts/make_testdata.R   # once; needs CytoML
Rscript -e 'shiny::runApp("webtool/app")' # desktop R, fastest way to develop
Rscript webtool/scripts/build_site.R      # bundle for the browser (~73 MB)
Rscript webtool/scripts/serve_site.R      # then open the URL it prints
```

Shinylive pages must be served over http(s). Opening `site/index.html` from
disk does not work.

## What the user does

1. In FlowJo, export the gate to analyse, one FCS per sample, with the same
   export settings for every sample.
2. Open the tool and select the FCS files. Optionally also select
   `samples.csv` (`file, sample, group`) and `panel.csv`
   (`channel, marker, use`).
3. Check the event counts against FlowJo. Untick the channels used for gating,
   then run.
4. Choose k, name the clusters, mark artefacts to exclude, and download the
   tables.

## Validation on the CD8 DP data (1,888 cells)

Run with `scripts/validate.R`:

- **FCS reader**: identical to `flowCore::read.FCS` (max |diff| = 0) on the 8
  test exports and on the full 39k- and 85k-event Aurora files.
- **Export path**: reapplying FlowJo's biex to the exported raw values gives
  back the original input. The largest deviation is 2.2, against a median
  marker SD of 53, which is float32 round-trip error.
- **Clustering**: on identical input at k = 6, ARI against the original is
  0.59. For comparison, the original FlowSOM method rerun here gives ARI 0.28–0.45
  against the saved result across seeds 1, 2, 3 and 42. The k = 6 partition is
  seed-sensitive in either implementation, so a seed/stability readout is worth
  adding.
- **Transform**: arcsinh with cofactor 1000–6000 gives ARI 0.43–0.49
  against the biex original, which is inside that seed noise.
- Timing for 1,888 cells: on desktop R, UMAP takes 5–6 s and clustering
  0.2 s. In the browser (headless Chrome, `scripts/browser_test.R`), first
  load with nothing cached takes ~10 s, reading the 8 FCS files 0.4 s, UMAP
  7.3 s and clustering 0.5 s, with no output errors.

## Hosting

`site/` is static files only, with no data in them. Wix cannot host them, so
put them on a static host (GitHub Pages, Netlify, Cloudflare Pages, Pitt web
space) and link to it from the Wix page in a new tab.
