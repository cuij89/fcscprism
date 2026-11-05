# fcscprism — browser-based UMAP and clustering for flow cytometry (prototype)

A Shiny app that runs entirely in the browser through Shinylive / webR. Users
pick their own FCS files, and the files are read and analysed on their own
computer. Nothing is uploaded, and nothing is stored on the server.

**Gating is not done here.** Each user gates in FlowJo and exports the
population to analyse as one FCS file per sample. The tool starts from those
files.

Live (temporary, public): https://cuij89.github.io/fcscprism/

## Layout

```
app/app.R            Shiny UI and server
app/R/fcs.R          FCS 3.0/3.1 reader in base R (flowCore has no WebAssembly build)
app/R/pipeline.R     arcsinh → z-score → UMAP (uwot) → SOM (kohonen) → consensus
                     hierarchical metaclustering (own implementation, same algorithm
                     as ConsensusClusterPlus clusterAlg = "hc") → markers, proportions
scripts/make_testdata.R   build testdata/ from the analysis project's Data/ with CytoML
                          (Q2 and CD8 gates per sample, raw values)
scripts/validate.R        compare against the original CytoML/FlowSOM pipeline
scripts/build_site.R      shinylive::export → site/
scripts/serve_site.R      preview site/ locally (prints the URL)
scripts/browser_test.R    headless-Chrome end-to-end test of site/
testdata/  site/  test_shots/   generated, git-ignored
```

## Use

Run from the repository root. The repo sits in the analysis project as
`webtool/`; `make_testdata.R` and `validate.R` read `../Data` and `../Out`.

```bash
Rscript scripts/make_testdata.R           # once; needs CytoML
Rscript -e 'shiny::runApp("app")'         # desktop R, fastest way to develop
Rscript scripts/build_site.R              # bundle for the browser (~73 MB)
Rscript scripts/serve_site.R              # then open the URL it prints
Rscript scripts/browser_test.R            # headless end-to-end test, screenshots in test_shots/
```

Shinylive pages must be served over http(s). Opening `site/index.html` from
disk does not work.

## What the user does

1. In FlowJo, export the gate to analyse, one FCS per sample, with the same
   export settings for every sample. Several gates can be exported and
   loaded together (for example Q2 and CD8 for every sample).
2. Open the tool and select the FCS files. Optionally also select
   `samples.csv` (`file, sample, population, group`) and `panel.csv`
   (`channel, marker, use`). Sample and population are read from FlowJo's
   `export_<sample>_<gate>.fcs` names and can be edited. The group (for
   example control or treated) is set once per sample, and proportions are
   averaged per group.
3. Choose which population to analyse.
4. Check the event counts against FlowJo. Untick the channels used for gating,
   then run.
5. Choose k, name the clusters, mark artefacts to exclude, and download the
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

## Size and speed in the browser

- Uploads are capped at 1 GB (`shiny.maxRequestSize`; Shiny's default is
  5 MB). Files are only read into the browser tab's memory. webR is 32-bit
  WebAssembly, so the tab has under 4 GB in total.
- Reading files is fast: 8 Lymphocytes files (53.7 MB, 454k events) take
  ~6 s.
- UMAP is the bottleneck. It costs about 20 s per 10,000 cells, measured in
  headless Chrome and single-threaded, so it depends on the user's CPU.

  | cells per sample | total cells | UMAP | clustering |
  |---|---|---|---|
  | all (Q2) | 1,888 | 7 s | 0.5 s |
  | 5,000 | 40,000 | 68 s | 3 s |
  | 20,000 | 160,000 | 298 s | 10 s |

- The default is 5,000 cells per sample. The Run page shows an estimate,
  which turns red above 3 minutes.

## Hosting

`site/` is static files only, with no data in them. Wix cannot host them, so
put them on a static host (GitHub Pages, Netlify, Cloudflare Pages, Pitt web
space) and link to it from the Wix page in a new tab.

The current deployment pushes the local build to the `gh-pages` branch of
this repo. GitHub Pages serves from that branch, and `.nojekyll` must be
present:

```bash
Rscript scripts/build_site.R
# copy site/ (plus an empty .nojekyll) into a checkout of gh-pages, commit, push
```

GitHub Pages on a free account needs a public repository, and a Pages site
is public even when the repository is private, except on Enterprise. To take
it down:

```bash
gh api -X DELETE repos/cuij89/fcscprism/pages
git push origin --delete gh-pages
gh repo edit cuij89/fcscprism --visibility private --accept-visibility-change-consequences
```
