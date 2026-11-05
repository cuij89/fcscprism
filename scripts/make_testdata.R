# Build test input for the web tool: what a user would get by exporting the
# Q2 (tdTomato+ YFP+) gate of each TIL sample from FlowJo as FCS.
#
# FlowJo exports raw (linear, unmixed) values, not the biex display values the
# original pipeline worked on, so the gated events are inverse-transformed
# before writing.
#
#   Rscript webtool/scripts/make_testdata.R      (needs CytoML / flowWorkspace)

suppressPackageStartupMessages({
  library(CytoML); library(flowWorkspace); library(flowCore)
})

root <- normalizePath(file.path(dirname(sub("--file=", "",
  grep("--file=", commandArgs(FALSE), value = TRUE))), "..", ".."))
data_dir <- file.path(root, "Data")
out_dir <- file.path(root, "webtool", "testdata")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

Q2 <- "/live t cells/CD8s/Q2: Comp-tdTomato-A+ , Comp-YFP-A+"
SAMPLES <- sprintf("export_T%s_Lymphocytes.fcs",
                   c("1", "2", "3", "4", "5-6", "7", "8", "9"))

ws <- open_flowjo_xml(file.path(data_dir, "exported 13-Oct-2026.wsp"))
gs <- flowjo_to_gatingset(ws, name = "TIL", path = data_dir, subset = SAMPLES)

# Two gates per sample, the way a user exporting several populations would:
# the double-positive cells and their parent CD8 gate.
GATES <- c(Q2 = Q2, CD8 = "/live t cells/CD8s")
for (sn in sampleNames(gs)) {
  label <- sub("^export_(.*?)_Lymphocytes\\.fcs.*$", "\\1", sn)
  for (g in names(GATES)) {
    ff <- gh_pop_get_data(gs[[sn]], GATES[[g]], inverse.transform = TRUE)
    ff <- cytoframe_to_flowFrame(ff)
    f <- file.path(out_dir, sprintf("export_%s_%s.fcs", label, g))
    write.FCS(ff, f)
    cat(sprintf("  %-5s %-4s %6d events -> %s\n", label, g, nrow(ff), basename(f)))
  }
}

# The marker names were read off the single-stain controls, not stored in the
# files, so they ship as a panel table - the same thing a user would upload.
markers <- c(
  "Comp-BUV395-A" = "CD8",   "Comp-BUV496-A" = "CD4",
  "Comp-BUV563-A" = "CXCR3", "Comp-BUV615-A" = "CD62L",
  "Comp-BUV661-A" = "CD69",  "Comp-BUV737-A" = "IFNG",
  "Comp-BUV805-A" = "LY108", "Comp-Pacific Blue-A" = "TBET",
  "Comp-BV510-A"  = "Thy1.2", "Comp-BV570-A" = "CD44",
  "Comp-BV605-A"  = "TNFa",  "Comp-BV650-A" = "PD1",
  "Comp-BV711-A"  = "CD127", "Comp-BV785-A" = "CD103",
  "Comp-FITC-A"   = "GZMB",  "Comp-YFP-A"   = "Ifng-YFP",
  "Comp-PerCP-eFluor 710-A" = "TIGIT", "Comp-PE-A" = "IL-2",
  "Comp-tdTomato-A" = "Lag3Cre-tdTomato", "Comp-PE-Dazzle594-A" = "CX3CR1",
  "Comp-PE-Cy5-A" = "CD25",  "Comp-PE-Cy7-A" = "LAG3",
  "Comp-APC-A"    = "KLRG1", "Comp-eFluor 660-A" = "TOX",
  "Comp-Alexa Fluor 700-A" = "CD39", "Comp-Zombie NIR-A" = "Zombie-NIR",
  "Comp-APC-Fire 750-A" = "TIM3"
)
gating <- c("Comp-BUV395-A", "Comp-BUV496-A", "Comp-BV510-A",
            "Comp-Zombie NIR-A", "Comp-tdTomato-A", "Comp-YFP-A")
write.csv(data.frame(channel = names(markers), marker = unname(markers),
                     use = !names(markers) %in% gating),
          file.path(out_dir, "panel.csv"), row.names = FALSE)

labels <- sub("^export_(.*?)_Lymphocytes\\.fcs.*$", "\\1", sampleNames(gs))
write.csv(data.frame(file = sprintf("export_%s_Q2.fcs", labels),
                     sample = labels, group = "unassigned", tissue = "TIL"),
          file.path(out_dir, "samples.csv"), row.names = FALSE)

# Also keep the transform FlowJo used, to check how far arcsinh is from biex.
tr <- gh_get_transformations(gs[[1]])
saveRDS(tr, file.path(out_dir, "flowjo_transforms.rds"))
cat("  wrote panel.csv, samples.csv\n")
