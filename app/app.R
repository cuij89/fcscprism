# Flow cytometry clustering web tool.
#
# Input: one FCS file per sample, each already gated in FlowJo to the
# population of interest and exported. No gating happens here.
# Runs in R on a desktop (shiny::runApp("webtool/app")) or entirely in the
# browser through Shinylive / webR - files never leave the user's machine.
# Helpers in R/ are sourced automatically by Shiny.

library(shiny)
library(bslib)

# Shiny caps uploads at 5 MB by default. In the browser build nothing is sent
# to a server - the file is only copied into webR's memory - so the real limit
# is the browser tab's memory (webR is 32-bit WebAssembly, < 4 GB in total).
options(shiny.maxRequestSize = 1024^3)

PAL <- c("#2a78d6", "#e8743b", "#19a979", "#d64e8a", "#945ecf", "#13a4b4",
         "#b8860b", "#525df4", "#bf399e", "#6c8893", "#ee6868", "#2f6497",
         "#6db33f", "#a0522d", "#7a7a7a")
RAMP <- c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")
DIV <- colorRampPalette(c("#2f6497", "#9ec5f4", "#f4f4f2", "#f2b28c", "#c4512b"))(41)

read_table_upload <- function(input_file) {
  if (is.null(input_file)) return(NULL)
  read.csv(input_file$datapath, stringsAsFactors = FALSE, check.names = FALSE)
}

safe_id <- function(prefix, i) paste0(prefix, "_", i)

ui <- page_navbar(
  id = "main",
  title = "Flow clustering",
  theme = bs_theme(version = 5, primary = "#2a78d6", base_font = font_collection("system-ui", "sans-serif")),
  fillable = FALSE,
  header = tags$style(".cl-name .form-group, .cl-excl .form-group{margin-bottom:0} .cl-excl .checkbox{margin:0}"),

  nav_panel("1 · Data",
    layout_columns(col_widths = c(5, 7),
      card(
        card_header("Upload"),
        p(class = "text-muted small",
          "One FCS file per sample, exported from FlowJo at the gate you want to analyse",
          "(right-click the gate → Export / Concatenate Populations → FCS, all compensated parameters).",
          "Files are read in your browser and are not uploaded anywhere."),
        fileInput("fcs", "FCS files", multiple = TRUE, accept = ".fcs"),
        fileInput("samples_csv", "Sample table (optional CSV: file, sample, group)", accept = ".csv"),
        fileInput("panel_csv", "Panel table (optional CSV: channel, marker, use)", accept = ".csv"),
        uiOutput("load_msg")
      ),
      card(
        card_header("Samples"),
        p(class = "text-muted small", "Check the event counts against FlowJo. Edit names and groups here."),
        uiOutput("sample_editor")
      )
    )
  ),

  nav_panel("2 · Channels",
    card(
      card_header("Channels used for UMAP and clustering"),
      p(class = "text-muted small",
        "Untick scatter, time, viability and the channels you gated on - they only tell the",
        "map what you already selected. Marker names are taken from the panel table, the FCS",
        "$PnS keyword, or the channel name, in that order."),
      div(class = "d-flex gap-2 mb-2",
        actionButton("ch_all", "Select all", class = "btn-sm btn-outline-secondary"),
        actionButton("ch_none", "Select none", class = "btn-sm btn-outline-secondary"),
        downloadButton("dl_panel", "Save as panel.csv", class = "btn-sm btn-outline-secondary")),
      uiOutput("channel_editor")
    )
  ),

  nav_panel("3 · Run",
    layout_columns(col_widths = c(4, 8),
      card(
        card_header("Parameters"),
        numericInput("cofactor", "arcsinh cofactor", 6000, min = 1, step = 500),
        numericInput("max_per", "Max cells per sample (downsample)", 5000, min = 50, step = 1000),
        numericInput("nn", "UMAP n_neighbors", 15, min = 2, max = 200),
        numericInput("md", "UMAP min_dist", 0.1, min = 0, max = 1, step = 0.05),
        layout_columns(
          numericInput("xdim", "SOM x", 8, min = 2, max = 20),
          numericInput("ydim", "SOM y", 8, min = 2, max = 20)),
        checkboxInput("balanced", "Train SOM on equal cells per sample", TRUE),
        numericInput("maxk", "Max k to compute", 15, min = 2, max = 30),
        numericInput("seed", "Random seed", 42),
        actionButton("run", "Run analysis", class = "btn-primary w-100"),
        uiOutput("run_estimate")
      ),
      card(
        card_header("Status"),
        verbatimTextOutput("run_log")
      )
    )
  ),

  nav_panel("4 · Results",
    layout_sidebar(
      sidebar = sidebar(width = 280,
        sliderInput("k", "Number of clusters (k)", min = 2, max = 15, value = 6, step = 1),
        selectInput("colour_by", "Colour UMAP by", choices = c("cluster", "sample")),
        uiOutput("cluster_editor")
      ),
      navset_card_underline(
        nav_panel("UMAP", plotOutput("umap_plot", height = "620px")),
        nav_panel("By sample", plotOutput("facet_plot", height = "700px")),
        nav_panel("Markers", plotOutput("heatmap", height = "560px"),
                  tableOutput("cluster_table")),
        nav_panel("Proportions",
                  plotOutput("prop_plot", height = "460px"),
                  p(class = "text-muted small",
                    "Bars are the mean across samples (± SEM) - not pooled cells, so one large sample cannot set the answer.",
                    "Clusters marked 'exclude' are dropped before the percentages are computed."),
                  tableOutput("prop_table")),
        nav_panel("Download",
          p("Everything is computed in this browser session. Download what you need before closing the tab."),
          downloadButton("dl_cells", "Per-cell table (UMAP, cluster, markers)"), br(), br(),
          downloadButton("dl_markers", "Cluster marker table"), br(), br(),
          downloadButton("dl_props", "Cluster proportions per sample"), br(), br(),
          downloadButton("dl_params", "Run parameters"))
      )
    )
  ),
  nav_spacer(),
  nav_item(tags$span(class = "text-muted small", "runs locally · no data leaves this computer"))
)

server <- function(input, output, session) {
  rv <- reactiveValues(fcs = NULL, ch = NULL, res = NULL, log = "No run yet.")

  # ---- data -----------------------------------------------------------------
  observeEvent(input$fcs, {
    f <- input$fcs
    out <- tryCatch({
      lst <- lapply(seq_len(nrow(f)), function(i) {
        x <- fcs_read(f$datapath[i]); x$file <- f$name[i]; x
      })
      names(lst) <- sub("\\.fcs$", "", f$name, ignore.case = TRUE)
      ch <- channel_table(lst)
      list(lst = lst, ch = ch)
    }, error = function(e) e)
    if (inherits(out, "error")) {
      rv$fcs <- NULL; rv$ch <- NULL
      output$load_msg <- renderUI(div(class = "alert alert-danger", conditionMessage(out)))
      return()
    }
    rv$fcs <- out$lst
    rv$ch <- apply_panel(out$ch, read_table_upload(input$panel_csv))
    rv$res <- NULL
    output$load_msg <- renderUI(div(class = "alert alert-success",
      sprintf("Read %d files, %s events, %d channels.", length(out$lst),
              format(sum(sapply(out$lst, function(x) nrow(x$exprs))), big.mark = ","),
              nrow(out$ch))))
  })

  observeEvent(input$panel_csv, {
    req(rv$ch)
    rv$ch <- apply_panel(rv$ch, read_table_upload(input$panel_csv))
  })

  sample_meta <- reactive({
    req(rv$fcs)
    files <- vapply(rv$fcs, `[[`, "", "file")
    st <- read_table_upload(input$samples_csv)
    sample <- names(rv$fcs); group <- rep("all", length(files))
    if (!is.null(st)) {
      names(st) <- tolower(names(st))
      i <- match(files, st$file)
      hit <- !is.na(i)
      if ("sample" %in% names(st)) sample[hit] <- st$sample[i[hit]]
      if ("group" %in% names(st)) group[hit] <- st$group[i[hit]]
    }
    data.frame(file = files, sample = sample, group = group,
               events = vapply(rv$fcs, function(x) nrow(x$exprs), 1L),
               stringsAsFactors = FALSE)
  })

  output$sample_editor <- renderUI({
    m <- sample_meta()
    rows <- lapply(seq_len(nrow(m)), function(i) tags$tr(
      tags$td(class = "small", m$file[i]),
      tags$td(class = "text-end small", format(m$events[i], big.mark = ",")),
      tags$td(textInput(safe_id("sname", i), NULL, m$sample[i], width = "120px")),
      tags$td(textInput(safe_id("sgroup", i), NULL, m$group[i], width = "120px"))))
    tags$table(class = "table table-sm align-middle",
      tags$thead(tags$tr(tags$th("File"), tags$th(class = "text-end", "Events"),
                         tags$th("Sample"), tags$th("Group"))),
      tags$tbody(rows))
  })

  samples_now <- reactive({
    m <- sample_meta()
    for (i in seq_len(nrow(m))) {
      v <- input[[safe_id("sname", i)]]; if (!is.null(v) && nzchar(v)) m$sample[i] <- v
      g <- input[[safe_id("sgroup", i)]]; if (!is.null(g) && nzchar(g)) m$group[i] <- g
    }
    m
  })

  # ---- channels ---------------------------------------------------------------
  output$channel_editor <- renderUI({
    ch <- rv$ch
    if (is.null(ch)) return(p(class = "text-muted", "Load FCS files first."))
    rows <- lapply(seq_len(nrow(ch)), function(i) tags$tr(
      tags$td(checkboxInput(safe_id("use", i), NULL, ch$use[i])),
      tags$td(class = "small", ch$channel[i]),
      tags$td(textInput(safe_id("mk", i), NULL, ch$marker[i], width = "180px"))))
    tags$table(class = "table table-sm align-middle",
      tags$thead(tags$tr(tags$th("Use"), tags$th("Channel"), tags$th("Marker"))),
      tags$tbody(rows))
  })
  observeEvent(input$ch_all, { req(rv$ch); for (i in seq_len(nrow(rv$ch))) updateCheckboxInput(session, safe_id("use", i), value = TRUE) })
  observeEvent(input$ch_none, { req(rv$ch); for (i in seq_len(nrow(rv$ch))) updateCheckboxInput(session, safe_id("use", i), value = FALSE) })

  channels_now <- reactive({
    ch <- rv$ch; req(ch)
    for (i in seq_len(nrow(ch))) {
      u <- input[[safe_id("use", i)]]; if (!is.null(u)) ch$use[i] <- u
      m <- input[[safe_id("mk", i)]]; if (!is.null(m) && nzchar(m)) ch$marker[i] <- m
    }
    ch$marker <- make.unique(ch$marker)
    ch
  })

  output$dl_panel <- downloadHandler("panel.csv", function(file)
    write.csv(channels_now(), file, row.names = FALSE))

  # ---- run ----------------------------------------------------------------------
  # Measured in Chrome via webR: UMAP ~20 s per 10,000 cells (7 s for 1,888,
  # 68 s for 40,000, 298 s for 160,000); everything else is small.
  run_size <- reactive({
    ev <- sample_meta()$events
    mp <- input$max_per
    if (is.null(mp) || is.na(mp)) mp <- Inf
    n <- sum(pmin(ev, mp))
    sec <- 5 + n / 10000 * 20
    label <- if (sec < 60) sprintf("%.0f s", max(5, round(sec, -1))) else sprintf("%.0f min", ceiling(sec / 60))
    list(n = n, sec = sec, label = label)
  })
  output$run_estimate <- renderUI({
    req(rv$fcs)
    e <- run_size()
    div(class = paste("small mt-2", if (e$sec > 180) "text-danger" else "text-muted"),
        sprintf("%s cells will be analysed - about %s in the browser.",
                format(e$n, big.mark = ","), e$label),
        if (e$sec > 180) " Lower 'Max cells per sample' to explore parameters first.")
  })

  observeEvent(input$run, {
    req(rv$fcs)
    ch <- channels_now()
    if (sum(ch$use) < 2) { rv$log <- "Select at least two channels."; return() }
    m <- samples_now()
    fcs <- rv$fcs; names(fcs) <- make.unique(m$sample)
    t0 <- proc.time()[["elapsed"]]
    est <- run_size()
    withProgress(message = "Running", value = 0, {
      res <- tryCatch(run_pipeline(fcs, ch,
        cofactor = input$cofactor, max_per_sample = input$max_per,
        n_neighbors = input$nn, min_dist = input$md,
        xdim = input$xdim, ydim = input$ydim, balanced = input$balanced,
        max_k = input$maxk, seed = input$seed,
        progress = function(msg) incProgress(0.25, detail = if (msg == "UMAP")
          sprintf("UMAP on %s cells - about %s, no progress is shown until it finishes",
                  format(est$n, big.mark = ","), est$label) else msg)),
        error = function(e) e)
    })
    if (inherits(res, "error")) { rv$log <- paste("Error:", conditionMessage(res)); return() }
    res$groups <- setNames(m$group, names(fcs))
    rv$res <- res
    kmax <- length(res$consensus) + 1
    updateSliderInput(session, "k", max = kmax, value = min(isolate(input$k), kmax))
    updateSelectInput(session, "colour_by", choices = c("cluster", "sample", res$params$markers))
    rv$log <- sprintf(paste0("Done in %.1f s  (UMAP %.1f s, clustering %.1f s)\n",
                             "%s cells from %d samples, %d markers\n",
                             "SOM %dx%d trained on %s cells%s\n\nGo to 4 · Results."),
      proc.time()[["elapsed"]] - t0, res$timing[["umap"]], res$timing[["cluster"]],
      format(length(res$sample), big.mark = ","), length(unique(res$sample)),
      ncol(res$expr), res$params$xdim, res$params$ydim,
      format(res$som$n_train, big.mark = ","),
      if (res$params$balanced) " (balanced per sample)" else "")
    nav_select("main", "4 · Results")
  })
  output$run_log <- renderText(rv$log)

  # ---- results ------------------------------------------------------------------
  cl <- reactive({ req(rv$res); clusters_at(rv$res, input$k) })

  cluster_names <- reactive({
    k <- input$k
    nm <- vapply(seq_len(k), function(i) {
      v <- input[[paste0("cname_", k, "_", i)]]
      if (is.null(v) || !nzchar(v)) paste0("C", i) else v
    }, "")
    nm
  })
  excluded <- reactive({
    k <- input$k
    which(vapply(seq_len(k), function(i) isTRUE(input[[paste0("cexcl_", k, "_", i)]]), TRUE))
  })

  output$cluster_editor <- renderUI({
    req(rv$res); k <- input$k
    tagList(tags$label(class = "form-label", "Cluster names / exclude"),
      lapply(seq_len(k), function(i) div(class = "d-flex align-items-center gap-2 mb-1",
        span(style = sprintf("width:12px;height:12px;border-radius:2px;background:%s;flex:none", PAL[(i - 1) %% length(PAL) + 1])),
        div(style = "flex:1;min-width:0", class = "cl-name", textInput(paste0("cname_", k, "_", i), NULL, width = "100%",
                                        isolate(input[[paste0("cname_", k, "_", i)]]) %||% paste0("C", i))),
        div(class = "cl-excl", checkboxInput(paste0("cexcl_", k, "_", i), "excl.",
                                             isolate(input[[paste0("cexcl_", k, "_", i)]]) %||% FALSE)))))
  })

  umap_lims <- reactive({ req(rv$res); apply(rv$res$umap, 2, range) })

  draw_umap <- function(xy, col, main = NULL, cex = 0.5) {
    lim <- umap_lims()
    plot(xy, col = col, pch = 16, cex = cex, asp = 1, xlim = lim[, 1], ylim = lim[, 2],
         axes = FALSE, xlab = "UMAP1", ylab = "UMAP2", main = main, font.main = 1)
    box(col = "#e4e4e0")
  }

  output$umap_plot <- renderPlot({
    res <- rv$res; req(res)
    par(mar = c(3, 3, 2, 10), mgp = c(1.5, 0.5, 0), xpd = NA)
    by <- input$colour_by
    if (by == "cluster") {
      cc <- cl(); nm <- cluster_names()
      draw_umap(res$umap, PAL[(cc - 1) %% length(PAL) + 1], sprintf("k = %d", input$k))
      sz <- tabulate(cc, input$k)
      legend("topright", inset = c(-0.32, 0), bty = "n", pch = 16, col = PAL[seq_len(input$k)],
             legend = sprintf("%s  (%d)", nm, sz), cex = 0.9)
    } else if (by == "sample") {
      s <- factor(res$sample)
      draw_umap(res$umap, PAL[(as.integer(s) - 1) %% length(PAL) + 1], "by sample")
      legend("topright", inset = c(-0.32, 0), bty = "n", pch = 16,
             col = PAL[seq_along(levels(s))], legend = levels(s), cex = 0.9)
    } else {
      v <- res$expr[, by]
      q <- quantile(v, c(0.02, 0.98)); s <- pmin(pmax((v - q[1]) / diff(q), 0), 1)
      o <- order(s)
      ramp <- colorRampPalette(RAMP)(100)
      draw_umap(res$umap[o, ], ramp[1 + floor(99 * s[o])], paste(by, "(2nd-98th percentile)"))
      legend("topright", inset = c(-0.32, 0), bty = "n", pch = 15, col = ramp[c(1, 50, 100)],
             legend = c("low", "", "high"), cex = 0.9)
    }
  })

  output$facet_plot <- renderPlot({
    res <- rv$res; req(res)
    s <- unique(res$sample); n <- length(s)
    nc <- ceiling(sqrt(n)); nr <- ceiling(n / nc)
    par(mfrow = c(nr, nc), mar = c(1, 1, 2, 1))
    cc <- cl()
    for (x in s) {
      i <- res$sample == x
      draw_umap(res$umap, "#dcdcd8", sprintf("%s  (n = %d)", x, sum(i)), cex = 0.3)
      points(res$umap[i, , drop = FALSE], pch = 16, cex = 0.5, col = PAL[(cc[i] - 1) %% length(PAL) + 1])
    }
  })

  markers_k <- reactive({ req(rv$res); cluster_markers(rv$res$expr, cl()) })

  output$heatmap <- renderPlot({
    cm <- markers_k(); z <- pmin(pmax(cm$z, -2), 2)
    nm <- cluster_names()[as.integer(rownames(z))]
    par(mar = c(7, 9, 2, 2))
    image(t(z[nrow(z):1, , drop = FALSE]), col = DIV, zlim = c(-2, 2), axes = FALSE,
          main = "Cluster median vs all cells (SD units, clipped at ±2)", font.main = 1)
    axis(1, at = seq(0, 1, length.out = ncol(z)), labels = colnames(z), las = 2, tick = FALSE, cex.axis = 0.85)
    axis(2, at = seq(0, 1, length.out = nrow(z)), labels = rev(sprintf("%s (%d)", nm, tabulate(cl(), input$k))),
         las = 1, tick = FALSE, cex.axis = 0.9)
  })

  output$cluster_table <- renderTable({
    cm <- markers_k(); s <- cluster_summary(cm, cl())
    s$cluster <- cluster_names()[as.integer(s$cluster)]
    names(s) <- c("Cluster", "Cells", "%", "Highest AUC vs rest", "Lowest AUC vs rest")
    s
  }, striped = TRUE, spacing = "s")

  props <- reactive({
    req(rv$res)
    cc <- factor(cluster_names()[cl()], levels = cluster_names())
    cluster_proportions(rv$res$sample, cc, rv$res$groups,
                        exclude = cluster_names()[excluded()])
  })

  output$prop_plot <- renderPlot({
    p <- props()$summary
    p <- p[p$cluster %in% setdiff(cluster_names(), cluster_names()[excluded()]), ]
    p$cluster <- factor(p$cluster, levels = cluster_names())
    mm <- tapply(p$mean, list(p$group, p$cluster), sum)
    se <- tapply(p$sem, list(p$group, p$cluster), sum)
    mm <- mm[, colSums(!is.na(mm)) > 0, drop = FALSE]
    se <- se[, colnames(mm), drop = FALSE]
    gcol <- colorRampPalette(c("#2a78d6", "#e8743b", "#19a979", "#945ecf"))(max(2, nrow(mm)))[seq_len(nrow(mm))]
    par(mar = c(5, 4, 2, 1))
    bp <- barplot(mm, beside = TRUE, col = gcol, border = NA, las = 2,
                  ylim = c(0, max(mm + se, na.rm = TRUE) * 1.15), ylab = "% of cells (mean across samples)",
                  legend.text = if (nrow(mm) > 1) rownames(mm) else NULL,
                  args.legend = list(bty = "n", border = NA))
    ok <- !is.na(se)
    arrows(bp[ok], (mm - se)[ok], bp[ok], (mm + se)[ok], angle = 90, code = 3, length = 0.03, col = "#52514e")
  })

  output$prop_table <- renderTable(props()$per_sample, striped = TRUE, spacing = "s", digits = 1)

  # ---- downloads ------------------------------------------------------------------
  output$dl_cells <- downloadHandler("cells.csv", function(file) {
    res <- rv$res
    d <- data.frame(sample = res$sample, UMAP1 = res$umap[, 1], UMAP2 = res$umap[, 2],
                    som_node = res$som$node, cluster = cluster_names()[cl()], res$expr,
                    check.names = FALSE)
    write.csv(d, file, row.names = FALSE)
  })
  output$dl_markers <- downloadHandler("cluster_markers.csv", function(file) {
    cm <- markers_k(); nm <- cluster_names()[as.integer(rownames(cm$z))]
    d <- rbind(data.frame(cluster = nm, stat = "median", cm$median, check.names = FALSE),
               data.frame(cluster = nm, stat = "z_vs_all", cm$z, check.names = FALSE),
               data.frame(cluster = nm, stat = "auc_vs_rest", cm$auc, check.names = FALSE))
    write.csv(d, file, row.names = FALSE)
  })
  output$dl_props <- downloadHandler("cluster_proportions.csv", function(file) {
    p <- props()
    write.csv(p$per_sample, file, row.names = FALSE)
  })
  output$dl_params <- downloadHandler("run_parameters.txt", function(file) {
    p <- rv$res$params
    lines <- c(sprintf("date: %s", format(Sys.time())),
               sprintf("k: %d", input$k),
               sprintf("%s: %s", names(p), vapply(p, function(x) paste(x, collapse = "; "), "")),
               sprintf("cluster names: %s", paste(cluster_names(), collapse = "; ")),
               sprintf("excluded: %s", paste(cluster_names()[excluded()], collapse = "; ")),
               sprintf("samples: %s", paste(sprintf("%s=%s", names(rv$res$groups), rv$res$groups), collapse = "; ")))
    writeLines(lines, file)
  })
}

shinyApp(ui, server)
