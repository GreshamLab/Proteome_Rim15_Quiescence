# Rim15 Proteomics Data Browser
# Companion app for Sun et al. (2026), PLOS Genetics
#
# Run from repo root:  shiny::runApp("shiny/")
# First-time setup:    Rscript shiny/prep_data.R

library(shiny)
library(bslib)
library(tidyverse)
library(plotly)
library(DT)
library(shinycssloaders)

# ── Load preprocessed data ────────────────────────────────────────────────────
# Paths are relative to shiny/ (Shiny sets wd to the app directory at startup)
protein_long  <- readRDS("data/protein_long.rds")
protein_de    <- readRDS("data/protein_de.rds")
psite_long    <- readRDS("data/psite_long.rds")
psite_de      <- readRDS("data/psite_de.rds")
gene_universe <- readRDS("data/gene_universe.rds")

ALL_GENES  <- gene_universe$genes
GO_SETS    <- gene_universe$go_gene_sets

TIME_ORDER   <- c(0L, 6L, 16L, 30L)
NUTRIENTS    <- c("C" = "Carbon (−C)", "P" = "Phosphorus (−P)")
GENO_COLORS  <- c("WT" = "#E69F00", "rim15δ" = "#56B4E9")
SIG_COLORS   <- c("Up in WT" = "#E69F00", "Down in WT" = "#56B4E9", "n.s." = "#CCCCCC")

# ── Helper functions ──────────────────────────────────────────────────────────

#' Build a plotly time-series for protein abundance.
#' d: filtered slice of protein_long (columns: gene, genotype, nutrient, time_h, ratio)
#' show_reps: overlay individual replicate points if TRUE
make_timeseries <- function(d, show_reps = FALSE) {
  summ <- d %>%
    group_by(gene, genotype, nutrient, time_h) %>%
    summarise(
      mean_l2r = mean(log2(ratio), na.rm = TRUE),
      sem_l2r  = sd(log2(ratio),   na.rm = TRUE) / sqrt(sum(!is.na(ratio))),
      .groups  = "drop"
    ) %>%
    mutate(
      nut_label  = NUTRIENTS[nutrient],
      ymin       = mean_l2r - sem_l2r,
      ymax       = mean_l2r + sem_l2r,
      hover_text = sprintf("%s | %s | %s | %d h<br>%.3f ± %.3f",
                           gene, genotype, nut_label, time_h, mean_l2r, sem_l2r)
    )

  p <- ggplot(summ,
              aes(x = time_h, y = mean_l2r, color = genotype, fill = genotype,
                  group = interaction(gene, genotype, nutrient),
                  text = hover_text)) +
    geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2, color = NA) +
    geom_line() +
    geom_point(size = 2.5) +
    scale_color_manual(values = GENO_COLORS) +
    scale_fill_manual(values  = GENO_COLORS) +
    scale_x_continuous(breaks = TIME_ORDER) +
    facet_grid(nut_label ~ gene) +
    labs(x = "Time (h)", y = "log₂ SILAC ratio",
         color = "Genotype", fill = "Genotype") +
    theme_bw(base_size = 13) +
    theme(strip.background = element_rect(fill = "#f2f2f2"),
          legend.position  = "bottom")

  if (show_reps) {
    reps <- d %>% mutate(nut_label = NUTRIENTS[nutrient])
    p <- p +
      geom_point(data = reps,
                 aes(x = time_h, y = log2(ratio), color = genotype),
                 size = 1.5, alpha = 0.5, shape = 1, inherit.aes = FALSE)
  }

  ggplotly(p, tooltip = "text") %>%
    layout(legend = list(orientation = "h", y = -0.18))
}

#' Build a plotly time-series for phosphosite abundance.
#' d: filtered slice of psite_long
make_psite_timeseries <- function(d) {
  summ <- d %>%
    group_by(psite, gene, genotype, nutrient, time_h) %>%
    summarise(
      mean_l2r = mean(log2(ratio), na.rm = TRUE),
      sem_l2r  = sd(log2(ratio),   na.rm = TRUE) / sqrt(sum(!is.na(ratio))),
      .groups  = "drop"
    ) %>%
    mutate(
      nut_label  = NUTRIENTS[nutrient],
      hover_text = sprintf("%s | %s | %s | %d h<br>%.3f ± %.3f",
                           psite, genotype, nut_label, time_h, mean_l2r, sem_l2r)
    )

  p <- ggplot(summ,
              aes(x = time_h, y = mean_l2r, color = genotype,
                  linetype = psite,
                  group = interaction(psite, genotype, nutrient),
                  text = hover_text)) +
    geom_line() +
    geom_point(size = 2) +
    scale_color_manual(values = GENO_COLORS) +
    scale_x_continuous(breaks = TIME_ORDER) +
    facet_grid(nut_label ~ gene) +
    labs(x = "Time (h)", y = "log₂ phosphosite ratio",
         color = "Genotype", linetype = "Phosphosite") +
    theme_bw(base_size = 13) +
    theme(strip.background = element_rect(fill = "#f2f2f2"),
          legend.position  = "bottom")

  ggplotly(p, tooltip = "text") %>%
    layout(legend = list(orientation = "h", y = -0.18))
}

# ── UI ─────────────────────────────────────────────────────────────────────────
ui <- page_navbar(
  title = "Rim15 Proteomics Browser",
  theme = bs_theme(bootswatch  = "flatly",
                   primary      = "#2C6FAC",
                   base_font    = font_google("Inter"),
                   heading_font = font_google("Inter")),

  # ── Tab 1: Gene Browser ───────────────────────────────────────────────────
  nav_panel(
    title = "Gene Browser",
    layout_sidebar(
      sidebar = sidebar(
        width = 290,
        selectizeInput(
          "genes", "Gene(s)",
          choices = NULL,    # populated server-side for speed
          multiple = TRUE,
          options  = list(placeholder = "e.g. IGO2, SIC1, ATG1")
        ),
        radioButtons(
          "data_type", "Data to display",
          choices  = c("Proteome" = "prot",
                       "Phosphoproteome" = "psite",
                       "Both" = "both"),
          selected = "prot"
        ),
        radioButtons(
          "nutrient_filter", "Starvation condition",
          choices  = c("All" = "all", "Carbon" = "C", "Phosphorus" = "P"),
          selected = "all"
        ),
        checkboxInput("show_reps", "Show individual replicates", value = FALSE),
        hr(),
        downloadButton("dl_prot", "Download protein data",  class = "btn-sm w-100"),
        downloadButton("dl_psite", "Download pSite data",   class = "btn-sm w-100 mt-1")
      ),
      uiOutput("browser_main")
    )
  ),

  # ── Tab 2: Volcano ────────────────────────────────────────────────────────
  nav_panel(
    title = "Volcano Plot",
    layout_sidebar(
      sidebar = sidebar(
        width = 270,
        radioButtons(
          "vol_type", "Data type",
          choices  = c("Proteome" = "prot", "Phosphoproteome" = "psite"),
          selected = "prot"
        ),
        selectInput(
          "vol_nutrient", "Condition",
          choices = c("Carbon (−C)" = "C", "Phosphorus (−P)" = "P")
        ),
        selectInput(
          "vol_time", "Timepoint",
          choices  = c("6 h" = "T06", "16 h" = "T16", "30 h" = "T30"),
          selected = "T30"
        ),
        sliderInput(
          "vol_q", "q-value threshold",
          min = 0.001, max = 0.2, value = 0.05, step = 0.001
        ),
        helpText(icon("info-circle"),
                 "Click a point to add the gene to the Gene Browser.")
      ),
      withSpinner(plotlyOutput("volcano_plot", height = "580px")),
      textOutput("volcano_summary")
    )
  ),

  # ── Tab 3: Heatmap ────────────────────────────────────────────────────────
  nav_panel(
    title = "Heatmap",
    layout_sidebar(
      sidebar = sidebar(
        width = 310,
        h6("Enter genes directly:", class = "mb-1"),
        selectizeInput(
          "hm_genes", "Gene(s)",
          choices  = NULL,
          multiple = TRUE,
          options  = list(placeholder = "Type gene names...")
        ),
        hr(),
        h6("...or expand a GO term (BP):", class = "mb-1"),
        selectizeInput(
          "go_term", "GO term",
          choices = NULL,
          options = list(placeholder = "Search GO terms...")
        ),
        actionButton("go_expand", "Add genes from GO term",
                     class = "btn-sm btn-secondary mb-3"),
        hr(),
        checkboxGroupInput(
          "hm_nutrients", "Conditions",
          choices  = c("Carbon (−C)" = "C", "Phosphorus (−P)" = "P"),
          selected = c("C", "P")
        ),
        checkboxGroupInput(
          "hm_genotypes", "Genotypes",
          choices  = c("WT", "rim15δ"),
          selected = c("WT", "rim15δ")
        ),
        checkboxInput("hm_cluster_rows", "Cluster genes (rows)",     value = TRUE),
        checkboxInput("hm_cluster_cols", "Cluster conditions (cols)", value = FALSE)
      ),
      withSpinner(plotlyOutput("heatmap_plot", height = "680px"))
    )
  ),

  # ── Tab 4: About ─────────────────────────────────────────────────────────
  nav_panel(
    title = "About",
    div(class = "container mt-4", style = "max-width:760px",
      h3("About this app"),
      p("Interactive data companion to:"),
      tags$blockquote(class = "blockquote",
        "Sun Y", tags$em(" et al."),
        " Rim15 orchestrates mitochondrial and cytoplasmic remodeling in",
        " quiescent ", tags$em("Saccharomyces cerevisiae"), " cells.",
        br(),
        tags$em("PLOS Genetics"), " (2026). PGENETICS-D-25-01150."
      ),
      h4("Experimental design"),
      tags$ul(
        tags$li("Quantitative SILAC proteomics and phosphoproteomics"),
        tags$li("Genotypes: wild-type (WT) and ", tags$em("rim15δ")),
        tags$li("Starvation conditions: carbon (", tags$em("−C"),
                ") and phosphorus (", tags$em("−P"), ")"),
        tags$li("Time-course: 0 / 6 / 16 / 30 h post-starvation, 3 biological replicates"),
        tags$li("y-axis shows log₂ SILAC ratio (protein or phosphosite abundance",
                " relative to the mixed-channel reference)"),
        tags$li("Positive values indicate enrichment in the heavy label (quiescent) channel")
      ),
      h4("Data availability"),
      tags$ul(
        tags$li("Raw mass spectrometry data: ProteomeXchange ",
                tags$a("PXD044239",
                       href   = "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD044239",
                       target = "_blank")),
        tags$li("Analysis code: ",
                tags$a("GitHub — GreshamLab/Proteome_Rim15_Quiescence",
                       href   = "https://github.com/GreshamLab/Proteome_Rim15_Quiescence",
                       target = "_blank"))
      ),
      h4("Contact"),
      p("David Gresham — ",
        tags$a("dgresham@nyu.edu", href = "mailto:dgresham@nyu.edu"),
        " — ",
        tags$a("Gresham Lab, NYU Biology",
               href = "https://greshamlab.bio.nyu.edu", target = "_blank"))
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # Populate server-side selectize inputs (avoids sending 2000+ choices to browser)
  updateSelectizeInput(session, "genes",    choices = ALL_GENES, server = TRUE)
  updateSelectizeInput(session, "hm_genes", choices = ALL_GENES, server = TRUE)
  updateSelectizeInput(session, "go_term",  choices = GO_SETS$TERM, server = TRUE)

  # ── Tab 1: Gene Browser ──────────────────────────────────────────────────

  # Reactive gene-filtered slices
  prot_data <- reactive({
    req(input$genes)
    d <- protein_long %>% filter(gene %in% input$genes)
    if (input$nutrient_filter != "all") d <- filter(d, nutrient == input$nutrient_filter)
    d
  })

  psite_data <- reactive({
    req(input$genes)
    d <- psite_long %>% filter(gene %in% input$genes)
    if (input$nutrient_filter != "all") d <- filter(d, nutrient == input$nutrient_filter)
    d
  })

  # Dynamic main panel: placeholder when no gene is selected
  output$browser_main <- renderUI({
    if (is.null(input$genes) || length(input$genes) == 0) {
      return(
        div(class = "d-flex justify-content-center align-items-center",
            style = "height:400px; color:#888;",
            div(class = "text-center",
                tags$i(class = "bi bi-search", style = "font-size:2rem"),
                p("Enter one or more gene names to explore their abundance profiles.",
                  style = "margin-top:1rem; font-size:1.1rem"),
                p(tags$em("Examples: IGO2, SIC1, ATG1, GDH1, TOR1"),
                  style = "font-size:0.9rem")))
      )
    }

    plots <- tagList()
    if (input$data_type %in% c("prot", "both")) {
      plots <- tagList(plots,
        h5("Protein abundance"),
        withSpinner(plotlyOutput("ts_prot", height = "380px")),
        br(),
        withSpinner(DTOutput("tbl_prot_de"))
      )
    }
    if (input$data_type %in% c("psite", "both")) {
      plots <- tagList(plots,
        h5("Phosphosite abundance", class = "mt-4"),
        withSpinner(plotlyOutput("ts_psite", height = "350px")),
        br(),
        withSpinner(DTOutput("tbl_psite_de"))
      )
    }
    plots
  })

  output$ts_prot <- renderPlotly({
    d <- prot_data()
    validate(need(nrow(d) > 0, "No protein-level data for the selected gene(s)."))
    make_timeseries(d, show_reps = isTRUE(input$show_reps))
  })

  output$ts_psite <- renderPlotly({
    d <- psite_data()
    validate(need(nrow(d) > 0,
                  "No phosphosite data for selected gene(s) in this condition."))
    make_psite_timeseries(d)
  })

  output$tbl_prot_de <- renderDT({
    req(input$genes)
    protein_de %>%
      filter(gene %in% input$genes) %>%
      dplyr::rename(Gene = gene, Condition = nutrient, Timepoint = time,
                    `p-value` = pvalue, `q-value` = qvalue, `log2FC` = log2FC) %>%
      mutate(across(where(is.numeric), ~ signif(.x, 3))) %>%
      datatable(
        caption  = "Differential protein abundance (WT vs rim15δ)",
        rownames = FALSE,
        options  = list(pageLength = 8, scrollX = TRUE, dom = "tip")
      )
  }, server = FALSE)

  output$tbl_psite_de <- renderDT({
    req(input$genes)
    psite_de %>%
      filter(gene %in% input$genes) %>%
      dplyr::rename(Gene = gene, Phosphosite = psite, Condition = nutrient,
                    Timepoint = time, `p-value` = pvalue, `q-value` = qvalue,
                    `log2FC` = log2FC) %>%
      mutate(across(where(is.numeric), ~ signif(.x, 3))) %>%
      datatable(
        caption  = "Differential phosphorylation (WT vs rim15δ)",
        rownames = FALSE,
        options  = list(pageLength = 8, scrollX = TRUE, dom = "tip")
      )
  }, server = FALSE)

  # Volcano click → add gene to browser tab
  observeEvent(event_data("plotly_click", source = "volcano"), {
    click <- event_data("plotly_click", source = "volcano")
    gene_clicked <- click$customdata
    if (!is.null(gene_clicked) && gene_clicked %in% ALL_GENES) {
      current <- isolate(input$genes)
      updateSelectizeInput(session, "genes",
                           choices  = ALL_GENES,
                           selected = unique(c(current, gene_clicked)),
                           server   = TRUE)
    }
  })

  output$dl_prot <- downloadHandler(
    filename = function() {
      paste0("rim15_protein_", paste(input$genes, collapse = "-"), ".csv")
    },
    content = function(file) write_csv(prot_data(), file)
  )

  output$dl_psite <- downloadHandler(
    filename = function() {
      paste0("rim15_psite_", paste(input$genes, collapse = "-"), ".csv")
    },
    content = function(file) write_csv(psite_data(), file)
  )

  # ── Tab 2: Volcano ────────────────────────────────────────────────────────

  vol_data <- reactive({
    if (input$vol_type == "prot") {
      protein_de %>%
        filter(nutrient == input$vol_nutrient,
               time     == input$vol_time) %>%
        mutate(
          label      = gene,
          neg_log10q = -log10(pmax(qvalue, 1e-10)),
          sig        = case_when(
            qvalue < input$vol_q & log2FC > 0 ~ "Up in WT",
            qvalue < input$vol_q & log2FC < 0 ~ "Down in WT",
            TRUE ~ "n.s."
          )
        )
    } else {
      psite_de %>%
        filter(nutrient == input$vol_nutrient,
               time     == input$vol_time) %>%
        mutate(
          label      = psite,
          neg_log10q = -log10(pmax(qvalue, 1e-10)),
          sig        = case_when(
            qvalue < input$vol_q & log2FC > 0 ~ "Up in WT",
            qvalue < input$vol_q & log2FC < 0 ~ "Down in WT",
            TRUE ~ "n.s."
          )
        )
    }
  })

  output$volcano_plot <- renderPlotly({
    d <- vol_data()
    validate(need(nrow(d) > 0, "No data for this condition and timepoint."))

    max_y  <- max(d$neg_log10q, na.rm = TRUE) * 1.05
    q_line <- -log10(input$vol_q)

    plot_ly(
      data          = d,
      x             = ~log2FC,
      y             = ~neg_log10q,
      color         = ~sig,
      colors        = SIG_COLORS,
      type          = "scatter",
      mode          = "markers",
      customdata    = ~gene,    # gene name passed to click handler
      text          = ~label,
      hovertemplate = paste0(
        "<b>%{text}</b><br>",
        "log₂FC: %{x:.2f}<br>",
        "-log₁₀(q): %{y:.2f}",
        "<extra></extra>"
      ),
      marker = list(size = 7, opacity = 0.75),
      source = "volcano"        # named source for event_data()
    ) %>%
      # Vertical reference line at log2FC = 0
      add_segments(x = 0, xend = 0, y = 0, yend = max_y,
                   line = list(dash = "dot", color = "grey60", width = 1),
                   showlegend = FALSE, inherit = FALSE) %>%
      # Horizontal significance threshold
      add_segments(x = min(d$log2FC, na.rm = TRUE),
                   xend = max(d$log2FC, na.rm = TRUE),
                   y = q_line, yend = q_line,
                   line = list(dash = "dot", color = "grey60", width = 1),
                   showlegend = FALSE, inherit = FALSE) %>%
      layout(
        xaxis  = list(title = "log₂ fold-change (WT / rim15δ)"),
        yaxis  = list(title = "-log₁₀(q-value)", range = c(0, max_y)),
        legend = list(orientation = "h", y = -0.12)
      )
  })

  output$volcano_summary <- renderText({
    d <- vol_data()
    sprintf(
      "%d proteins/sites up in WT, %d down in WT (q < %.3f, %s %s)",
      sum(d$sig == "Up in WT",   na.rm = TRUE),
      sum(d$sig == "Down in WT", na.rm = TRUE),
      input$vol_q,
      ifelse(input$vol_nutrient == "C", "Carbon", "Phosphorus"),
      input$vol_time
    )
  })

  # ── Tab 3: Heatmap ────────────────────────────────────────────────────────

  # Expand GO term into gene list and add to heatmap gene selectize
  observeEvent(input$go_expand, {
    req(input$go_term)
    go_genes <- GO_SETS %>%
      filter(TERM == input$go_term) %>%
      pull(genes)
    if (length(go_genes) == 0) return()
    go_genes <- intersect(go_genes[[1]], ALL_GENES)
    if (length(go_genes) > 200) {
      showNotification(
        sprintf("GO term contains %d genes; showing first 200.", length(go_genes)),
        type = "warning", duration = 6
      )
      go_genes <- go_genes[seq_len(200)]
    }
    current <- isolate(input$hm_genes)
    updateSelectizeInput(session, "hm_genes",
                         choices  = ALL_GENES,
                         selected = unique(c(current, go_genes)),
                         server   = TRUE)
  })

  output$heatmap_plot <- renderPlotly({
    req(input$hm_genes, length(input$hm_genes) >= 2)
    req(length(input$hm_nutrients) >= 1, length(input$hm_genotypes) >= 1)

    # Mean log2 ratio per gene × genotype × nutrient × timepoint
    hm <- protein_long %>%
      filter(gene     %in% input$hm_genes,
             nutrient %in% input$hm_nutrients,
             genotype %in% input$hm_genotypes) %>%
      group_by(gene, genotype, nutrient, time_h) %>%
      summarise(mean_l2r = mean(log2(ratio), na.rm = TRUE), .groups = "drop") %>%
      mutate(
        col_key = sprintf("%s | %s | %dh", genotype, NUTRIENTS[nutrient], time_h)
      )

    validate(need(nrow(hm) > 0, "No data for the selected genes and conditions."))

    # Deterministic column order: nutrient, genotype, time
    col_order <- hm %>%
      distinct(col_key, nutrient, genotype, time_h) %>%
      arrange(nutrient, genotype, time_h) %>%
      pull(col_key)

    # Pivot to matrix
    mat <- hm %>%
      pivot_wider(id_cols = gene, names_from = col_key, values_from = mean_l2r) %>%
      tibble::column_to_rownames("gene")

    col_order <- intersect(col_order, names(mat))
    mat <- as.matrix(mat[, col_order, drop = FALSE])

    # Row z-score (scale each gene independently)
    row_m   <- rowMeans(mat, na.rm = TRUE)
    row_sd  <- apply(mat, 1, sd, na.rm = TRUE)
    mat_z   <- sweep(sweep(mat, 1, row_m, "-"), 1, pmax(row_sd, 1e-6), "/")
    mat_z[is.na(mat_z)] <- 0   # impute remaining NAs as row mean (= 0 after z-score)

    # Optional hierarchical clustering
    row_idx <- if (input$hm_cluster_rows && nrow(mat_z) >= 2) {
      hclust(dist(mat_z))$order
    } else seq_len(nrow(mat_z))

    col_idx <- if (input$hm_cluster_cols && ncol(mat_z) >= 2) {
      hclust(dist(t(mat_z)))$order
    } else seq_len(ncol(mat_z))

    mat_out <- mat_z[row_idx, col_idx, drop = FALSE]

    plot_ly(
      z             = mat_out,
      x             = colnames(mat_out),
      y             = rownames(mat_out),
      type          = "heatmap",
      colorscale    = "RdBu",
      reversescale  = TRUE,
      zmid          = 0,
      hovertemplate = "Gene: %{y}<br>Sample: %{x}<br>z-score: %{z:.2f}<extra></extra>"
    ) %>%
      layout(
        xaxis  = list(tickangle = -40, title = "", tickfont = list(size = 10)),
        yaxis  = list(title = "", autorange = "reversed",
                      tickfont = list(size = 10)),
        margin = list(l = 90, b = 130, t = 20)
      )
  })
}

shinyApp(ui, server)
