#  ──────────────────────────────────────────────────────────────────────────
# PPI Network Analysis Pipeline
# Supports: STRING PPI download, degree-based network filtering, chord plot
# Input: candidate gene list
# Outputs:
#   1. STRING interaction table
#   2. cleaned PPI edge table
#   3. PPI chord plot PDF
#   4. PPI chord plot PNG
#  ──────────────────────────────────────────────────────────────────────────

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(tidyverse)
  library(circlize)
  library(ComplexHeatmap)
  library(Cairo)
  library(grid)
  library(curl)
  library(readr)
})


# ══════════════════════════════════════════════════════════════════════════════
# ██  USER CONFIGURATION — ONLY EDIT THIS BLOCK                             ██
# ══════════════════════════════════════════════════════════════════════════════
# 说明：
#   1. 标记为 [AI] 的参数表示可以根据不同项目进行修改。
#   2. 没有标记 [AI] 的函数主体和 main workflow 一般不建议修改。
#   3. 如果让 AI 修改脚本，优先只修改本配置区。
#   4. 本脚本用于 STRING PPI 网络分析，包括 PPI 数据下载、网络过滤和弦图绘制


# ── A. Project paths ───────────────────────────────────────────────────────
# [AI] 项目根目录
project_dir <- "D:/R/sci/04_ppi/01_data"                                     # [AI]

# [AI] PPI 结果输出目录
outdir     <- "D:/R/sci/04_ppi/03_result"                              # [AI]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# [AI] 候选基因文件
# 一般来自 Venn 交集、差异基因结果、WGCNA 模块基因或机器学习筛选结果
local_ppi_file <- file.path(outdir, "Target_PPI.csv")

# target_gene_file <- file.path(project_dir, "02_Venn", "Target_PPI.csv")      # [AI]


# ── B. Candidate gene column settings ──────────────────────────────────────
# [AI] 候选基因文件中可能出现的基因列名
# 如果输入文件的基因列名不在这里，可以继续添加
possible_gene_cols <- c(                                                  # [AI]
  "Gene", "gene", "GENE",
  "Symbol", "SYMBOL", "symbol",
  "gene_name", "Gene_Name", "gene_symbol",
  "GeneSymbol", "geneSymbol",
  "external_gene_name",
  "hgnc_symbol", "mgi_symbol"
)


# ── C. STRING database parameters ──────────────────────────────────────────
# [AI] STRING 物种编号
#   9606  = human / Homo sapiens
#   10090 = mouse / Mus musculus
#   10116 = rat / Rattus norvegicus
string_species_id <- 3635                                                # [AI]

# [AI] STRING 互作置信度阈值
#   150 = low confidence
#   400 = medium confidence
#   700 = high confidence
#   900 = highest confidence
string_min_score <- 400                                                  # [AI]


# ── D. PPI network filtering parameters ────────────────────────────────────
# [AI] 最小节点连接度
#   1 = 保留至少有 1 条互作边的基因，网络较完整
#   2 = 去除低连接节点，网络更简洁
#   3 或更高 = 网络更严格，但可能丢失较多基因
min_node_degree <- 1                                                     # [AI]


# ── E. Plot settings ───────────────────────────────────────────────────────
# [AI] PPI 网络图标题
plot_title <- "Protein-Protein Interaction Network"                       # [AI]

# [AI] 绘图字体
# Linux 服务器如果 Times 不可用，可改为 "serif"
plot_family <- "Times"                                                   # [AI]

# [AI] 图片宽度、高度和分辨率
plot_width  <- 9                                                         # [AI]
plot_height <- 6.6                                                       # [AI]
plot_dpi    <- 300                                                       # [AI]


# ── F. Output file names ───────────────────────────────────────────────────
# [AI] STRING 原始互作结果文件
string_output_file <- file.path(outdir, "string_interactions.tsv")        # [AI]

# [AI] 清洗和过滤后的 PPI 边文件
ppi_edge_file <- file.path(outdir, "ppi_edges_cleaned.csv")               # [AI]

# [AI] PPI 网络图 PDF 输出文件
ppi_pdf_file <- file.path(outdir, "ppi.pdf")                              # [AI]

# [AI] PPI 网络图 PNG 输出文件
ppi_png_file <- file.path(outdir, "ppi.png")                              # [AI]


# ══════════════════════════════════════════════════════════════════════════════
# ██  END OF USER CONFIGURATION                                             ██
# ══════════════════════════════════════════════════════════════════════════════
# [AI_KEEP] Do not modify the functions and main workflow below unless needed.



# ── 1. read_candidate_genes ────────────────────────────────────────────────
# Purpose:
#   Read candidate genes from a CSV file.
#   The function automatically detects common gene column names.
#   If no standard gene column is found, it re-reads the file without header
#   and uses the column with the largest number of non-empty values.

read_candidate_genes <- function(target_gene_file, possible_gene_cols) {
  
  if (!file.exists(target_gene_file)) {
    stop("Candidate gene file does not exist: ", target_gene_file)
  }
  
  df_header <- readr::read_csv(
    target_gene_file,
    show_col_types = FALSE
  )
  
  cat("Columns detected in candidate gene file:\n")
  print(colnames(df_header))
  
  matched_cols <- intersect(possible_gene_cols, colnames(df_header))
  
  if (length(matched_cols) > 0) {
    
    gene_col <- matched_cols[1]
    cat("Detected standard gene column:", gene_col, "\n")
    genes <- df_header[[gene_col]]
    
  } else {
    
    cat("No standard gene column detected. Re-reading without header.\n")
    
    df_no_header <- readr::read_csv(
      target_gene_file,
      col_names = FALSE,
      show_col_types = FALSE
    )
    
    non_empty_counts <- sapply(df_no_header, function(x) {
      x <- as.character(x)
      x <- trimws(x)
      sum(!is.na(x) & x != "")
    })
    
    print(non_empty_counts)
    
    gene_col <- names(which.max(non_empty_counts))
    cat("Using auto-selected gene column:", gene_col, "\n")
    genes <- df_no_header[[gene_col]]
  }
  
  genes <- genes %>%
    as.character() %>%
    trimws()
  
  genes <- genes[!is.na(genes)]
  genes <- genes[genes != ""]
  genes <- genes[!genes %in% c(possible_gene_cols, "X", "X1")]
  genes <- unique(genes)
  
  cat("Final candidate gene number:", length(genes), "\n")
  
  if (length(genes) == 0) {
    stop("No candidate genes detected from target gene file.")
  }
  
  return(genes)
}


# ── 2. download_string_ppi ─────────────────────────────────────────────────
# Purpose:
#   Query STRING database using gene symbols.
#   First, genes are mapped to STRING preferred names.
#   Then, the STRING interaction network is downloaded.

download_string_ppi <- function(
    gene_symbols,
    species_id,
    min_score,
    output_file
) {
  
  options(timeout = 600)
  
  string_api <- "https://string-db.org/api/tsv"
  
  gene_symbols <- unique(na.omit(gene_symbols))
  gene_symbols <- gene_symbols[gene_symbols != ""]
  
  if (length(gene_symbols) == 0) {
    stop("Input gene symbol vector is empty.")
  }
  
  gene_query <- paste(gene_symbols, collapse = "%0d")
  
  mapping_url <- paste0(
    string_api,
    "/get_string_ids?species=", species_id,
    "&identifiers=", gene_query
  )
  
  mapping_file <- tempfile(fileext = ".tsv")
  
  cat("Downloading STRING gene mapping table...\n")
  
  curl::curl_download(
    url = mapping_url,
    destfile = mapping_file
  )
  
  mapping_table <- read.delim(
    mapping_file,
    stringsAsFactors = FALSE
  )
  
  unlink(mapping_file)
  
  if (!"preferredName" %in% colnames(mapping_table)) {
    stop("STRING mapping result does not contain preferredName column.")
  }
  
  mapped_genes <- mapping_table %>%
    dplyr::pull(preferredName) %>%
    unique() %>%
    na.omit()
  
  cat("Mapped gene number in STRING:", length(mapped_genes), "\n")
  
  if (length(mapped_genes) == 0) {
    stop("No genes were successfully mapped by STRING.")
  }
  
  mapped_query <- paste(mapped_genes, collapse = "%0d")
  
  network_url <- paste0(
    string_api,
    "/network?species=", species_id,
    "&required_score=", min_score,
    "&identifiers=", mapped_query
  )
  
  cat("Downloading STRING PPI network...\n")
  
  curl::curl_download(
    url = network_url,
    destfile = output_file
  )
  
  cat("STRING PPI data downloaded:\n")
  cat(output_file, "\n")
  
  return(output_file)
}


# ── 3. read_string_network ─────────────────────────────────────────────────
# Purpose:
#   Read STRING interaction table and keep key columns:
#   preferredName_A, preferredName_B and score.

read_string_network <- function(string_output_file) {
  
  if (!file.exists(string_output_file)) {
    stop("STRING interaction file does not exist: ", string_output_file)
  }
  
  ppi_data_raw <- readr::read_tsv(
    string_output_file,
    show_col_types = FALSE
  )
  
  required_cols <- c("preferredName_A", "preferredName_B", "score")
  
  if (!all(required_cols %in% colnames(ppi_data_raw))) {
    stop(
      "STRING output does not contain required columns: ",
      paste(required_cols, collapse = ", ")
    )
  }
  
  ppi_edges <- ppi_data_raw %>%
    dplyr::select(
      protein_a = preferredName_A,
      protein_b = preferredName_B,
      score
    ) %>%
    dplyr::filter(
      !is.na(protein_a),
      !is.na(protein_b),
      protein_a != "",
      protein_b != "",
      protein_a != protein_b
    ) %>%
    dplyr::distinct()
  
  cat("Raw PPI edge number:", nrow(ppi_edges), "\n")
  cat(
    "Raw PPI node number:",
    length(unique(c(ppi_edges$protein_a, ppi_edges$protein_b))),
    "\n"
  )
  
  if (nrow(ppi_edges) == 0) {
    stop("No valid PPI edges found after reading STRING output.")
  }
  
  return(ppi_edges)
}


# ── 4. filter_ppi_by_degree ────────────────────────────────────────────────
# Purpose:
#   Iteratively remove nodes whose degree is smaller than min_degree.
#   This helps simplify the network and remove isolated or weakly connected nodes.

filter_ppi_by_degree <- function(ppi_edges, min_degree = 1) {
  
  if (nrow(ppi_edges) == 0) {
    stop("Input PPI edge table is empty.")
  }
  
  filtered_edges <- ppi_edges
  
  repeat {
    
    degree_table <- table(c(filtered_edges$protein_a, filtered_edges$protein_b))
    
    retained_nodes <- names(degree_table[degree_table >= min_degree])
    
    updated_edges <- filtered_edges %>%
      dplyr::filter(
        protein_a %in% retained_nodes,
        protein_b %in% retained_nodes
      )
    
    if (nrow(updated_edges) == nrow(filtered_edges)) {
      break
    }
    
    filtered_edges <- updated_edges
    
    if (nrow(filtered_edges) == 0) {
      stop("All PPI edges were removed by degree filtering. Please lower min_degree.")
    }
  }
  
  node_count <- length(unique(c(filtered_edges$protein_a, filtered_edges$protein_b)))
  edge_count <- nrow(filtered_edges)
  
  cat("Filtered PPI node number:", node_count, "\n")
  cat("Filtered PPI edge number:", edge_count, "\n")
  
  return(filtered_edges)
}


# ── 5. plot_ppi_chord ──────────────────────────────────────────────────────
# Purpose:
#   Draw a chord diagram for the PPI network.
#   Node color represents node degree.
#   Link color represents STRING interaction score.

plot_ppi_chord <- function(
    ppi_edges,
    min_degree,
    main_title,
    font_family
) {
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  
  par(family = font_family)
  
  plot_edges <- filter_ppi_by_degree(
    ppi_edges = ppi_edges,
    min_degree = min_degree
  )
  
  chord_data <- plot_edges %>%
    dplyr::transmute(
      from = protein_a,
      to = protein_b,
      value = score
    )
  
  all_nodes <- unique(c(chord_data$from, chord_data$to))
  
  node_degree <- table(c(chord_data$from, chord_data$to))
  node_degree <- sort(node_degree[all_nodes], decreasing = TRUE)
  node_order <- names(node_degree)
  
  degree_range <- range(as.numeric(node_degree), na.rm = TRUE)
  score_range <- range(chord_data$value, na.rm = TRUE)
  
  if (degree_range[1] == degree_range[2]) {
    degree_range <- c(degree_range[1], degree_range[1] + 1)
  }
  
  if (score_range[1] == score_range[2]) {
    score_range <- c(score_range[1], score_range[1] + 0.001)
  }
  
  degree_color_fun <- circlize::colorRamp2(
    c(degree_range[1], mean(degree_range), degree_range[2]),
    c("#4575B4", "#FFFFBF", "#D73027")
  )
  
  score_color_fun <- circlize::colorRamp2(
    c(score_range[1], mean(score_range), score_range[2]),
    c("#A6CEE3", "#984EA3", "#B2182B")
  )
  
  node_colors <- stats::setNames(
    degree_color_fun(as.numeric(node_degree)),
    names(node_degree)
  )
  
  circos.clear()
  
  circos.par(
    start.degree = 90,
    gap.degree = max(0.5, min(2.5, 360 / length(all_nodes) / 3)),
    canvas.xlim = c(-1.05, 1.05),
    canvas.ylim = c(-1.15, 1.15)
  )
  
  chordDiagram(
    x = chord_data,
    order = node_order,
    grid.col = node_colors,
    col = score_color_fun(chord_data$value),
    transparency = 0.45,
    link.lwd = 1.3,
    link.sort = TRUE,
    link.decreasing = TRUE,
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.06)
  )
  
  circos.trackPlotRegion(
    track.index = 1,
    bg.border = NA,
    panel.fun = function(x, y) {
      
      gene_name <- get.cell.meta.data("sector.index")
      x_center <- mean(get.cell.meta.data("xlim"))
      y_pos <- get.cell.meta.data("ylim")[1]
      
      circos.text(
        x = x_center,
        y = y_pos,
        labels = gene_name,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        cex = 0.65,
        family = font_family
      )
    }
  )
  
  title(
    main = main_title,
    line = 0.0,
    cex.main = 1.15,
    font.main = 2,
    family = font_family
  )
  
  degree_legend <- ComplexHeatmap::Legend(
    title = "Node degree",
    col_fun = degree_color_fun,
    at = unique(round(seq(degree_range[1], degree_range[2], length.out = 5))),
    direction = "horizontal",
    title_position = "topcenter",
    title_gp = grid::gpar(
      fontsize = 10,
      fontface = "bold",
      fontfamily = font_family
    ),
    labels_gp = grid::gpar(
      fontsize = 8,
      fontfamily = font_family
    )
  )
  
  score_legend <- ComplexHeatmap::Legend(
    title = "STRING score",
    col_fun = score_color_fun,
    at = round(seq(score_range[1], score_range[2], length.out = 5), 3),
    direction = "horizontal",
    title_position = "topcenter",
    title_gp = grid::gpar(
      fontsize = 10,
      fontface = "bold",
      fontfamily = font_family
    ),
    labels_gp = grid::gpar(
      fontsize = 8,
      fontfamily = font_family
    )
  )
  
  legend_group <- ComplexHeatmap::packLegend(
    degree_legend,
    score_legend,
    direction = "vertical",
    gap = grid::unit(4, "mm")
  )
  
  grid::pushViewport(
    grid::viewport(
      x = 0.84,
      y = 0.15,
      width = 0.32,
      height = 0.28
    )
  )
  
  ComplexHeatmap::draw(legend_group)
  grid::popViewport()
  
  circos.clear()
}


# ── 6. save_ppi_plots ──────────────────────────────────────────────────────
# Purpose:
#   Save the PPI chord plot as PDF and PNG.

save_ppi_plots <- function(
    ppi_edges,
    pdf_file,
    png_file,
    min_degree,
    main_title,
    font_family,
    plot_width,
    plot_height,
    plot_dpi
) {
  
  Cairo::CairoPDF(
    file = pdf_file,
    width = plot_width,
    height = plot_height,
    family = font_family
  )
  
  plot_ppi_chord(
    ppi_edges = ppi_edges,
    min_degree = min_degree,
    main_title = main_title,
    font_family = font_family
  )
  
  dev.off()
  
  Cairo::CairoPNG(
    filename = png_file,
    width = plot_width,
    height = plot_height,
    units = "in",
    res = plot_dpi
  )
  
  plot_ppi_chord(
    ppi_edges = ppi_edges,
    min_degree = min_degree,
    main_title = main_title,
    font_family = font_family
  )
  
  dev.off()
  
  cat("PPI plots saved:\n")
  cat(pdf_file, "\n")
  cat(png_file, "\n")
}

# ── 7. main workflow ───────────────────────────────────────────────────────

setwd(project_dir)

cat("\nPPI network analysis started.\n")
cat("Results will be saved to:\n")
cat(outdir, "\n\n")

# ------------------------------------------------------------
# 不再从单独的基因文件读取候选基因，因为你的 PPI 文件已经包含所有节点
# 如果你仍想保存一个节点列表，可以从 PPI 边表中提取
# ------------------------------------------------------------

# 读取本地 PPI 文件（跳过 STRING 下载）
cat("Reading local PPI file:", local_ppi_file, "\n")
ppi_data_local <- readr::read_csv(local_ppi_file, show_col_types = FALSE)

# 检查必要的列
required_local_cols <- c("from_symbol", "to_symbol", "combined_score")
if (!all(required_local_cols %in% colnames(ppi_data_local))) {
  stop("Your PPI file must contain columns: from_symbol, to_symbol, combined_score")
}

# 转换为脚本内部标准格式：protein_a, protein_b, score
ppi_edges_raw <- ppi_data_local %>%
  dplyr::select(
    protein_a = from_symbol,
    protein_b = to_symbol,
    score = combined_score
  ) %>%
  dplyr::filter(
    !is.na(protein_a),
    !is.na(protein_b),
    protein_a != "",
    protein_b != "",
    protein_a != protein_b
  ) %>%
  dplyr::distinct()

cat("Raw PPI edge number:", nrow(ppi_edges_raw), "\n")
cat("Raw PPI node number:",
    length(unique(c(ppi_edges_raw$protein_a, ppi_edges_raw$protein_b))),
    "\n")

if (nrow(ppi_edges_raw) == 0) {
  stop("No valid PPI edges found in your local file.")
}

# 可选：将原始边表保存一份（相当于原脚本中的 string_interactions.tsv）
readr::write_csv(ppi_edges_raw, file.path(outdir, "raw_ppi_edges.csv"))

# 使用 degree 过滤（与原来一样）
ppi_edges_filtered <- filter_ppi_by_degree(
  ppi_edges = ppi_edges_raw,
  min_degree = min_node_degree
)

# 保存过滤后的边表
readr::write_csv(ppi_edges_filtered, ppi_edge_file)

# 绘图（与原来一样）
save_ppi_plots(
  ppi_edges = ppi_edges_filtered,
  pdf_file = ppi_pdf_file,
  png_file = ppi_png_file,
  min_degree = min_node_degree,
  main_title = plot_title,
  font_family = plot_family,
  plot_width = plot_width,
  plot_height = plot_height,
  plot_dpi = plot_dpi
)

cat("\nPPI network analysis finished.\n")
cat("Results saved to:\n")
cat(outdir, "\n")
