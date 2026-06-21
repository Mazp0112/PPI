# PPI
可以生成一份精美的ppi图
# High-Quality PPI Network Visualization Code
自研R脚本用于蛋白互作网络(PPI)构建与美化绘图，出图效果优于STRING、Cytoscape默认可视化方案，可直接产出期刊高清配图。
也可以直接联系我
## 项目亮点
1. 支持STRING互作数据导入，根据STRING置信分数调整连线粗细
2. 节点Degree值梯度配色，直观区分核心Hub基因
3. 自定义布局、线条、图例，无杂乱锯齿，排版精致美观
4. 一键导出PDF/PNG/TIFF高分辨率图片，满足期刊投稿要求

## 文件说明
- 全部.R文件：完整PPI分析+可视化绘图代码
- ppi.pdf：测试数据集生成的成品精美PPI网络图
  本测试集选取的物种是陆地棉，当然不同的品种用不同的string来生成不同的输入文件
## 效果图展示
[ppi.pdf](https://github.com/user-attachments/files/29175601/ppi.pdf)
## 使用说明
1. 准备基因列表与STRING互作表格
2. 修改代码内文件路径、配色、尺寸参数
3. 直接运行脚本即可生成美化版PPI网络
# 下面是对代码的描述
---
name: ppi-network-analysis
description: Generate, modify, explain, and standardize R scripts for STRING PPI network analysis, including candidate gene input, STRING download, score/degree filtering, top hub genes, chord plots, and AI-editable configuration blocks.
---

# PPI 网络分析 Skill 简版

## 使用场景

用户需要写、修改或解释 STRING PPI 网络分析 R 脚本时使用本 reference，尤其是：

- 候选基因做 STRING PPI。
- 调整物种、STRING score、degree、top hub genes。
- 图太乱，只展示核心基因。
- 绘制 PPI 弦图并输出 PDF/PNG。
- 把普通 PPI 脚本整理成 AI 易修改格式。

## 修改原则

- 修改已有脚本时，保留原始分析逻辑。
- 如果有 `USER CONFIGURATION - ONLY EDIT THIS BLOCK`，优先只改配置区。
- `# [AI]` 表示可根据项目修改。
- 函数主体默认不改，除非用户要求、代码报错、输入格式变化或需要新增 top hub 功能。
- 物种默认先用人类 `string_species_id <- 9606`；如果 STRING 映射匹配率明显偏低，再结合基因名特征、上下文说明或试跑 mouse/rat 的匹配率辅助判断。
- 如果依据不足、多个物种匹配率接近，先询问用户。
- 不要擅自改 STRING 阈值等分析参数；如需调整，应说明原因。

## 参考代码

完整 PPI 脚本模板在：

```text
skills/rnaseq-workflow-agent/scripts/04_PPI.R
```

需要生成完整代码、查看函数实现或对齐脚本结构时，优先读取该文件作为参考；本 reference 只保留决策规则，不重复粘贴完整代码。

## 新脚本结构

最小结构：

```text
packages -> user config -> read genes -> download STRING -> clean edges ->
filter by degree/top hub -> chord plot -> save PDF/PNG -> main workflow
```

推荐包：`tidyverse`, `circlize`, `ComplexHeatmap`, `Cairo`, `grid`, `curl`, `readr`。

如果 `Cairo` 报错，用 `pdf()` / `png()`；如果 `Times` 字体不可用，用 `serif`。

## 配置区

生成或整理脚本时，配置区保留这些变量即可：

```r
project_dir <- "/path/to/project"                                      # [AI]
outdir <- file.path(project_dir, "03_PPI")                             # [AI]
target_gene_file <- file.path(project_dir, "02_Venn", "overlap.csv")   # [AI]

possible_gene_cols <- c(
  "Gene", "gene", "GENE", "Symbol", "SYMBOL", "symbol",
  "gene_name", "gene_symbol", "GeneSymbol",
  "external_gene_name", "hgnc_symbol", "mgi_symbol"
)                                                                      # [AI]

string_species_id <- 9606                                              # [AI]
string_min_score <- 400                                                # [AI]
min_node_degree <- 1                                                   # [AI]
top_degree_n <- NA                                                     # [AI]

plot_title <- "Protein-Protein Interaction Network"                    # [AI]
plot_family <- "Times"                                                 # [AI]
plot_width <- 9; plot_height <- 6.6; plot_dpi <- 300                   # [AI]

string_output_file <- file.path(outdir, "string_interactions.tsv")      # [AI]
ppi_edge_file <- file.path(outdir, "ppi_edges_cleaned.csv")             # [AI]
ppi_pdf_file <- file.path(outdir, "ppi.pdf")                            # [AI]
ppi_png_file <- file.path(outdir, "ppi.png")                            # [AI]
```

## 参数速查

物种：

```text
9606 = human, 默认先用
10090 = mouse
10116 = rat
```

默认按 human 运行；只有当基因映射率明显偏低或上下文提示不是人类时，再考虑改成 mouse/rat。

STRING score：

```text
150 = low
400 = medium, 默认
700 = high, 适合作图展示
900 = highest
```

过滤：

```text
min_node_degree = 1  网络完整
min_node_degree = 2  去掉低连接节点
min_node_degree >= 3 更严格
top_degree_n = NA    展示所有过滤后节点
top_degree_n = 50    只展示连接度最高的前 50 个节点
```

## 输入和清洗

- 候选基因可来自 Venn、DEG、WGCNA、机器学习或手动列表。
- 优先按 `possible_gene_cols` 自动识别基因列。
- 找不到标准列名时，可按无表头读取，选择非空值最多的一列。
- 清洗基因名：去空值、去重复、去明显异常值。
- PPI 边表统一为：`protein_a`, `protein_b`, `score`。

## PPI 过滤

推荐顺序：

1. 下载 STRING 时使用 `string_min_score`。
2. 清洗边表。
3. 用 `min_node_degree` 过滤节点。
4. 如果 `top_degree_n` 不是 `NA`，再保留 degree 最高的前 N 个节点。

top hub 核心逻辑：

```r
degree_table <- table(c(filtered_edges$protein_a, filtered_edges$protein_b))
top_nodes <- names(head(sort(degree_table, decreasing = TRUE), top_degree_n))
filtered_edges <- filtered_edges %>%
  dplyr::filter(protein_a %in% top_nodes, protein_b %in% top_nodes)
```

新增 `top_degree_n` 时，同步更新配置区、`filter_ppi_by_degree()`、`plot_ppi_chord()`、`save_ppi_plots()` 和 main workflow。

注意：`top_degree_n` 主要用于作图展示，完整 STRING 结果和清洗后边表仍应保存。

## 作图和输出

- 使用 `circlize::chordDiagram()` 绘制弦图。
- 节点颜色表示 degree，连线颜色表示 STRING score。
- 用 `ComplexHeatmap::Legend()` 添加图例。

标准输出：

```text
03_PPI/string_interactions.tsv
03_PPI/ppi_edges_cleaned.csv
03_PPI/ppi.pdf
03_PPI/ppi.png
03_PPI/candidate_genes_cleaned.csv
```

## 图太乱时

按这个顺序调整：

1. `string_min_score <- 700`，必要时 `900`。
2. `min_node_degree <- 2`，必要时 `3`。
3. `top_degree_n <- 50`。

## 常见报错

- 没有候选基因：检查文件路径、基因列名、表头和空值。
- STRING 映射失败：检查 gene symbol、物种编号和网络。
- degree 过滤后无边：降低 `min_node_degree` 或 `string_min_score`。
- top hub 过滤后无边：增大 `top_degree_n` 或设为 `NA`。
- Cairo/字体报错：换 `pdf()` / `png()` 或 `plot_family <- "serif"`。

## 回复规范

- 用户要完整代码：给完整 R 脚本，保留配置区和 `# [AI]`，不要擅自改生物学参数。
- 用户要解释代码：说明输入、输出、主流程和可安全修改的参数。
- 重点解释：`string_species_id`, `string_min_score`, `min_node_degree`, `top_degree_n`。
