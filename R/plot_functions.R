# 加载必要的库----------
library(ggplot2)
library(dplyr)
library(rlang)

# basic--------------------
# 基础颜色，直接调用的全局变量
low_color = '#76c07c' # #76c07c
high_color = '#f77680'# #f77680
mid_color = '#fff6f6'
no_sig_color = '#d9d9d9'  
single_color = '#abddf1'

# theme for text color and size
my_theme <-   theme(axis.text = element_text(colour = 'black', size = 12),
                    axis.title = element_text(colour = 'black', size = 15),
                    plot.title = element_text(colour = 'black', size = 16),
                    legend.text = element_text(color = 'black',size = 12),
                    legend.title = element_text(colour = 'black', size = 15))



# single_dataset-----------------
# 箱式图和表格
generate_gene_expression_plot <- function(GSE, given_gene) {
  
  
  # 获取数据
  qc_tpm_loged <- tpm_loged_list[[GSE]]
  qc_sample_info <- sample_info_list[[GSE]]
  deg_list <- deg_list_all[[GSE]]
  
  # 生成箱式图数据
  df_boxplot <- qc_tpm_loged[given_gene, , drop = FALSE] %>%
    t() %>%
    data.frame() %>%
    mutate(Group = qc_sample_info$Group)
  
  # 绘制箱式图
  p <- ggplot(df_boxplot, aes(x = Group, y = !!sym(given_gene), fill = Group)) +
    geom_boxplot() +
    geom_dotplot(binaxis = 'y', stackdir = 'center') +
    scale_fill_brewer(palette = "Blues") +
    theme_bw() +
    labs(x = '', y = paste0(given_gene, " expression (log TPM)"), title = GSE) +
    theme(legend.position = 'none', plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 30, hjust = 1)) +
    my_theme
  
  # 保存图形（可选）
  # ggsave(filename = paste0(GSE, "_", given_gene, "_boxplot.png"), plot = p)
  
  # 打印P值和logFC
  gene_pval <- list()
  gene_pval[["Note"]] <- "DESeq2 is used for differential expression analysis."  # 添加说明
  
  if (is.null(deg_list) || length(deg_list) == 0) {
    gene_pval[["Result"]] <- "No differential analysis results available."
  } else {
    for (i in names(deg_list)) {
      df <- deg_list[[i]]
      if (given_gene %in% df$Gene) {
        gene_pval[[i]] <- paste0('P.Value = ', df$pvalue[df$Gene == given_gene], 
                                 ', logFC = ', df$log2FoldChange[df$Gene == given_gene])
      } else {
        gene_pval[[i]] <- paste0(given_gene, ' not in the DEA results')
      }
    }
  }
  
  
  # 提取GSE信息
  GSE_info_df <- sample_info_all_final[sample_info_all_final$GSE == GSE, ]
  
  # 返回结果
  return(list(gene_boxplot = p, gene_pval = gene_pval, GSE_info_df = GSE_info_df))
}

# single_class-------------------

# bulk
# 定义函数 plot_bar_diff_median-----------------------
plot_bar_diff_median <- function(median_diff_results, given_gene, 
                                 Class, sort_order = "desc") {
  
  # 修改标题和Y轴标签，确保只有第一个单词首字母大写
  title_text <- paste("Difference from first group for", given_gene, "in", Class, "class")
  y_axis_label <- "Difference from first group (Median)"  # Y轴标签
  
  # 根据 sort_order 设置 decreasing 参数
  decreasing_order <- ifelse(sort_order == "desc", FALSE, TRUE)
  
  # 绘制条形图
  p <- ggplot(median_diff_results, aes(x = reorder(GSE, Difference_from_First_Group, 
                                                   decreasing = decreasing_order), 
                                       y = Difference_from_First_Group)) +
    geom_bar(stat = "identity", fill = "steelblue") +  # 条形图
    geom_text(aes(label = sprintf("%.2f", Difference_from_First_Group)),  # 显示两位小数
              hjust = ifelse(median_diff_results$Difference_from_First_Group < 0, 
                             -0.5, 1.1),  # 负值靠右，正值靠左
              colour = "black") +  # 文本颜色为黑色
    coord_flip() +  # 将条形图水平显示，便于阅读
    labs(title = title_text,  # 使用修改后的标题
         x = "GSE",
         y = y_axis_label) +  # 修改Y轴标签
    theme_minimal() +  # 使用简约主题
    theme(axis.text = element_text(colour = 'black'),  # 坐标轴文字颜色
          plot.title = element_text(hjust = 0.5))  # 标题居中
  
  return(p)
}
# 定义函数 plot_boxplot_by_gse_list----------------
plot_boxplot_by_gse_list <- function(median_diff_results, given_gene) {
  # 初始化一个空列表来存储绘图结果
  plot_list <- list()
  
  # 获取 GSE 列表并按照 median_diff_results 中的顺序处理
  gse_list <- median_diff_results$GSE
  
  # 遍历每个 GSE 号，提取数据并绘制箱式图
  for (GSE in gse_list) {
    # 提取对应的表达数据和分组信息
    expression_data <- tpm_loged_list[[GSE]][given_gene, , drop = FALSE]  # 获取指定基因的表达数据
    sample_info <- sample_info_list[[GSE]]  # 获取对应的样本信息
    
    # 确保数据中存在该基因
    if (!is.null(expression_data)) {
      # 创建一个数据框，将表达数据和分组信息合并
      plot_data <- data.frame(Expression = as.numeric(expression_data),
                              Group = sample_info$Group)
      
      # 绘制箱式图和点图
      p <- ggplot(plot_data, aes(x = Group, y = Expression)) +
        geom_boxplot(aes(fill = Group), show.legend = FALSE, outlier.shape = NA) +  # 绘制箱式图，隐藏异常点
        geom_dotplot(binaxis = 'y', stackdir = 'center',
                     color = 'black',
                     dotsize = 1,alpha = 0.7) +  # 增加点图
        scale_fill_brewer(palette = "Set2") +  # 使用更好看的配色方案
        labs(title = paste("Gene:", given_gene, "GSE:", GSE),
             x = "",
             y = "Expression (log TPM)") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # X轴标签倾斜，方便阅读
        theme(axis.text = element_text(colour = 'black'))
      
      # 将绘图结果保存到列表中
      plot_list[[GSE]] <- p
    } else {
      message(paste("Gene", given_gene, "not found in GSE", GSE))
    }
  }
  
  return(plot_list)
}

# array 
plot_array_bar_diff_median <- function(median_diff_results, given_gene, 
                                       Class, sort_order = "desc") {
  
  # 修改标题和Y轴标签，确保只有第一个单词首字母大写
  title_text <- paste("Difference from first group for", given_gene, "in", Class, "class")
  y_axis_label <- "Difference from first group (Median)"  # Y轴标签
  
  # 根据 sort_order 设置 decreasing 参数
  decreasing_order <- ifelse(sort_order == "desc", FALSE, TRUE)
  
  # 绘制条形图
  p <- ggplot(median_diff_results, aes(x = reorder(GSE, Difference_from_First_Group, 
                                                   decreasing = decreasing_order), 
                                       y = Difference_from_First_Group)) +
    geom_bar(stat = "identity", fill = "steelblue") +  # 条形图
    geom_text(aes(label = sprintf("%.2f", Difference_from_First_Group)),  # 显示两位小数
              hjust = ifelse(median_diff_results$Difference_from_First_Group < 0, 
                             -0.5, 1.1),  # 负值靠右，正值靠左
              colour = "black") +  # 文本颜色为黑色
    coord_flip() +  # 将条形图水平显示，便于阅读
    labs(title = title_text,  # 使用修改后的标题
         x = "GSE",
         y = y_axis_label) +  # 修改Y轴标签
    theme_minimal() +  # 使用简约主题
    theme(axis.text = element_text(colour = 'black'),  # 坐标轴文字颜色
          plot.title = element_text(hjust = 0.5))  # 标题居中
  
  return(p)
}

plot_array_boxplot_by_gse_list <- function(median_diff_results, 
                                           given_gene) {
  # 初始化一个空列表来存储绘图结果
  plot_list <- list()
  
  # 获取 GSE 列表并按照 median_diff_results 中的顺序处理
  gse_list <- median_diff_results$GSE
  
  # 遍历每个 GSE 号，提取数据并绘制箱式图
  for (GSE in gse_list) {
    # 提取对应的表达数据和分组信息
    expression_data <- array_exp_list[[GSE]][given_gene, , drop = FALSE]  # 获取指定基因的表达数据
    sample_info <- array_info_list[[GSE]]  # 获取对应的样本信息
    
    # 确保数据中存在该基因
    if (!is.null(expression_data)) {
      # 创建一个数据框，将表达数据和分组信息合并
      plot_data <- data.frame(Expression = as.numeric(expression_data),
                              Group = sample_info$Group)
      
      # 绘制箱式图和点图
      p <- ggplot(plot_data, aes(x = Group, y = Expression)) +
        geom_boxplot(aes(fill = Group), show.legend = FALSE, outlier.shape = NA) +  # 绘制箱式图，隐藏异常点
        geom_dotplot(binaxis = 'y', stackdir = 'center',
                     color = 'black',
                     dotsize = 1,alpha = 0.7) +  # 增加点图
        scale_fill_brewer(palette = "Set2") +  # 使用更好看的配色方案
        labs(title = paste("Gene:", given_gene, "GSE:", GSE),
             x = "",
             y = "Expression (log level)") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # X轴标签倾斜，方便阅读
        theme(axis.text = element_text(colour = 'black'))
      
      # 将绘图结果保存到列表中
      plot_list[[GSE]] <- p
    } else {
      message(paste("Gene", given_gene, "not found in GSE", GSE))
    }
  }
  
  return(plot_list)
}





# all_datasets ---------------------
# 函数：汇总基因图-----------------
plot_gene_summary <- function(gene_summary_result,
                              total_studies = total_studies,
                              filter_threshold = filter_threshold
                              ) {
  
  gene_summary_plot_list <- list()  # 创建一个列表存储绘图结果
  
  # 计算总的研究数和cutoff
  # total_studies <- as.numeric(length(correlation_results))
  cutoff <- floor(total_studies * filter_threshold)
  
  # 绘制直方图1: 过滤前
  p0 <- gene_summary_result %>% 
    ggplot(aes(x = Count, fill = Filter)) +
    geom_histogram(binwidth = 1, color = "black") +
    geom_vline(xintercept = cutoff + 0.5, linetype = "dashed", color = "red", size = 1, alpha = 0.5) +
    scale_fill_manual(values = c("Filter" = "grey", "Keep" = single_color)) +
    annotate("text", x = cutoff + 0.5, y = Inf, label = paste("Count = ", cutoff), 
             vjust = 2, hjust = -0.1, color = "red") +
    labs(x = "Dataset count", y = "Gene frequency") +
    theme_minimal() + 
    theme(legend.position = "none") +
    my_theme
  
  gene_summary_plot_list[["Histogram_Count"]] <- p0  # 保存图形
  
  # 绘制Sample_Count_Sum直方图
  p1 <- gene_summary_result %>% 
    filter(Filter == 'Keep') %>% 
    ggplot(aes(x = Total_Sample_Count)) +
    geom_histogram(binwidth = 50, position = "identity", fill = single_color, color = 'black') +
    xlab("Total sample count") +
    ylab("Gene frequency") +
    theme_minimal() + 
    theme(legend.position = "none") +
    my_theme
  
  gene_summary_plot_list[["Histogram_Total_Sample_Count"]] <- p1  # 保存图形
  
  # 绘制Weight直方图
  p2 <- gene_summary_result %>% 
    filter(Filter == 'Keep') %>% 
    ggplot(aes(x = Weight, fill = Confidence)) +
    geom_histogram(binwidth = 0.05, color = 'black', size = 0.1) +
    scale_fill_manual(values = c("HC" = high_color, "LC" = low_color)) +
    labs(x = "Weight", y = "Gene frequency") +
    theme_minimal() + 
    my_theme +
    theme(legend.position = c(0.8, 0.8))
  
  gene_summary_plot_list[["Histogram_Weight"]] <- p2  # 保存图形
  
  # 绘制Confidence直条图
  p3 <- gene_summary_result %>% 
    filter(Filter == 'Keep') %>% 
    ggplot(aes(x = Confidence, fill = Confidence)) +
    geom_bar(color = "black") +
    geom_text(stat = 'count', aes(label = ..count..), vjust = -0.5, color = "black", size = 5) +
    scale_fill_manual(values = c("HC" = high_color, "LC" = low_color)) +
    labs(x = "Confidence", y = "Gene frequency") +
    theme_minimal() + 
    my_theme +
    ylim(NA, max(table(gene_summary_result$Confidence)) + 200) +
    theme(legend.position = "none")
  
  gene_summary_plot_list[["Barplot_Confidence"]] <- p3  # 保存图形
  
  # 绘制top 20的棒棒糖图
  sorted_data <- gene_summary_result %>%
    arrange(desc(Weight)) %>%
    slice(1:20)  # 选择前20个基因
  
  p4 <- ggplot(sorted_data, aes(x = reorder(Gene, Weight), y = Weight)) +
    geom_segment(aes(xend = Gene, yend = 0), color = "#6565fe") +
    geom_point(size = 4, color = "#6565fe") +
    coord_flip() +
    theme_bw() +
    xlab("") +
    ylab("Weight") +
    scale_y_continuous(limits = c(0, NA)) +
    theme(axis.text.y = element_text(face = "italic"),
          plot.title = element_text(hjust = 0.5)) +
    my_theme
  
  gene_summary_plot_list[["Top20_Lollipop"]] <- p4  # 保存图形
  
  return(gene_summary_plot_list)  # 返回绘图结果列表
}

# array
array_plot_gene_summary <- function(gene_summary_result,
                                    total_studies = total_studies,
                                    filter_threshold = filter_threshold) {
  
  gene_summary_plot_list <- list()  # 创建一个列表存储绘图结果
  
  # 计算总的研究数和cutoff
  cutoff <- floor(total_studies * filter_threshold)
  
  # 绘制直方图1: 过滤前
  p0 <- gene_summary_result %>% 
    ggplot(aes(x = Count, fill = Filter)) +
    geom_histogram(binwidth = 1, color = "black") +
    geom_vline(xintercept = cutoff + 0.5, linetype = "dashed", color = "red", size = 1, alpha = 0.5) +
    scale_fill_manual(values = c("Filter" = "grey", "Keep" = single_color)) +
    annotate("text", x = cutoff + 0.5, y = Inf, label = paste("Count = ", cutoff), 
             vjust = 2, hjust = -0.1, color = "red") +
    labs(x = "Dataset count", y = "Gene frequency") +
    theme_minimal() + 
    theme(legend.position = "none") +
    my_theme
  
  gene_summary_plot_list[["Histogram_Count"]] <- p0  # 保存图形
  
  # 绘制Sample_Count_Sum直方图
  p1 <- gene_summary_result %>% 
    filter(Filter == 'Keep') %>% 
    ggplot(aes(x = Total_Sample_Count)) +
    geom_histogram(binwidth = 50, position = "identity", fill = single_color, color = 'black') +
    xlab("Total sample count") +
    ylab("Gene frequency") +
    theme_minimal() + 
    theme(legend.position = "none") +
    my_theme
  
  gene_summary_plot_list[["Histogram_Total_Sample_Count"]] <- p1  # 保存图形
  
  # 绘制Weight直方图
  p2 <- gene_summary_result %>% 
    filter(Filter == 'Keep') %>% 
    ggplot(aes(x = Weight, fill = Confidence)) +
    geom_histogram(binwidth = 0.05, color = 'black', size = 0.1) +
    scale_fill_manual(values = c("HC" = high_color, "LC" = low_color)) +
    labs(x = "Weight", y = "Gene frequency") +
    theme_minimal() + 
    my_theme +
    theme(legend.position = c(0.8, 0.8))
  
  gene_summary_plot_list[["Histogram_Weight"]] <- p2  # 保存图形
  
  # 绘制Confidence直条图
  p3 <- gene_summary_result %>% 
    filter(Filter == 'Keep') %>% 
    ggplot(aes(x = Confidence, fill = Confidence)) +
    geom_bar(color = "black") +
    geom_text(stat = 'count', aes(label = ..count..), vjust = -0.5, color = "black", size = 5) +
    scale_fill_manual(values = c("HC" = high_color, "LC" = low_color)) +
    labs(x = "Confidence", y = "Gene frequency") +
    theme_minimal() + 
    my_theme +
    ylim(NA, max(table(gene_summary_result$Confidence)) + 200) +
    theme(legend.position = "none")
  
  gene_summary_plot_list[["Barplot_Confidence"]] <- p3  # 保存图形
  
  # 绘制top 20的棒棒糖图
  sorted_data <- gene_summary_result %>%
    arrange(desc(Weight)) %>%
    slice(1:20)  # 选择前20个基因
  
  p4 <- ggplot(sorted_data, aes(x = reorder(Gene, Weight), y = Weight)) +
    geom_segment(aes(xend = Gene, yend = 0), color = "#6565fe") +
    geom_point(size = 4, color = "#6565fe") +
    coord_flip() +
    theme_bw() +
    xlab("") +
    ylab("Weight") +
    scale_y_continuous(limits = c(0, NA)) +
    theme(axis.text.y = element_text(face = "italic"),
          plot.title = element_text(hjust = 0.5)) +
    my_theme
  
  gene_summary_plot_list[["Top20_Lollipop"]] <- p4  # 保存图形
  
  return(gene_summary_plot_list)  # 返回绘图结果列表
}


# 函数：富集分析结果条图，前10------------------
plot_enrich_bar_list <- function(ora_results, 
                                 n = 10, 
                                 to_sentence = TRUE, 
                                 High = '#2171b5', 
                                 Low = '#9ecae1') {
  
  plot_list <- list()  # 创建一个列表存储绘图结果
  
  # 定义内部绘图函数
  plot_enrich_bar <- function(em_df, n, to_sentence, database, High, Low) {
    if (nrow(em_df) == 0) {
      return(ggplot() + labs(title = paste0('No significant pathways found in ', database, ' analysis')))  # 返回一个空图并提示没有找到通路
    }
    
    pathway <- em_df[1:min(n, nrow(em_df)), ]  # 选择最多n行
    pathway <- pathway[pathway$p.adjust < 0.05, ]  # 筛选p.adjust小于0.05的通路
    
    if (nrow(pathway) == 0) {
      return(ggplot() + labs(title = paste0('No significant pathways after filtering in ', database, ' analysis')))  # 返回一个提示没有显著通路的空图
    }
    
    if (to_sentence) {
      pathway$Description <- gsub('_', ' ', pathway$Description)
      pathway$Description <- str_to_sentence(pathway$Description)
    }
    
    porder <- factor(1:nrow(pathway), labels = rev(pathway$Description))
    
    ggplot(pathway, aes(x = Description, y = -1 * log10(p.adjust))) +
      geom_bar(stat = "identity", aes(x = rev(porder), fill = -1 * log10(p.adjust))) +
      scale_fill_gradient(low = Low, high = High) +
      theme_bw() +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 40)) +
      coord_flip() +
      theme(legend.position = "none",
            axis.text.x = element_text(colour = 'black', size = 12),
            axis.text.y = element_text(colour = 'black', size = 10),
            plot.title = element_text(hjust = 0.5)) +
      labs(x = '', y = '-Log10 (p.adjust)', title = paste0('Top ', min(n, nrow(em_df)), ' pathways in\n', database, ' analysis'))
  }
  
  # 绘制条形图并存储在列表中
  for (gmt in names(ora_results)) {
    if (nrow(ora_results[[gmt]]) > 0) {
      plot_list[[gmt]] <- plot_enrich_bar(ora_results[[gmt]], n, to_sentence, database = gmt, High = High, Low = Low)
    } else {
      plot_list[[gmt]] <- ggplot() + labs(title = paste0('No significant pathways found in ', gmt, ' analysis'))  # 没有数据时的空图
    }
  }
  
  return(plot_list)  # 返回包含所有绘图结果的列表
}


# array 
array_plot_enrich_bar_list <- function(ora_results, 
                                       n = 10, 
                                       to_sentence = TRUE, 
                                       High = '#2171b5', 
                                       Low = '#9ecae1') {
  
  plot_list <- list()  # 创建一个列表存储绘图结果
  
  # 定义内部绘图函数
  array_plot_enrich_bar <- function(em_df, n, to_sentence, database, High, Low) {
    if (nrow(em_df) == 0) {
      return(ggplot() + labs(title = paste0('No significant pathways found in ', database, ' analysis')))  # 返回空图并提示没有找到通路
    }
    
    pathway <- em_df[1:min(n, nrow(em_df)), ]  # 选择最多n行
    pathway <- pathway[pathway$p.adjust < 0.05, ]  # 筛选p.adjust小于0.05的通路
    
    if (nrow(pathway) == 0) {
      return(ggplot() + labs(title = paste0('No significant pathways after filtering in ', database, ' analysis')))  # 返回一个提示没有显著通路的空图
    }
    
    if (to_sentence) {
      pathway$Description <- gsub('_', ' ', pathway$Description)
      pathway$Description <- str_to_sentence(pathway$Description)
    }
    
    porder <- factor(1:nrow(pathway), labels = rev(pathway$Description))
    
    ggplot(pathway, aes(x = Description, y = -1 * log10(p.adjust))) +
      geom_bar(stat = "identity", aes(x = rev(porder), fill = -1 * log10(p.adjust))) +
      scale_fill_gradient(low = Low, high = High) +
      theme_bw() +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 40)) +
      coord_flip() +
      theme(legend.position = "none",
            axis.text.x = element_text(colour = 'black', size = 12),
            axis.text.y = element_text(colour = 'black', size = 10),
            plot.title = element_text(hjust = 0.5)) +
      labs(x = '', y = '-Log10 (p.adjust)', title = paste0('Top ', min(n, nrow(em_df)), ' pathways in\n', database, ' analysis'))
  }
  
  # 绘制条形图并存储在列表中
  for (gmt in names(ora_results)) {
    if (nrow(ora_results[[gmt]]) > 0) {
      plot_list[[gmt]] <- array_plot_enrich_bar(ora_results[[gmt]], n, to_sentence, database = gmt, High = High, Low = Low)
    } else {
      plot_list[[gmt]] <- ggplot() + labs(title = paste0('No significant pathways found in ', gmt, ' analysis'))  # 没有数据时的空图
    }
  }
  
  return(plot_list)  # 返回包含所有绘图结果的列表
}



# altas_plot-----------
# 封装 plot_altas 函数
plot_altas <- function(Gene, V1, V2) {
  # 第一部分：生成四个图像列表用于第二行轮播
  plot_list_1 <- list()
  
  # 绘制 CellDimPlot
  plot_list_1[[1]] <- CellDimPlot(muscle_core_altas,
                                  group.by = V1, 
                                  reduction = "umap",
                                  split.by = V2,
                                  label = TRUE, 
                                  label_insitu = TRUE, 
                                  label_repel = TRUE, 
                                  label.size = 4,
                                  label.fg = "black",
                                  label.bg = "white",
                                  label_segment_color = "black",
                                  cells.highlight = TRUE,
                                  theme_use = "theme_blank", 
                                  legend.position = "none"
  )
  
  # 绘制 CellStatPlot
  plot_list_1[[2]] <- CellStatPlot(
    srt = muscle_core_altas,
    stat.by = V1, 
    group.by = V2, 
    split.by = NULL
  )
  
  # 绘制 FeatureDimPlot
  plot_list_1[[3]] <- FeatureDimPlot(muscle_core_altas, 
                                     features = Gene, 
                                     reduction = "umap",
                                     split.by = V2, 
                                     label = TRUE,
                                     cells.highlight = TRUE,
                                     theme_use = 'theme_blank'
  )
  
  # 绘制 FeatureStatPlot
  plot_list_1[[4]] <- FeatureStatPlot(muscle_core_altas, 
                                      stat.by = Gene, 
                                      group.by = V1, 
                                      split.by = NULL, 
                                      legend.position = 'none',
                                      xlab = NULL
  )
  
  # 第二部分：生成单个研究的图像列表用于第三行轮播
  data_names <- unique(muscle_core_altas$Data_name)
  plot_list_2 <- list()
  
  for (Data_name in data_names) {
    # 子集数据
    altas_sub <- muscle_core_altas[, muscle_core_altas$Data_name == Data_name]
    altas_sub$Group <- droplevels(altas_sub$Group)
    
    # 绘制 FeatureStatPlot
    p <- FeatureStatPlot(
      altas_sub,  
      stat.by = Gene,  
      group.by = V1,  
      split.by = "Group",  
      title = Data_name  
    )
    
    # 保存到列表
    plot_list_2[[Data_name]] <- p
  }
  
  return(list(plot_list_1, plot_list_2))  # 返回两个列表
}

# gtex---------------
# 趋势分析----------------
# 定义绘图函数
plot_trend <- function(data, gene, center, title, fill_color) {
  if (center == "mean") {
    plot <- ggplot(data, aes(x = AGE, y = Expression)) +
      stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "pointrange", color = fill_color, size = 0.8) +
      geom_dotplot(aes(fill = AGE), binaxis = 'y', stackdir = 'center', dotsize = 0.4, position = position_dodge(width = 0.3)) +
      theme_bw() +
      scale_fill_brewer(palette = ifelse(fill_color == '#00bfc4', "Blues", "Reds")) +
      labs(x = 'Age group', y = paste0(gene, ' expression (TPM)'), title = title) +
      my_theme +
      theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))
  } else if (center == "median") {
    plot <- ggplot(data, aes(x = AGE, y = Expression)) +
      stat_summary(fun = median, geom = "pointrange", fun.min = function(x) quantile(x, 0.25), fun.max = function(x) quantile(x, 0.75), color = fill_color, size = 0.8) +
      geom_dotplot(aes(fill = AGE), binaxis = 'y', stackdir = 'center', dotsize = 0.4, position = position_dodge(width = 0.3)) +
      theme_bw() +
      scale_fill_brewer(palette = ifelse(fill_color == '#00bfc4', "Blues", "Reds")) +
      labs(x = 'Age group', y = paste0(gene, ' expression (TPM)'), title = title) +
      my_theme +
      theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))
  } else {
    stop("Invalid value for 'center'. Choose 'mean' or 'median'.")
  }
  return(plot)
}
# 定义趋势分析函数，增加 seed 参数
gtex_trend_analysis_with_plots <- function(gene, gtex_muscle_info, gtex_muscle_tpm, 
                                      age_group = c('40-49', '50-59', '60-69'), 
                                      center = "median", 
                                      seed = 2024) {
  
  library(ggplot2)
  library(dplyr)
  library(clinfun)
  
  # 设置随机种子
  set.seed(seed = seed)
  
  # Step 1: 检查基因是否存在
  if (!(gene %in% colnames(gtex_muscle_tpm))) {
    stop(paste("Gene", gene, "not found in gtex_muscle_tpm data."))
  }
  
  # Step 2: 合并 gtex_muscle_info 和基因表达数据
  expression_data <- gtex_muscle_info %>%
    mutate(Expression = gtex_muscle_tpm[, gene])
  
  # Step 3: 将年龄组设置为有序因子，保留所有年龄组用于绘图
  expression_data$AGE <- factor(expression_data$AGE, ordered = TRUE)
  
  # Step 4: 分别提取男性（SEX == 1）和女性（SEX == 2）的数据
  male_data <- expression_data %>% filter(SEX == 1)
  female_data <- expression_data %>% filter(SEX == 2)
  
  # Step 5: 过滤掉缺失值的行，确保数据完整
  male_data <- male_data %>% filter(!is.na(Expression) & !is.na(AGE))
  female_data <- female_data %>% filter(!is.na(Expression) & !is.na(AGE))
  
  # Step 6: 执行趋势检验仅针对指定的 age_group
  male_trend_result <- NULL
  female_trend_result <- NULL
  
  if (!is.null(age_group)) {
    # 筛选指定年龄组数据
    male_data_for_test <- male_data %>% filter(AGE %in% age_group)
    female_data_for_test <- female_data %>% filter(AGE %in% age_group)
    
    if (nrow(male_data_for_test) > 0) {
      male_trend_result <- jonckheere.test(male_data_for_test$Expression, male_data_for_test$AGE, nperm = 1000)
    } else {
      warning("No valid data for male samples in the selected age group.")
    }
    
    if (nrow(female_data_for_test) > 0) {
      female_trend_result <- jonckheere.test(female_data_for_test$Expression, female_data_for_test$AGE, nperm = 1000)
    } else {
      warning("No valid data for female samples in the selected age group.")
    }
  }
  
  # Step 7: 绘制趋势图
  male_plot <- plot_trend(male_data, gene, center, 'Male samples', '#00bfc4')
  female_plot <- plot_trend(female_data, gene, center, 'Female samples', '#f8766d')
  
  # Step 8: 返回结果列表
  return(list(
    male_plot = male_plot,
    female_plot = female_plot,
    male_trend_result = male_trend_result,
    female_trend_result = female_trend_result
  ))
}

# 棒棒糖图------------
plot_gtex_top20_lollipop <- function(result_df) {
  # 对结果根据相关系数 r 进行排序并选择前 20 个基因
  sorted_data <- result_df %>%
    arrange(desc(r)) %>%
    slice(1:20)
  
  # 绘制棒棒糖图
  p <- ggplot(sorted_data, aes(x = reorder(Gene, r), y = r)) +
    geom_segment(aes(xend = Gene, yend = 0), color = "#6565fe") +
    geom_point(size = 4, color = "#6565fe") +
    coord_flip() +
    theme_bw() +
    xlab("") +
    ylab("Correlation coefficient (r)") +
    scale_y_continuous(limits = c(0, NA)) +
    theme(axis.text.y = element_text(face = "italic"),
          plot.title = element_text(hjust = 0.5)) +
    my_theme+
    ggtitle("Top 20 genes by correlation (r)")  # 添加标题
  
  return(p)
}

# enrich bar--------
plot_gtex_enrich_bar_list <- function(ora_results, 
                                      n = 10, 
                                      to_sentence = TRUE, 
                                      High = '#2171b5', 
                                      Low = '#9ecae1') {
  
  plot_list <- list()  # 创建一个列表存储绘图结果
  
  # 定义内部绘图函数
  plot_enrich_bar <- function(em_df, n, to_sentence, database, High, Low) {
    if (nrow(em_df) == 0) {
      return(ggplot() + labs(title = paste0('No significant pathways found in ', database, ' analysis')))  # 返回空图并提示没有找到通路
    }
    
    pathway <- em_df[1:min(n, nrow(em_df)), ]  # 选择最多n行
    pathway <- pathway[pathway$p.adjust < 0.05, ]  # 筛选p.adjust小于0.05的通路
    
    if (nrow(pathway) == 0) {
      return(ggplot() + labs(title = paste0('No significant pathways after filtering in ', database, ' analysis')))  # 返回一个提示没有显著通路的空图
    }
    
    if (to_sentence) {
      pathway$Description <- gsub('_', ' ', pathway$Description)
      pathway$Description <- str_to_sentence(pathway$Description)
    }
    
    porder <- factor(1:nrow(pathway), labels = rev(pathway$Description))
    
    ggplot(pathway, aes(x = Description, y = -1 * log10(p.adjust))) +
      geom_bar(stat = "identity", aes(x = rev(porder), fill = -1 * log10(p.adjust))) +
      scale_fill_gradient(low = Low, high = High) +
      theme_bw() +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 40)) +
      coord_flip() +
      theme(legend.position = "none",
            axis.text.x = element_text(colour = 'black', size = 12),
            axis.text.y = element_text(colour = 'black', size = 10),
            plot.title = element_text(hjust = 0.5)) +
      labs(x = '', y = '-Log10 (p.adjust)', title = paste0('Top ', min(n, nrow(em_df)), ' pathways in\n', database, ' analysis'))
  }
  
  # 绘制条形图并存储在列表中
  for (gmt in names(ora_results)) {
    if (nrow(ora_results[[gmt]]) > 0) {
      plot_list[[gmt]] <- plot_enrich_bar(ora_results[[gmt]], n, to_sentence, database = gmt, High = High, Low = Low)
    } else {
      plot_list[[gmt]] <- ggplot() + labs(title = paste0('No significant pathways found in ', gmt, ' analysis'))  # 没有数据时的空图
    }
  }
  
  return(plot_list)  # 返回包含所有绘图结果的列表
}











