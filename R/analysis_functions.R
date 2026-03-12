# single_class----------------------
# 组间中位数差值
# bulk 
calculate_group_median_difference_by_class <- function(Class, 
                                                       given_gene, 
                                                       sort_order = "desc") {
  # 加载必要的库
  library(dplyr)
  
  # 获取该Class下的所有GSE
  class_gse <- sample_info_all_final$GSE[sample_info_all_final$Class == Class] %>% unique()
  
  # 初始化结果列表
  results <- list()
  
  # 遍历每个GSE
  for (GSE in class_gse) {
    # 从tpm_loged_list和sample_info_list获取矩阵和样本信息
    tpm_loged_gse <- tpm_loged_list[[GSE]]
    sample_info_gse <- sample_info_list[[GSE]]
    
    # 确保给定的基因在表达矩阵中
    if (!(given_gene %in% rownames(tpm_loged_gse))) {
      next  # 如果基因不在矩阵中，则跳过
    }
    
    # 获取分组信息
    group_info <- sample_info_gse$Group
    
    # 计算每个分组的表达中位数
    expression_data <- tpm_loged_gse[given_gene, , drop = TRUE]
    expression_df <- data.frame(Expression = expression_data, Group = group_info)
    
    median_results <- expression_df %>%
      group_by(Group) %>%
      summarize(Median_Expression = median(Expression, na.rm = TRUE), .groups = 'drop')
    
    # 获取第一组的中位数
    first_group_value <- median_results$Median_Expression[1]
    
    # 计算每组中位数与第一组中位数的差值
    diff_values <- median_results$Median_Expression - first_group_value
    median_results$Difference_from_First_Group <- diff_values
    
    # 找到绝对值最大的差值及其对应的分组
    max_diff_row <- median_results %>%
      mutate(Abs_Difference = abs(Difference_from_First_Group)) %>%
      filter(Abs_Difference == max(Abs_Difference)) %>%
      select(Group, Median_Expression, Difference_from_First_Group) %>%
      slice(1)  # 选择第一行，如果有多个最大值则取第一条
    
    # 保存结果
    results[[GSE]] <- c(GSE, max_diff_row)
  }
  
  # 将结果转换为数据框并命名列
  results_df <- do.call(rbind, results)
  colnames(results_df) <- c("GSE", "Group", "Median_Expression", "Difference_from_First_Group")
  
  # 转换为数据框，并将GSE列转换为字符类型
  results_df <- as.data.frame(results_df)
  results_df$GSE <- as.character(results_df$GSE)  # 确保GSE列为字符类型
  results_df$Difference_from_First_Group <- as.numeric(as.character(results_df$Difference_from_First_Group))
  
  # 排序逻辑
  if (sort_order == "asc") {
    results_df <- results_df %>% arrange(Difference_from_First_Group)
  } else {
    results_df <- results_df %>% arrange(desc(Difference_from_First_Group))
  }
  
  # 去重后再将 "Keywords" 和 "Keywords_zh" 添加到结果中
  keywords_info <- sample_info_all_final %>%
    select(GSE, Keywords, Keywords_zh) %>%
    distinct(GSE, .keep_all = TRUE)  # 确保每个 GSE 只保留一条记录
  
  results_df <- results_df %>%
    left_join(keywords_info, by = "GSE")
  
  results_df$Group <- unlist(results_df$Group)
  results_df$Median_Expression <- unlist(results_df$Median_Expression)
  
  return(results_df)
}

# array 
calculate_array_group_median_difference_by_class <- function(Class, 
                                                             given_gene, 
                                                             sort_order = "desc") {
  # 加载必要的库
  library(dplyr)
  
  # 获取该Class下的所有GSE
  class_gse <- array_sample_info_all$GSE[array_sample_info_all$Class == Class] %>% unique()
  
  # 初始化结果列表
  results <- list()
  
  # 遍历每个GSE
  for (GSE in class_gse) {
    # 从array_exp_list和array_info_list获取矩阵和样本信息
    tpm_loged_gse <- array_exp_list[[GSE]]
    sample_info_gse <- array_info_list[[GSE]]
    
    # 确保给定的基因在表达矩阵中
    if (!(given_gene %in% rownames(tpm_loged_gse))) {
      next  # 如果基因不在矩阵中，则跳过
    }
    
    # 获取分组信息
    group_info <- sample_info_gse$Group
    
    print(GSE)
    # 计算每个分组的表达中位数
    expression_data <- tpm_loged_gse[given_gene, , drop = TRUE]
    expression_df <- data.frame(Expression = expression_data, Group = group_info)
    
    median_results <- expression_df %>%
      group_by(Group) %>%
      summarize(Median_Expression = median(Expression, na.rm = TRUE), .groups = 'drop')
    
    # 获取第一组的中位数
    first_group_value <- median_results$Median_Expression[1]
    
    # 计算每组中位数与第一组中位数的差值
    diff_values <- median_results$Median_Expression - first_group_value
    median_results$Difference_from_First_Group <- diff_values
    
    # 找到绝对值最大的差值及其对应的分组
    max_diff_row <- median_results %>%
      mutate(Abs_Difference = abs(Difference_from_First_Group)) %>%
      filter(Abs_Difference == max(Abs_Difference)) %>%
      select(Group, Median_Expression, Difference_from_First_Group) %>%
      slice(1)  # 选择第一行，如果有多个最大值则取第一条
    
    # 保存结果
    results[[GSE]] <- c(GSE, max_diff_row)
  }
  
  # 将结果转换为数据框并命名列
  results_df <- do.call(rbind, results)
  colnames(results_df) <- c("GSE", "Group", "Median_Expression", "Difference_from_First_Group")
  
  # 转换为数据框，并将GSE列转换为字符类型
  results_df <- as.data.frame(results_df)
  results_df$GSE <- as.character(results_df$GSE)  # 确保GSE列为字符类型
  results_df$Difference_from_First_Group <- as.numeric(as.character(results_df$Difference_from_First_Group))
  
  # 排序逻辑
  if (sort_order == "asc") {
    results_df <- results_df %>% arrange(Difference_from_First_Group)
  } else {
    results_df <- results_df %>% arrange(desc(Difference_from_First_Group))
  }
  
  # 去重后再将 "Keywords" 和 "Keywords_zh" 添加到结果中
  keywords_info <- array_sample_info_all %>%
    select(GSE, Keywords, Keywords_zh) %>%
    distinct(GSE, .keep_all = TRUE)  # 确保每个 GSE 只保留一条记录
  
  results_df <- results_df %>%
    left_join(keywords_info, by = "GSE")
  
  results_df$Group <- unlist(results_df$Group)
  results_df$Median_Expression <- unlist(results_df$Median_Expression)
  
  return(results_df)
}




# all_datasets analysis-----------------------
# RNA step1 (fast + robust) ------------------------
calculate_gene_correlations <- function(given_gene, Class = "All",
                                        threshold = 1,
                                        cor_threshold = 0.5,
                                        pvalue_threshold = 0.05) {
  library(dplyr)
  
  # 获取GSE列表
  if (Class == "All") {
    gse_list <- names(tpm_samples_list)
  } else {
    class_gse <- sample_info_all_final$GSE[sample_info_all_final$Class == Class] %>% unique()
    gse_list <- class_gse
  }
  
  results_list <- list()
  
  for (GSE in gse_list) {
    exp_qc <- as.matrix(tpm_samples_list[[GSE]])
    if (is.null(exp_qc) || nrow(exp_qc) == 0 || ncol(exp_qc) < 3) next
    if (is.null(rownames(exp_qc))) next
    
    # 检查 given_gene
    if (!(given_gene %in% rownames(exp_qc))) {
      message(paste("Given gene", given_gene, "not found in GSE", GSE, "- Skipping."))
      next
    }
    
    # 过滤低表达（na.rm + drop=FALSE 防报错）
    exp_qc_high <- exp_qc[rowMeans(exp_qc, na.rm = TRUE) > threshold, , drop = FALSE]
    if (nrow(exp_qc_high) < 2) next
    
    if (!(given_gene %in% rownames(exp_qc_high))) {
      message(paste("Given gene", given_gene, "not found in high expression genes for GSE", GSE, "- Skipping."))
      next
    }
    
    # given gene 向量
    x <- as.numeric(exp_qc_high[given_gene, ])
    ok_x <- !is.na(x)
    if (sum(ok_x) < 3) next
    
    # 候选基因（去掉自身 + 名字过滤）
    genes <- rownames(exp_qc_high)
    genes <- genes[genes != given_gene & nchar(genes) > 1]
    if (length(genes) == 0) next
    
    mat <- exp_qc_high[genes, ok_x, drop = FALSE]
    x2  <- x[ok_x]
    if (ncol(mat) < 3 || nrow(mat) == 0) next
    
    # 每行有效样本数（该行非NA；x2 已无 NA）
    ok_mat <- !is.na(mat)
    n_i <- rowSums(ok_mat)
    keep_n <- which(n_i >= 3)
    if (length(keep_n) == 0) next
    
    mat2 <- mat[keep_n, , drop = FALSE]
    ok2  <- ok_mat[keep_n, , drop = FALSE]
    n2   <- n_i[keep_n]
    gene_names2 <- rownames(mat2)
    
    # ---- 向量化 Pearson r（pairwise：每行按自身非NA）----
    x_centered <- x2 - mean(x2, na.rm = TRUE)
    
    row_means <- rowSums(mat2, na.rm = TRUE) / n2
    mat_centered <- mat2 - row_means
    mat_centered[!ok2] <- 0
    
    cov_num <- as.vector(mat_centered %*% x_centered)
    ss_y <- rowSums(mat_centered^2)
    ss_x <- sum(x_centered^2)
    
    r <- cov_num / sqrt(pmax(1e-12, ss_y) * pmax(1e-12, ss_x))
    r[!is.finite(r)] <- NA_real_
    
    # p值（Pearson t 检验）
    tstat <- r * sqrt((n2 - 2) / pmax(1e-12, 1 - r^2))
    pvals <- 2 * pt(-abs(tstat), df = n2 - 2)
    
    keep <- which(!is.na(r) & r > cor_threshold & !is.na(pvals) & pvals < pvalue_threshold)
    if (length(keep) == 0) {
      message(paste("No significant correlations found for GSE", GSE, "with gene", given_gene, "- Skipping."))
      next
    }
    
    results_df <- data.frame(
      Gene = gene_names2[keep],
      Correlation = as.numeric(r[keep]),
      PValue = as.numeric(pvals[keep]),
      stringsAsFactors = FALSE
    )
    rownames(results_df) <- results_df$Gene
    
    results_list[[GSE]] <- results_df
  }
  
  save(results_list, file = paste0("data/", given_gene, "_", Class, "_cor_list.rdata"))
  results_list
}
# array step1 (fast + robust) --------------
array_calculate_gene_correlations <- function(given_gene, Class = "All",
                                              threshold = 0.5,
                                              cor_threshold = 0.5,
                                              pvalue_threshold = 0.05) {
  library(dplyr)
  
  # 获取GSE列表
  if (Class == "All") {
    gse_list <- names(array_exp_list)
  } else {
    gse_list <- array_sample_info_all$GSE[array_sample_info_all$Class == Class] %>% unique()
  }
  
  results_list <- list()
  
  for (GSE in gse_list) {
    exp_qc <- as.matrix(array_exp_list[[GSE]])
    if (is.null(exp_qc) || nrow(exp_qc) == 0 || ncol(exp_qc) < 3) next
    if (is.null(rownames(exp_qc))) next
    
    # 检查给定基因是否存在
    if (!(given_gene %in% rownames(exp_qc))) {
      message(paste("Given gene", given_gene, "not found in GSE", GSE, "- Skipping."))
      next
    }
    
    # 过滤低表达（na.rm 防止 rowMeans 遇 NA 变 NA；drop=FALSE 防止维度坍塌）
    exp_qc_high <- exp_qc[rowMeans(exp_qc, na.rm = TRUE) > threshold, , drop = FALSE]
    if (nrow(exp_qc_high) < 2) next
    
    if (!(given_gene %in% rownames(exp_qc_high))) {
      message(paste("Given gene", given_gene, "not found in high expression genes for GSE", GSE, "- Skipping."))
      next
    }
    
    # 给定基因表达向量
    x <- as.numeric(exp_qc_high[given_gene, ])
    ok_x <- !is.na(x)
    if (sum(ok_x) < 3) next
    
    # 候选基因矩阵（去掉自身；顺便过滤掉奇怪的空/1字符名）
    genes <- rownames(exp_qc_high)
    genes <- genes[genes != given_gene & nchar(genes) > 1]
    if (length(genes) == 0) next
    
    mat <- exp_qc_high[genes, ok_x, drop = FALSE]
    x2  <- x[ok_x]
    if (ncol(mat) < 3 || nrow(mat) == 0) next
    
    # 每个基因的有效样本数（该行非NA的列数；x2 已无 NA）
    ok_mat <- !is.na(mat)
    n_i <- rowSums(ok_mat)
    keep_n <- which(n_i >= 3)
    if (length(keep_n) == 0) next
    
    mat2 <- mat[keep_n, , drop = FALSE]
    ok2  <- ok_mat[keep_n, , drop = FALSE]
    n2   <- n_i[keep_n]
    gene_names2 <- rownames(mat2)
    
    # ---- 向量化计算 Pearson r（pairwise：按每行自己的非NA）----
    x_centered <- x2 - mean(x2, na.rm = TRUE)
    
    row_means <- rowSums(mat2, na.rm = TRUE) / n2
    mat_centered <- mat2 - row_means
    mat_centered[!ok2] <- 0
    
    cov_num <- as.vector(mat_centered %*% x_centered)
    ss_y <- rowSums(mat_centered^2)
    ss_x <- sum(x_centered^2)
    
    r <- cov_num / sqrt(pmax(1e-12, ss_y) * pmax(1e-12, ss_x))
    r[!is.finite(r)] <- NA_real_
    
    # p值（与 cor.test 的 Pearson 对应：t = r*sqrt((n-2)/(1-r^2))）
    tstat <- r * sqrt((n2 - 2) / pmax(1e-12, 1 - r^2))
    pvals <- 2 * pt(-abs(tstat), df = n2 - 2)
    
    keep <- which(!is.na(r) & r > cor_threshold & !is.na(pvals) & pvals < pvalue_threshold)
    if (length(keep) == 0) {
      message(paste("No significant correlations found for GSE", GSE, "with gene", given_gene, "- Skipping."))
      next
    }
    
    results_df <- data.frame(
      Gene = gene_names2[keep],
      Correlation = as.numeric(r[keep]),
      PValue = as.numeric(pvals[keep]),
      stringsAsFactors = FALSE
    )
    rownames(results_df) <- results_df$Gene
    
    results_list[[GSE]] <- results_df
  }
  
  save(results_list, file = paste0("data/", given_gene, "_", Class, "_array_cor_list.rdata"))
  results_list
}
# step2函数 汇总相关--------------
create_gene_summary <- function(
                                correlation_results, 
                                filter_threshold = 0.1, 
                                confidence_threshold = 0.75,
                                given_gene) {
  combined_results <- data.frame(Gene = character(), 
                                 Correlation = numeric(), 
                                 PValue = numeric(), 
                                 GSE = character(), 
                                 Sample_Count = numeric(), 
                                 stringsAsFactors = FALSE)
  
  for (GSE in names(correlation_results)) {
    results_df <- correlation_results[[GSE]]
    sample_count <- ncol(tpm_samples_list[[GSE]])
    
    results_df$GSE <- GSE
    results_df$Sample_Count <- sample_count
    combined_results <- rbind(combined_results, results_df)
  }
  
  gene_summary <- combined_results %>%
    group_by(Gene) %>%
    summarize(Count = n(),                     
              Total_Sample_Count = sum(Sample_Count),  
              Correlation_PValue_Combo = paste0(GSE, ": ", 
                                                Correlation, ", ", 
                                                PValue, collapse = "; "),
              .groups = 'drop')
  
  total_studies <- length(correlation_results)
  included_studies <- names(correlation_results)
  
  total_samples <- sample_info_all_final %>%
    filter(GSE %in% included_studies) %>%
    group_by(GSE) %>%
    summarise(n = n()) %>%
    summarise(Total_Sample_Count = sum(n)) %>%
    pull(Total_Sample_Count)
  
  gene_summary <- gene_summary %>%
    mutate(Filter = ifelse(Count > total_studies * filter_threshold, "Keep", "Filter"),
           Weight = (Count / total_studies) * (Total_Sample_Count / total_samples))
  
  quantiles <- gene_summary %>%
    filter(Filter == "Keep") %>%
    summarise(Q3 = quantile(Weight, confidence_threshold)) %>%
    pull(Q3)
  
  gene_summary <- gene_summary %>%
    mutate(Confidence = ifelse(Filter == "Keep", 
                               ifelse(Weight > quantiles, "HC", "LC"), 
                               NA)) %>% 
    arrange(desc(Weight)) %>% 
    select(Gene,Count, Weight,everything())
  
  print(total_studies)
  print(total_samples)
  print(table(gene_summary$Filter))
  print(table(gene_summary$Confidence))
  save(gene_summary, file = paste0('data/',given_gene,'_icn_summary.rdata'))
  
  # save(gene_summary, file = 'data/gene_summary.rdata')
  return(gene_summary)
  
  
}

# array汇总-------------
array_create_gene_summary <- function(
                                      correlation_results, 
                                      filter_threshold = 0.1, 
                                      confidence_threshold = 0.5,
                                      given_gene) {
  combined_results <- data.frame(Gene = character(), 
                                 Correlation = numeric(), 
                                 PValue = numeric(), 
                                 GSE = character(), 
                                 Sample_Count = numeric(), 
                                 stringsAsFactors = FALSE)
  
  # 遍历每个GSE并组合结果
  for (GSE in names(correlation_results)) {
    results_df <- correlation_results[[GSE]]
    sample_count <- ncol(array_exp_list[[GSE]])
    
    results_df$GSE <- GSE
    results_df$Sample_Count <- sample_count
    combined_results <- rbind(combined_results, results_df)
  }
  
  # 汇总基因信息
  gene_summary <- combined_results %>%
    group_by(Gene) %>%
    summarize(Count = n(),                     
              Total_Sample_Count = sum(Sample_Count),  
              Correlation_PValue_Combo = paste0(GSE, ": ", 
                                                Correlation, ", ", 
                                                PValue, collapse = "; "),
              .groups = 'drop')
  
  # 计算研究总数和样本总数
  total_studies <- length(correlation_results)
  included_studies <- names(correlation_results)
  
  total_samples <- array_sample_info_all %>%
    filter(GSE %in% included_studies) %>%
    group_by(GSE) %>%
    summarise(n = n()) %>%
    summarise(Total_Sample_Count = sum(n)) %>%
    pull(Total_Sample_Count)
  
  # 计算基因的过滤标记和权重
  gene_summary <- gene_summary %>%
    mutate(Filter = ifelse(Count > total_studies * filter_threshold, "Keep", "Filter"),
           Weight = (Count / total_studies) * (Total_Sample_Count / total_samples))
  
  # 根据设置的 `confidence_threshold` 计算对应分位数
  quantile_value <- gene_summary %>%
    filter(Filter == "Keep") %>%
    summarise(Q = quantile(Weight, confidence_threshold)) %>%
    pull(Q)
  
  # 基于分位数设置 Confidence 列为 "HC" 或 "LC"
  gene_summary <- gene_summary %>%
    mutate(Confidence = ifelse(Filter == "Keep", 
                               ifelse(Weight > quantile_value, "HC", "LC"), 
                               NA)) %>% 
    arrange(desc(Weight)) %>% 
    select(Gene, Count, Weight, everything())
  
  # 输出统计信息
  print(total_studies)
  print(total_samples)
  print(table(gene_summary$Filter))
  print(table(gene_summary$Confidence))
  
  # 保存基因汇总结果
  # save(gene_summary, file = 'data/array_gene_summary.rdata')
  save(gene_summary, file = paste0('data/',given_gene,'_array_icn_summary.rdata'))
  
  return(gene_summary)
}

# step3 基因功能注释--------------------
perform_ora <- function(
                        gene_summary_result,
                        given_gene) {
  library(clusterProfiler)
  library(stringr)
  
  # 加载GMT文件
  load('data/gmt_lists.rdata')
  
  # 初始化结果列表
  ora_list <- list()
  geneset <- gene_summary_result$Gene[gene_summary_result$Confidence == 'HC']
  
  for (gmt in names(gmt_lists)) {
    em <- enricher(geneset,
                   TERM2GENE = gmt_lists[[gmt]],
                   pvalueCutoff = 0.1,
                   maxGSSize = 100000)
    
    if (!is.null(em)) {
      em_df <- data.frame(em) %>% filter(pvalue < 0.05)
    } else {
      em_df <- NULL
    }
    
    
    ora_list[[gmt]] <- em_df
  }
  save(ora_list, file = paste0('data/',given_gene,'_ora_list.rdata'))
  
  # save(ora_list, file = 'data/ora_list.rdata')
  return(ora_list)
  
}

# array 注释-----------
array_perform_ora <- function(
                              gene_summary_result,
                              given_gene) {
  library(clusterProfiler)
  library(stringr)
  
  # 加载GMT文件
  load('data/gmt_lists.rdata')
  
  # 初始化结果列表
  ora_list <- list()
  geneset <- gene_summary_result$Gene[gene_summary_result$Confidence == 'HC']
  
  for (gmt in names(gmt_lists)) {
    em <- enricher(geneset,
                   TERM2GENE = gmt_lists[[gmt]],
                   pvalueCutoff = 0.1,
                   maxGSSize = 100000)
    
    if (!is.null(em)) {
      em_df <- data.frame(em) %>% filter(pvalue < 0.05)
    } else {
      em_df <- NULL
    }
    
    ora_list[[gmt]] <- em_df
  }
  save(ora_list, file = paste0('data/',given_gene,'_array_ora_list.rdata'))
  
  # save(ora_list, file = 'data/array_ora_list.rdata')
  return(ora_list)
}


# gtex-------------
# 相关性分析--------------
compute_gtex_gene_correlation <- function(gene, 
                                     threshold = 2, 
                                     gtex_muscle_tpm, 
                                     cor_threshold = 0.5, 
                                     p_threshold = 0.05) {
  # Step 1: 检查阈值是否在可接受范围（0.5, 1, 2）
  if (!threshold %in% c(0.5, 1, 2)) {
    stop("Threshold must be one of the following values: 0.5, 1, 2.")
  }
  
  # Step 2: 检查目标基因是否存在于原始的gtex_muscle_tpm数据中
  if (!(gene %in% colnames(gtex_muscle_tpm))) {
    stop(paste("Gene", gene, "not found in the provided gtex_muscle_tpm data."))
  }
  
  # Step 3: 计算每个基因的均值
  gene_means <- colMeans(gtex_muscle_tpm)
  
  # Step 4: 过滤低表达基因，保留均值大于阈值的基因
  gtex_muscle_tpm_filtered <- gtex_muscle_tpm[, gene_means > threshold]
  
  # Step 5: 检查输入基因是否存在于过滤后的基因列表中
  if (!(gene %in% colnames(gtex_muscle_tpm_filtered))) {
    stop(paste("Gene", gene, "not found in the filtered dataset after applying the threshold."))
  }
  
  # Step 6: 提取目标基因表达数据
  gene_expression <- gtex_muscle_tpm_filtered[, gene]
  
  # Step 7: 初始化存储结果的数据框
  results <- data.frame(Gene = character(), r = numeric(), p = numeric(), stringsAsFactors = FALSE)
  
  # Step 8: 计算目标基因与其他基因的皮尔逊相关系数和显著性
  for (other_gene in colnames(gtex_muscle_tpm_filtered)) {
    if (other_gene != gene) {
      # 计算皮尔逊相关系数和p值
      cor_test <- try(cor.test(gene_expression, gtex_muscle_tpm_filtered[, other_gene], method = "pearson"), silent = TRUE)
      
      # 检查cor_test是否成功
      if (!inherits(cor_test, "try-error") && !is.null(cor_test$estimate) && !is.null(cor_test$p.value)) {
        # 将结果添加到数据框中
        results <- rbind(results, data.frame(Gene = other_gene, r = cor_test$estimate, p = cor_test$p.value))
      }
    }
  }
  
  # Step 9: 对结果进行过滤，基于相关性阈值和显著性阈值，并按r降序排序
  results_filtered <- results %>%
    filter(r > cor_threshold & p < p_threshold) %>%
    arrange(desc(r))  # 按r降序排序
  
  # Step 10: 将Gene列中的点替换为横线
  results_filtered$Gene <- gsub("\\.", "-", results_filtered$Gene)
  
  # Step 11: 设置行名为基因名
  rownames(results_filtered) <- results_filtered$Gene
  
  return(results_filtered)
}

# ora-----------
perform_ora_gtex <- function(cor_result_df) {
  # 加载必要的包
  library(clusterProfiler)
  library(dplyr)
  
  # 加载GMT文件
  load('data/gmt_lists.rdata')
  
  # 初始化结果列表
  ora_list <- list()
  
  # 提取相关性显著的基因集
  geneset <- cor_result_df$Gene
  
  # 针对每个GMT集执行 ORA 分析
  for (gmt in names(gmt_lists)) {
    em <- enricher(geneset,
                   TERM2GENE = gmt_lists[[gmt]],
                   pvalueCutoff = 0.1,
                   maxGSSize = 100000)
    
    # 筛选 pvalue < 0.05 的结果
    if (!is.null(em)) {
      em_df <- data.frame(em) %>% filter(pvalue < 0.05)
    } else {
      em_df <- NULL
    }
    
    
    # 将结果保存到列表中
    ora_list[[gmt]] <- em_df
  }
  
  # 将ORA结果保存为文件 ora_list_gtex.rdata
  save(ora_list, file = 'data/ora_list_gtex.rdata')
  
  return(ora_list)
}




