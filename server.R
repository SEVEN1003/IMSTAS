server <- function(input, output, session) {
  # load deg summary-----------
  # 在应用启动时自动加载 bulk_deg_summary 数据
  output$bulk_deg_summary_table <- DT::renderDataTable({
    DT::datatable(bulk_deg_summary[,-5],  filter = "top",  # 启用列筛选功能，筛选框显示在每列的上方
                  options = list(scrollX = TRUE, pageLength = 10))
  })
  
  # 下载 bulk_deg_summary 数据
  output$download_bulk_deg_summary <- downloadHandler(
    filename = function() {
      paste("bulk_deg_summary_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(bulk_deg_summary, file, row.names = FALSE)
    }
  )
  
  # 在应用启动时自动加载 array_deg_summary 数据
  output$array_deg_summary_table <- DT::renderDataTable({
    DT::datatable(array_deg_summary[,-5], filter = "top",
                  options = list(scrollX = TRUE, pageLength = 10))
  })
  
  # 下载 array_deg_summary 数据
  output$download_array_deg_summary <- downloadHandler(
    filename = function() {
      paste("array_deg_summary_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(array_deg_summary, file, row.names = FALSE)
    }
  )
  
  # 加载 GTEx 趋势分析摘要表格
  output$gtex_trend_summary_table <- DT::renderDataTable({
    DT::datatable(gtex_trend_summary,  filter = "top",
                  options = list(scrollX = TRUE, pageLength = 10))  # 启用列筛选功能
  })
  
  # single_class_tab ----------------
  observeEvent(input$run_trend_analysis, {
    
    # 检查是否输入了基因，并且基因是否在 bulk_gene 中
    if (input$gene_trend_input == "" || !input$gene_trend_input %in% bulk_gene) {
      showModal(modalDialog(
        title = "Input Error",
        "Please enter a valid gene from the bulk gene list.",
        easyClose = TRUE,
        footer = NULL
      ))
      return(NULL)  # 停止执行
    }
    
    # 提示用户分析可能花费 30 秒时间
    showNotification("Processing... This may take around 5 seconds to refresh the results.", 
                     type = "message", 
                     duration = 5)  # 不自动关闭
    
    
    # 获取用户输入的值
    Class <- input$class_trend_input
    given_gene <- input$gene_trend_input
    sort_order <- input$order_trend_input
    
    # 计算中位数差异
    median_diff_results <- calculate_group_median_difference_by_class(Class, 
                                                                      given_gene, 
                                                                      sort_order = sort_order)
    
    # 绘制排序条形图
    bar_diff_plot <- plot_bar_diff_median(median_diff_results, 
                                          given_gene, 
                                          Class,
                                          sort_order = sort_order)
    
    # 将条形图输出到 UI
    output$bar_trend_output <- renderPlot({
      print(bar_diff_plot)
    })
    
    # 添加PDF下载功能
    output$download_barplot_pdf <- downloadHandler(
      filename = function() {
        paste0(given_gene, "_barplot.pdf")
      },
      content = function(file) {
        ggsave(file, plot = bar_diff_plot, device = "pdf", 
               width = input$barplot_width, height = input$barplot_height)
      }
    )
    
    # 生成箱式图列表
    boxplot_list <- plot_boxplot_by_gse_list(median_diff_results, given_gene)
    
    # 初始化 GSE 选择器
    updateSelectInput(session, "gse_boxplot_selector", 
                      choices = names(boxplot_list), 
                      selected = names(boxplot_list)[1])
    
    # 显示第一个箱式图
    output$boxplot_trend_output <- renderPlot({
      print(boxplot_list[[input$gse_boxplot_selector]])
    })
    
    # 添加PDF下载功能
    output$download_boxplot_pdf <- downloadHandler(
      filename = function() {
        paste0(given_gene, "_boxplot.pdf")
      },
      content = function(file) {
        ggsave(file, plot = boxplot_list[[input$gse_boxplot_selector]], device = "pdf", 
               width = input$boxplot_width, height = input$boxplot_height)
      }
    )
    
    # 当前选择的 GSE
    current_gse <- reactiveVal(names(boxplot_list)[1])  # 初始化为第一个 GSE
    
    # 上一个 GSE 按钮
    observeEvent(input$prev_gse, {
      current_idx <- match(current_gse(), names(boxplot_list))
      if (current_idx > 1) {
        new_gse <- names(boxplot_list)[current_idx - 1]
        current_gse(new_gse)
        updateSelectInput(session, "gse_boxplot_selector", selected = new_gse)
      }
    })
    
    # 下一个 GSE 按钮
    observeEvent(input$next_gse, {
      current_idx <- match(current_gse(), names(boxplot_list))
      if (current_idx < length(names(boxplot_list))) {
        new_gse <- names(boxplot_list)[current_idx + 1]
        current_gse(new_gse)
        updateSelectInput(session, "gse_boxplot_selector", selected = new_gse)
      }
    })
    
    # 根据 GSE 选择器切换箱式图
    observeEvent(input$gse_boxplot_selector, {
      output$boxplot_trend_output <- renderPlot({
        print(boxplot_list[[input$gse_boxplot_selector]])
      })
      
      # 更新提示信息
      output$gse_keywords_output <- renderText({
        paste(
          "Keywords:", median_diff_results$Keywords[median_diff_results$GSE == input$gse_boxplot_selector], 
          "\nKeywords (ZH):", median_diff_results$Keywords_zh[median_diff_results$GSE == input$gse_boxplot_selector]
        )
      })
    })
    
    # 输出表格
    output$table_trend_output <- DT::renderDataTable({
      DT::datatable(median_diff_results, 
                    options = list(pageLength = 10, scrollX = TRUE))
    })
    
    # 表格下载功能
    output$download_table <- downloadHandler(
      filename = function() {
        paste0("median_diff_results_", given_gene, ".csv")
      },
      content = function(file) {
        write.csv(median_diff_results, file, row.names = FALSE)
      }
    )
    
  })
  
  # micro_single_class_tab ----------------
  observeEvent(input$run_micro_trend_analysis, {
    
    # 检查是否输入了基因，并且基因是否在 array_gene 中
    if (input$micro_gene_trend_input == "" || !input$micro_gene_trend_input %in% array_gene) {
      showModal(modalDialog(
        title = "Input Error",
        "Please enter a valid gene from the array gene list.",
        easyClose = TRUE,
        footer = NULL
      ))
      return(NULL)  # 停止执行
    }
    
    # 提示用户分析可能花费 30 秒时间
    showNotification("Processing... This may take around 5 seconds to refresh the results.", 
                     type = "message", 
                     duration = 5)  # 不自动关闭
    
    
    # 获取用户输入的值
    Class <- input$micro_class_trend_input
    given_gene <- input$micro_gene_trend_input
    sort_order <- input$micro_order_trend_input
    
    # 计算中位数差异
    median_diff_results <- calculate_array_group_median_difference_by_class(Class, 
                                                                            given_gene, 
                                                                            sort_order = sort_order)
    
    # 绘制排序条形图
    bar_diff_plot <- plot_array_bar_diff_median(median_diff_results, 
                                                given_gene, 
                                                Class,
                                                sort_order = sort_order)
    
    # 将条形图输出到 UI
    output$micro_bar_trend_output <- renderPlot({
      print(bar_diff_plot)
    })
    
    # 添加PDF下载功能
    output$micro_download_barplot_pdf <- downloadHandler(
      filename = function() {
        paste0(given_gene, "_micro_barplot.pdf")
      },
      content = function(file) {
        ggsave(file, plot = bar_diff_plot, device = "pdf", 
               width = input$micro_barplot_width, height = input$micro_barplot_height)
      }
    )
    
    # 生成箱式图列表
    boxplot_list <- plot_array_boxplot_by_gse_list(median_diff_results, given_gene)
    
    # 初始化 GSE 选择器
    updateSelectInput(session, "micro_gse_boxplot_selector", 
                      choices = names(boxplot_list), 
                      selected = names(boxplot_list)[1])
    
    # 显示第一个箱式图
    output$micro_boxplot_trend_output <- renderPlot({
      print(boxplot_list[[input$micro_gse_boxplot_selector]])
    })
    
    # 添加PDF下载功能
    output$micro_download_boxplot_pdf <- downloadHandler(
      filename = function() {
        paste0(given_gene, "_micro_boxplot.pdf")
      },
      content = function(file) {
        ggsave(file, plot = boxplot_list[[input$micro_gse_boxplot_selector]], device = "pdf", 
               width = input$micro_boxplot_width, height = input$micro_boxplot_height)
      }
    )
    
    # 当前选择的 GSE
    current_gse <- reactiveVal(names(boxplot_list)[1])  # 初始化为第一个 GSE
    
    # 上一个 GSE 按钮
    observeEvent(input$micro_prev_gse, {
      current_idx <- match(current_gse(), names(boxplot_list))
      if (current_idx > 1) {
        new_gse <- names(boxplot_list)[current_idx - 1]
        current_gse(new_gse)
        updateSelectInput(session, "micro_gse_boxplot_selector", selected = new_gse)
      }
    })
    
    # 下一个 GSE 按钮
    observeEvent(input$micro_next_gse, {
      current_idx <- match(current_gse(), names(boxplot_list))
      if (current_idx < length(names(boxplot_list))) {
        new_gse <- names(boxplot_list)[current_idx + 1]
        current_gse(new_gse)
        updateSelectInput(session, "micro_gse_boxplot_selector", selected = new_gse)
      }
    })
    
    # 根据 GSE 选择器切换箱式图
    observeEvent(input$micro_gse_boxplot_selector, {
      output$micro_boxplot_trend_output <- renderPlot({
        print(boxplot_list[[input$micro_gse_boxplot_selector]])
      })
      
      # 更新提示信息
      output$micro_gse_keywords_output <- renderText({
        paste(
          "Keywords:", median_diff_results$Keywords[median_diff_results$GSE == input$micro_gse_boxplot_selector], 
          "\nKeywords (ZH):", median_diff_results$Keywords_zh[median_diff_results$GSE == input$micro_gse_boxplot_selector]
        )
      })
    })
    
    # 输出表格
    output$micro_table_trend_output <- DT::renderDataTable({
      DT::datatable(median_diff_results, 
                    options = list(pageLength = 10, scrollX = TRUE))
    })
    
    # 表格下载功能
    output$micro_download_table <- downloadHandler(
      filename = function() {
        paste0("micro_median_diff_results_", given_gene, ".csv")
      },
      content = function(file) {
        write.csv(median_diff_results, file, row.names = FALSE)
      }
    )
    
  })
  
  
  
  
  # all_datasets_tab --------------
  observeEvent(input$run_all_datasets_analysis, {
    if (input$gene_cor_input == "" || !input$gene_cor_input %in% bulk_gene) {
      showModal(modalDialog(
        title = "Input Error",
        "Please enter a gene name before running the analysis.",
        easyClose = TRUE,
        footer = NULL
      ))
      return(NULL)  # 停止执行
    }
    
    given_gene <- input$gene_cor_input
    Class <- input$class_cor_input
    
    # 强制转换为数值型
    filter_threshold <- as.numeric(input$filter_threshold_cor_input)
    confidence_threshold <- as.numeric(input$confidence_threshold_cor_input)
    
    # 构造文件路径
    file_path <- paste0("data/", given_gene, '_', Class, "_cor_list.rdata")
    
    if (file.exists(file_path)) {
      # 提示使用已经计算好的结果
      showNotification("Using precomputed results...", type = "message", duration = 5)
      
      # 加载现有的计算结果
      load(file_path)
      correlation_results <- results_list  # 假设文件中有该对象
    } else {
      # 提示开始计算相关网络
      showNotification("Starting calculation of the correlation network...This may take approximately 10 minutes.", type = "warning", duration = 60)
      
      # 开始计算
      correlation_results <- calculate_gene_correlations(given_gene, Class = Class)
      
      # 计算完成后提示
      showNotification("Correlation network calculation completed.", type = "message", duration = 5)
    }
    
    gene_summary_result <- create_gene_summary(correlation_results,
                                               filter_threshold = filter_threshold,
                                               confidence_threshold = confidence_threshold)
    
    ora_results <- perform_ora(gene_summary_result)
    ora_barplot_list <- plot_enrich_bar_list(ora_results)
    
    total_studies <- as.numeric(length(correlation_results))
    gene_summary_plot_list <- plot_gene_summary(gene_summary_result,
                                                total_studies = total_studies,
                                                filter_threshold = filter_threshold)
    
    # 表格输出
    output$cor_table_output <- DT::renderDataTable({
      DT::datatable(gene_summary_result[,-5],
                    options = list(scrollX = TRUE, pageLength = 10)
      )
    })
    
    # 表格下载功能
    output$download_cor_table <- downloadHandler(
      filename = function() {
        paste("Correlation_Table_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(gene_summary_result, file, row.names = FALSE)
      }
    )
    
    # 更新图形输出，使用 ggplot2
    output$Histogram_Count_output <- renderPlot({
      gene_summary_plot_list$Histogram_Count
    })
    
    output$Barplot_Confidence_output <- renderPlot({
      gene_summary_plot_list$Barplot_Confidence
    })
    
    output$Top20_Lollipop_output <- renderPlot({
      gene_summary_plot_list$Top20_Lollipop
    })
    
    # 下载所有图为一个PDF
    output$download_plots_pdf <- downloadHandler(
      filename = function() {
        paste("Plots_", Sys.Date(), ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, height = input$pdf_height_plots, width = input$pdf_width_plots)
        print(gene_summary_plot_list$Histogram_Count)
        print(gene_summary_plot_list$Barplot_Confidence)
        print(gene_summary_plot_list$Top20_Lollipop)
        dev.off()
      }
    )
    
    # Pathway Outputs (BP, KEGG, Reactome)
    output$BP_output <- renderPlot({
      ora_barplot_list$BP
    })
    
    output$KEGG_LEGACY_output <- renderPlot({
      ora_barplot_list$KEGG_LEGACY
    })
    
    output$REACTOME_output <- renderPlot({
      ora_barplot_list$REACTOME
    })
    
    # BP Output 下载功能
    output$download_BP <- downloadHandler(
      filename = function() {
        paste("BP_output_", Sys.Date(), ".pdf", sep = "")
      },
      content = function(file) {
        ggsave(file, plot = ora_barplot_list$BP, device = "pdf", 
               height = input$pdf_height_BP, width = input$pdf_width_BP)
      }
    )
    
    # KEGG Legacy Output 下载功能
    output$download_KEGG <- downloadHandler(
      filename = function() {
        paste("KEGG_output_", Sys.Date(), ".pdf", sep = "")
      },
      content = function(file) {
        ggsave(file, plot = ora_barplot_list$KEGG_LEGACY, device = "pdf", 
               height = input$pdf_height_KEGG, width = input$pdf_width_KEGG)
      }
    )
    
    # Reactome Output 下载功能
    output$download_REACTOME <- downloadHandler(
      filename = function() {
        paste("Reactome_output_", Sys.Date(), ".pdf", sep = "")
      },
      content = function(file) {
        ggsave(file, plot = ora_barplot_list$REACTOME, device = "pdf", 
               height = input$pdf_height_REACTOME, width = input$pdf_width_REACTOME)
      }
    )
  })
  
  
  # micro_all_datasets_tab --------------
  observeEvent(input$micro_run_all_datasets_analysis, {
    if (input$micro_gene_cor_input == "" || !input$micro_gene_cor_input %in% array_gene) {
      showModal(modalDialog(
        title = "Input Error",
        "Please enter a gene name before running the analysis.",
        easyClose = TRUE,
        footer = NULL
      ))
      return(NULL)  # 停止执行
    }
    
    given_gene <- input$micro_gene_cor_input
    Class <- input$micro_class_cor_input
    
    # 强制转换为数值型
    filter_threshold <- as.numeric(input$micro_filter_threshold_cor_input)
    confidence_threshold <- as.numeric(input$micro_confidence_threshold_cor_input)
    
    # 构造文件路径
    file_path <- paste0("data/", given_gene, '_', Class, "_array_cor_list.rdata")
    
    if (file.exists(file_path)) {
      # 提示使用已经计算好的结果
      showNotification("Using precomputed results...", type = "message", duration = 5)
      
      # 加载现有的计算结果
      load(file_path)
      correlation_results <- results_list  # 假设文件中有该对象
    } else {
      # 提示开始计算相关网络
      showNotification("Starting calculation of the correlation network...This may take approximately 10 minutes.", type = "warning", duration = 60)
      
      # 开始计算
      correlation_results <- array_calculate_gene_correlations(given_gene, 
                                                               Class = Class)
      
      # 计算完成后提示
      showNotification("Correlation network calculation completed.", type = "message", duration = 5)
    }
    
    gene_summary_result <- array_create_gene_summary(correlation_results,
                                                     filter_threshold = filter_threshold,
                                                     confidence_threshold = confidence_threshold)
    
    ora_results <- array_perform_ora(gene_summary_result)
    ora_barplot_list <- array_plot_enrich_bar_list(ora_results)
    
    total_studies <- as.numeric(length(correlation_results))
    gene_summary_plot_list <- array_plot_gene_summary(gene_summary_result,
                                                      total_studies = total_studies,
                                                      filter_threshold = filter_threshold)
    
    # 表格输出
    output$micro_cor_table_output <- DT::renderDataTable({
      DT::datatable(gene_summary_result[,-5],
                    options = list(scrollX = TRUE, pageLength = 10)
      )
    })
    
    # 表格下载功能
    output$micro_download_cor_table <- downloadHandler(
      filename = function() {
        paste("Correlation_Table_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(gene_summary_result, file, row.names = FALSE)
      }
    )
    
    # 更新图形输出，使用 ggplot2
    output$micro_Histogram_Count_output <- renderPlot({
      gene_summary_plot_list$Histogram_Count
    })
    
    output$micro_Barplot_Confidence_output <- renderPlot({
      gene_summary_plot_list$Barplot_Confidence
    })
    
    output$micro_Top20_Lollipop_output <- renderPlot({
      gene_summary_plot_list$Top20_Lollipop
    })
    
    # 下载所有图为一个PDF
    output$micro_download_plots_pdf <- downloadHandler(
      filename = function() {
        paste("Micro_Plots_", Sys.Date(), ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, height = input$micro_pdf_height_plots, width = input$micro_pdf_width_plots)
        print(gene_summary_plot_list$Histogram_Count)
        print(gene_summary_plot_list$Barplot_Confidence)
        print(gene_summary_plot_list$Top20_Lollipop)
        dev.off()
      }
    )
    
    # Pathway Outputs (BP, KEGG, Reactome)
    output$micro_BP_output <- renderPlot({
      ora_barplot_list$BP
    })
    
    output$micro_KEGG_LEGACY_output <- renderPlot({
      ora_barplot_list$KEGG_LEGACY
    })
    
    output$micro_REACTOME_output <- renderPlot({
      ora_barplot_list$REACTOME
    })
    
    # BP Output 下载功能
    output$micro_download_BP <- downloadHandler(
      filename = function() {
        paste("Micro_BP_output_", Sys.Date(), ".pdf", sep = "")
      },
      content = function(file) {
        ggsave(file, plot = ora_barplot_list$BP, device = "pdf", 
               height = input$micro_pdf_height_BP, width = input$micro_pdf_width_BP)
      }
    )
    
    # KEGG Legacy Output 下载功能
    output$micro_download_KEGG <- downloadHandler(
      filename = function() {
        paste("Micro_KEGG_output_", Sys.Date(), ".pdf", sep = "")
      },
      content = function(file) {
        ggsave(file, plot = ora_barplot_list$KEGG_LEGACY, device = "pdf", 
               height = input$micro_pdf_height_KEGG, width = input$micro_pdf_width_KEGG)
      }
    )
    
    # Reactome Output 下载功能
    output$micro_download_REACTOME <- downloadHandler(
      filename = function() {
        paste("Micro_Reactome_output_", Sys.Date(), ".pdf", sep = "")
      },
      content = function(file) {
        ggsave(file, plot = ora_barplot_list$REACTOME, device = "pdf", 
               height = input$micro_pdf_height_REACTOME, width = input$micro_pdf_width_REACTOME)
      }
    )
  })
  
  
  
  
  # atlas_tab ----------------
  observeEvent(input$run_altas_analysis, {
    
    # 检查 Gene 是否在 altas_gene 中
    if (input$gene_input == "" || !input$gene_input %in% altas_gene) {
      showModal(modalDialog(
        title = "Input Error",
        "Please enter a valid gene from the atlas.",
        easyClose = TRUE,
        footer = NULL
      ))
      return(NULL)  # 停止执行
    }
    
    # 提示用户分析可能花费 30 秒时间
    showNotification("Processing... This may take around 30 seconds to refresh the results.", 
                     type = "message", 
                     duration = 30)  # 不自动关闭
    
    # 获取输入值
    Gene <- input$gene_input
    V1 <- input$v1_input
    V2 <- if (input$v2_input == "All") NULL else input$v2_input
    
    # 生成 Atlas 图像列表
    altas_plots <- plot_altas(Gene, V1, V2)
    
    # UMAP Annotation 输出
    output$atlas_UMAP_annotation_output <- renderPlot({
      print(altas_plots[[1]][[1]])  # 显示 UMAP Annotation 图
    })
    
    # UMAP Feature 输出
    output$atlas_UMAP_feature_output <- renderPlot({
      print(altas_plots[[1]][[3]])  # 显示 UMAP Feature 图
    })
    
    # CellStat 输出
    output$atlas_CellStat_output <- renderPlot({
      print(altas_plots[[1]][[2]])  # 显示 CellStat 图
    })
    
    # FeatureStat 输出
    output$atlas_FeatureStat_output <- renderPlot({
      print(altas_plots[[1]][[4]])  # 新增 FeatureStat 图输出
    })
    
    # Data_name Feature 输出
    data_names <- unique(muscle_core_altas$Data_name)
    current_data_name <- reactiveVal(data_names[1])  # 初始化为第一个 Data_name
    
    output$data_name_feature_plot_output <- renderPlot({
      print(altas_plots[[2]][[current_data_name()]])  # 根据 current_data_name 显示相应的 Data_name 图
    })
    
    # 切换到上一个 Data_name
    observeEvent(input$prev_data_name, {
      current_idx <- match(current_data_name(), data_names)
      if (current_idx > 1) {
        current_data_name(data_names[current_idx - 1])
      }
    })
    
    # 切换到下一个 Data_name
    observeEvent(input$next_data_name, {
      current_idx <- match(current_data_name(), data_names)
      if (current_idx < length(data_names)) {
        current_data_name(data_names[current_idx + 1])
      }
    })
    
    # 下载功能实现
    output$download_annotation_plot <- downloadHandler(
      filename = function() {
        paste("UMAP_Annotation_", Sys.Date(), ".pdf", sep = "")
      },
      content = function(file) {
        ggsave(file, plot = altas_plots[[1]][[1]], device = "pdf", height = input$pdf_height_annotation, width = input$pdf_width_annotation)
      }
    )
    
    output$download_feature_plot <- downloadHandler(
      filename = function() {
        paste("UMAP_Feature_", Sys.Date(), ".pdf", sep = "")
      },
      content = function(file) {
        ggsave(file, plot = altas_plots[[1]][[3]], device = "pdf", height = input$pdf_height_feature, width = input$pdf_width_feature)
      }
    )
    
    output$download_CellStat_plot <- downloadHandler(
      filename = function() {
        paste("CellStat_Plot_", Sys.Date(), ".pdf", sep = "")
      },
      content = function(file) {
        ggsave(file, plot = altas_plots[[1]][[2]], device = "pdf", height = input$pdf_height_CellStat, width = input$pdf_width_CellStat)
      }
    )
    
    output$download_FeatureStat_plot <- downloadHandler(
      filename = function() {
        paste("FeatureStat_Plot_", Sys.Date(), ".pdf", sep = "")
      },
      content = function(file) {
        ggsave(file, plot = altas_plots[[1]][[4]], device = "pdf", height = input$pdf_height_FeatureStat, width = input$pdf_width_FeatureStat)
      }
    )
    
    output$download_data_name_plot <- downloadHandler(
      filename = function() {
        paste("Data_name_Plot_", Sys.Date(), ".pdf", sep = "")
      },
      content = function(file) {
        ggsave(file, plot = altas_plots[[2]][[current_data_name()]], device = "pdf", height = input$pdf_height_data_name, width = input$pdf_width_data_name)
      }
    )
  })
  
  
  
  # gtex tab -------------
  observeEvent(input$run_gtex_analysis, {
    if (input$gtex_gene_input == "" || !input$gtex_gene_input %in% gtex_gene) {
      showModal(modalDialog(
        title = "Input Error",
        "Please enter a valid gene from the GTEx dataset.",
        easyClose = TRUE,
        footer = NULL
      ))
      return(NULL)  # 停止执行
    }
    
    # 提示用户分析可能需要时间
    showNotification("Processing... This may take around 30 seconds to refresh the results.", 
                     type = "message", duration = 30)
    
    # 获取输入值
    gene <- input$gtex_gene_input
    threshold <- as.numeric(input$gtex_threshold_input)
    cor_threshold <- as.numeric(input$gtex_cor_threshold_input)
    
    # 执行趋势分析
    result <- gtex_trend_analysis_with_plots(gene = gene, 
                                        gtex_muscle_info = gtex_muscle_info, 
                                        gtex_muscle_tpm = gtex_muscle_tpm, 
                                        age_group = c('40-49', '50-59', '60-69'), 
                                        center = "median", seed = 2024)
    
    # 更新趋势分析图
    output$gtex_male_trend_plot <- renderPlot({
      result$male_plot
    })
    
    output$gtex_female_trend_plot <- renderPlot({
      result$female_plot
    })
    
    # 更新趋势检验结果，提取并格式化为文本
    output$gtex_male_trend_result <- renderText({
      male_trend <- result$male_trend_result
      
      # 提取各个元素
      statistic <- male_trend$statistic
      p_value <- male_trend$p.value
      alternative <- male_trend$alternative
      method <- male_trend$method
      
      # 格式化输出文本
      result_text <- paste(
        "Trend Analysis Result:",
        paste("Statistic (JT):", round(statistic, 2)),
        paste("P-value:", format.pval(p_value, digits = 3)),
        paste("Alternative Hypothesis:", alternative),
        paste("Method:", method),
        sep = "\n"
      )
      
      return(result_text)
    })
    
    # 更新女性趋势检验结果，同样的处理方式
    output$gtex_female_trend_result <- renderText({
      female_trend <- result$female_trend_result
      
      # 提取各个元素
      statistic <- female_trend$statistic
      p_value <- female_trend$p.value
      alternative <- female_trend$alternative
      method <- female_trend$method
      
      # 格式化输出文本
      result_text <- paste(
        "Trend Analysis Result:",
        paste("Statistic (JT):", round(statistic, 2)),
        paste("P-value:", format.pval(p_value, digits = 3)),
        paste("Alternative Hypothesis:", alternative),
        paste("Method:", method),
        sep = "\n"
      )
      
      return(result_text)
    })
    
    
    # 下载趋势分析图
    output$download_gtex_male_trend <- downloadHandler(
      filename = function() {
        paste("Male_Trend_Plot_", Sys.Date(), ".pdf", sep = "")
      },
      content = function(file) {
        ggsave(file, plot = result$male_plot, device = "pdf", 
               height = input$gtex_pdf_height_male, width = input$gtex_pdf_width_male)
      }
    )
    
    output$download_gtex_female_trend <- downloadHandler(
      filename = function() {
        paste("Female_Trend_Plot_", Sys.Date(), ".pdf", sep = "")
      },
      content = function(file) {
        ggsave(file, plot = result$female_plot, device = "pdf", 
               height = input$gtex_pdf_height_female, width = input$gtex_pdf_width_female)
      }
    )
    
    # 相关分析结果
    result_df <- compute_gtex_gene_correlation(gene = gene, 
                                          threshold = threshold, 
                                          gtex_muscle_tpm = gtex_muscle_tpm, 
                                          cor_threshold = cor_threshold, 
                                          p_threshold = 0.05)
    
    # 相关分析图
    top20_lollipop_plot <- plot_gtex_top20_lollipop(result_df)
    
    output$gtex_top20_lollipop_plot <- renderPlot({
      top20_lollipop_plot
    })
    
    # 下载相关分析图
    output$download_gtex_top20_lollipop <- downloadHandler(
      filename = function() {
        paste("Top20_Lollipop_Plot_", Sys.Date(), ".pdf", sep = "")
      },
      content = function(file) {
        ggsave(file, plot = top20_lollipop_plot, device = "pdf", 
               height = input$gtex_pdf_height_lollipop, width = input$gtex_pdf_width_lollipop)
      }
    )
    
    # 相关分析表
    output$gtex_correlation_table <- DT::renderDataTable({
      DT::datatable(result_df, options = list(pageLength = 5, scrollX = TRUE))
    })
    
    # 下载相关分析表
    output$download_gtex_correlation_table <- downloadHandler(
      filename = function() {
        paste("Correlation_Table_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(result_df, file, row.names = FALSE)
      }
    )
    
    # 富集分析结果
    ora_results <- perform_ora_gtex(cor_result_df = result_df)
    ora_barplot_list <- plot_gtex_enrich_bar_list(ora_results)
    
    output$gtex_BP_output <- renderPlot({
      ora_barplot_list$BP
    })
    
    output$gtex_KEGG_LEGACY_output <- renderPlot({
      ora_barplot_list$KEGG_LEGACY
    })
    
    output$gtex_REACTOME_output <- renderPlot({
      ora_barplot_list$REACTOME
    })
    
    # 下载 BP, KEGG, Reactome 的图
    output$gtex_download_BP <- downloadHandler(
      filename = function() {
        paste("BP_output_", Sys.Date(), ".pdf", sep = "")
      },
      content = function(file) {
        ggsave(file, plot = ora_barplot_list$BP, device = "pdf", 
               height = input$gtex_pdf_height_BP, width = input$gtex_pdf_width_BP)
      }
    )
    
    output$gtex_download_KEGG <- downloadHandler(
      filename = function() {
        paste("KEGG_output_", Sys.Date(), ".pdf", sep = "")
      },
      content = function(file) {
        ggsave(file, plot = ora_barplot_list$KEGG_LEGACY, device = "pdf", 
               height = input$gtex_pdf_height_KEGG, width = input$gtex_pdf_width_KEGG)
      }
    )
    
    output$gtex_download_REACTOME <- downloadHandler(
      filename = function() {
        paste("Reactome_output_", Sys.Date(), ".pdf", sep = "")
      },
      content = function(file) {
        ggsave(file, plot = ora_barplot_list$REACTOME, device = "pdf", 
               height = input$gtex_pdf_height_REACTOME, width = input$gtex_pdf_width_REACTOME)
      }
    )
  })
  
  
  
}
