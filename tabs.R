# introduction_tab------------------
introduction_tab = bs4TabItem(
  tabName = "introduction",
  fluidRow(
    column(
      12,
      tags$img(
        src = "图片1.png",
        height = "400px",
        style = "display: block; margin: auto;"
      )
    ),
    column(
      12,
      h3("IMSTAS Introduction"),
      p(
        "IMSTAS的设计目标在于整合并检索预计算分析结果，从而快速筛选骨骼肌稳态维持及萎缩相关调控因子，并系统解析其在多种生理及疾病状态下的表达特征、ICN HC网络，以及其所参与的关键生物学功能与信号通路。依托上述模块化且连贯的一体化分析流程，IMSTAS可完整支持从关键调控因子筛选、年龄相关基因表达趋势分析、跨状态DEA查询和可视化、单细胞/单核水平组间DEA、基因-基因及变量相关性分析，到ICN构建与功能富集分析的全链条研究，为系统性挖掘骨骼肌稳态维持与萎缩驱动因子及其下游调控网络提供了高效、统一的分析框架。"
      ),
      p(
        "IMSTAS is designed to integrate and retrieve precomputed analytical results to rapidly identify regulatory factors involved in skeletal muscle homeostasis maintenance and atrophy. It systematically characterizes their expression patterns across diverse physiological and pathological conditions, as well as their associated ICN/HC networks, biological functions, and signaling pathways. Through a modular and coherent analytical framework, IMSTAS supports the full workflow from key regulator screening, age-related gene expression trend analysis, cross-condition DEA querying and visualization, single-cell/single-nucleus group-wise DEA, gene–gene and variable correlation analysis, to ICN construction and functional enrichment analysis. This integrated platform provides an efficient and unified framework for systematically uncovering drivers of skeletal muscle homeostasis and atrophy and their downstream regulatory networks."
      )
    )
  )
)




# single_class_tab UI --------------
single_class_tab <- bs4TabItem(
  tabName = "single_class",
  
  # 第一行：输入组件
  bs4Card(
    title = "Analysis Inputs",
    status = "primary",
    solidHeader = TRUE,
    width = 12,            # 占满页面宽度
    collapsible = TRUE,     # 允许收缩
    collapsed = FALSE,      # 默认展开
    fluidRow(
      column(3,
             selectInput("class_trend_input", "Select Class", 
                         choices = c("Drug", "Exercise", "Aging", "Others", "Metabolism"))
      ),
      column(3, 
             textInput("gene_trend_input", "Enter Gene Name", placeholder = "e.g., TP53")
      ),
      column(3, 
             selectInput("order_trend_input", "Select Order", 
                         choices = c("desc" = "desc", "asc" = "asc"), selected = "desc")
      )
    ),
    fluidRow(
      column(12,            # 将开始分析按钮放在第二行，占满整行宽度
             actionButton("run_trend_analysis", "Start Trend Analysis", class = "btn-primary")
      )
    )
  ),
  
  # 第二行：加载并显示 bulk_deg_summary 数据框
  bs4Card(
    title = "DEG Summary",
    status = "info",
    solidHeader = TRUE,
    width = 12,
    collapsible = TRUE,   
    collapsed = FALSE,    
    fluidRow(
      column(12,
             DT::dataTableOutput("bulk_deg_summary_table"),  # 显示数据框
             downloadButton("download_bulk_deg_summary", "Download Bulk DEG Summary as CSV")  # 下载按钮
      )
    )
  ),
  
  # 第二行：图形输出
  bs4Card(
    title = "DEA Table & Plots",
    status = "info",
    solidHeader = TRUE,
    width = NULL,
    collapsible = TRUE,  # 允许收缩
    collapsed = FALSE,    # 默认展开
    fluidRow(
      column(6, 
             plotOutput("bar_trend_output"),  # 左侧排序条图输出
             
             # 增加空白行
             # br(), br(), br(), # 用于对齐右侧的文本输出和选择框
             
             downloadButton("download_barplot_pdf", "Download Barplot as PDF"),  # 下载按钮
             numericInput("barplot_width", "Width (inches)", 8),  # 用户输入宽度
             numericInput("barplot_height", "Height (inches)", 6)  # 用户输入高度
      ),
      column(6, 
             selectInput("gse_boxplot_selector", "Select", choices = NULL),  # GSE下拉选择
             
             actionButton("prev_gse", "Previous"),  # 上一个GSE按钮
             actionButton("next_gse", "Next"),  # 下一个GSE按钮
             
             plotOutput("boxplot_trend_output"),  # 右侧箱式图输出
             
             textOutput("gse_keywords_output"),  # 显示GSE的提示信息
             
             
             downloadButton("download_boxplot_pdf", "Download Boxplot as PDF"),  # 下载按钮
             numericInput("boxplot_width", "Width (inches)", 8),  # 用户输入宽度
             numericInput("boxplot_height", "Height (inches)", 6)  # 用户输入高度
             
      )
    )
  ),
  
  
  # 第三行：表格输出
  bs4Card(
    title = "Correlation Analysis Plot",
    status = "success",
    solidHeader = TRUE,
    width = NULL,
    collapsible = TRUE,  # 允许收缩
    collapsed = FALSE,    # 默认展开
    fluidRow(
      column(12, 
             DT::dataTableOutput("table_trend_output"),  # 表格输出
             downloadButton("download_table", "Download Table as CSV")  # 下载表格按钮
      )
    )
  )
  
)



# all_datasets_tab------------
all_datasets_tab = bs4TabItem(
  tabName = "all_datasets",
  
  # 第一行：输入组件
  bs4Card(
    title = "Analysis Inputs",
    status = "primary",
    solidHeader = TRUE,
    width = 12,            # 占满页面宽度
    collapsible = TRUE,     # 允许收缩
    collapsed = FALSE,      # 默认展开
    fluidRow(
      column(3,             
             textInput("gene_cor_input", "Enter Gene Name:", placeholder = "e.g., NFE2L1")
      ),
      column(3,
             selectInput("class_cor_input", "Select Class:", 
                         choices = c('All', "Drug", "Exercise", "Aging", "Others", "Metabolism"),
                         selected = 'All')
      ),
      column(3,
             selectInput("filter_threshold_cor_input", "Select Filter Threshold:", 
                         choices = c(0.1, 0.2),
                         selected = 0.1)
      ),
      column(3,
             selectInput("confidence_threshold_cor_input", "Select Confidence Threshold:", 
                         choices = c(0.5, 0.75, 0.9,0.95),
                         selected = 0.9)
      )
    ),
    fluidRow(
      column(12,            
             actionButton("run_all_datasets_analysis", "Run Analysis", class = "btn-primary")
      )
    )
  ),
  
  # 第二行：图片输出
  bs4Card(
    title = "Plots",
    status = "info",
    solidHeader = TRUE,
    width = NULL,
    collapsible = TRUE,  
    collapsed = FALSE,    
    fluidRow(
      column(4,
             plotOutput("Histogram_Count_output")  
      ),
      column(4,
             plotOutput("Barplot_Confidence_output")  
      ),
      column(4,
             plotOutput("Top20_Lollipop_output")  
      )
    ),
    fluidRow(
      column(12, 
             numericInput("pdf_height_plots", "PDF Height", value = 8, min = 1),  # 设置PDF高度
             numericInput("pdf_width_plots", "PDF Width", value = 10, min = 1),   # 设置PDF宽度
             downloadButton("download_plots_pdf", "Download All Plots as PDF")    # 下载按钮
      )
    )
  ),
  
  # 第三行：Pathway Outputs
  bs4Card(
    title = "Pathway Outputs",
    status = "warning",
    solidHeader = TRUE,
    width = NULL,
    collapsible = TRUE,  
    collapsed = FALSE,   
    fluidRow(
      column(12,
             tabsetPanel(
               id = "pathway_tabs",
               tabPanel("BP Output", 
                        plotOutput("BP_output"),  
                        numericInput("pdf_height_BP", "PDF Height", value = 6, min = 1),
                        numericInput("pdf_width_BP", "PDF Width", value = 6, min = 1),
                        downloadButton("download_BP", "Download BP as PDF")  
               ),
               tabPanel("KEGG Legacy Output", 
                        plotOutput("KEGG_LEGACY_output"),  
                        numericInput("pdf_height_KEGG", "PDF Height", value = 6, min = 1),
                        numericInput("pdf_width_KEGG", "PDF Width", value = 6, min = 1),
                        downloadButton("download_KEGG", "Download KEGG as PDF")  
               ),
               tabPanel("Reactome Output", 
                        plotOutput("REACTOME_output"),  
                        numericInput("pdf_height_REACTOME", "PDF Height", value = 6, min = 1),
                        numericInput("pdf_width_REACTOME", "PDF Width", value = 6, min = 1),
                        downloadButton("download_REACTOME", "Download Reactome as PDF")  
               )
             )
      )
    )
  ),
  
  # 第四行：表格输出
  bs4Card(
    title = "Correlation Table",
    status = "success",
    solidHeader = TRUE,
    width = NULL,
    collapsible = TRUE,  
    collapsed = FALSE,   
    fluidRow(
      column(12,
             DT::dataTableOutput("cor_table_output"),  # 表格输出
             downloadButton("download_cor_table", "Download Table as CSV")  # 表格下载按钮
      )
    )
  )
)









# micro_single_class_tab-----------------
micro_single_class_tab = bs4TabItem(
  tabName = "micro_single_class",
  
  # 第一行：输入组件
  bs4Card(
    title = "Analysis Inputs",
    status = "primary",
    solidHeader = TRUE,
    width = 12,            # 占满页面宽度
    collapsible = TRUE,     # 允许收缩
    collapsed = FALSE,      # 默认展开
    fluidRow(
      column(3,
             selectInput("micro_class_trend_input", "Select Class", 
                         choices = c("Drug", "Exercise", "Aging", "Others", "Metabolism"))
      ),
      column(3, 
             textInput("micro_gene_trend_input", "Enter Gene Name", placeholder = "e.g., TP53")
      ),
      column(3, 
             selectInput("micro_order_trend_input", "Select Order", 
                         choices = c("desc" = "desc", "asc" = "asc"), selected = "desc")
      )
    ),
    fluidRow(
      column(12,            
             actionButton("run_micro_trend_analysis", "Start Microarray Trend Analysis", class = "btn-primary")
      )
    )
  ),
  
  # 第二行：加载并显示 array_deg_summary 数据框
  bs4Card(
    title = "DEG Summary",
    status = "info",
    solidHeader = TRUE,
    width = 12,
    collapsible = TRUE,   
    collapsed = FALSE,    
    fluidRow(
      column(12,
             DT::dataTableOutput("array_deg_summary_table"),  # 显示数据框
             downloadButton("download_array_deg_summary", "Download Microarray DEG Summary as CSV")  # 下载按钮
      )
    )
  ),
  
  # 第二行：图形输出
  bs4Card(
    title = "DEA Table & Plots",
    status = "info",
    solidHeader = TRUE,
    width = NULL,
    collapsible = TRUE, 
    collapsed = FALSE,    
    fluidRow(
      column(6, 
             plotOutput("micro_bar_trend_output"),
             downloadButton("micro_download_barplot_pdf", "Download Barplot as PDF"),
             numericInput("micro_barplot_width", "Width (inches)", 8),
             numericInput("micro_barplot_height", "Height (inches)", 6)
      ),
      column(6, 
             selectInput("micro_gse_boxplot_selector", "Select", choices = NULL),
             actionButton("micro_prev_gse", "Previous"),
             actionButton("micro_next_gse", "Next"),
             plotOutput("micro_boxplot_trend_output"),
             textOutput("micro_gse_keywords_output"),
             downloadButton("micro_download_boxplot_pdf", "Download Boxplot as PDF"),
             numericInput("micro_boxplot_width", "Width (inches)", 8),
             numericInput("micro_boxplot_height", "Height (inches)", 6)
      )
    )
  ),
  
  # 第三行：表格输出
  bs4Card(
    title = "Correlation Analysis Plot",
    status = "success",
    solidHeader = TRUE,
    width = NULL,
    collapsible = TRUE,  
    collapsed = FALSE,   
    fluidRow(
      column(12, 
             DT::dataTableOutput("micro_table_trend_output"),
             downloadButton("micro_download_table", "Download Table as CSV")
      )
    )
  )
)


# micro_all_datasets_tab---------------
micro_all_datasets_tab = bs4TabItem(
  tabName = "micro_all_datasets",
  
  # 第一行：输入组件
  bs4Card(
    title = "Analysis Inputs",
    status = "primary",
    solidHeader = TRUE,
    width = 12,            # 占满页面宽度
    collapsible = TRUE,     # 允许收缩
    collapsed = FALSE,      # 默认展开
    fluidRow(
      column(3,             
             textInput("micro_gene_cor_input", "Enter Gene Name:", placeholder = "e.g., NFE2L1")
      ),
      column(3,
             selectInput("micro_class_cor_input", "Select Class:", 
                         choices = c('All', "Drug", "Exercise", "Aging", "Others", "Metabolism"),
                         selected = 'All')
      ),
      column(3,
             selectInput("micro_filter_threshold_cor_input", "Select Filter Threshold:", 
                         choices = c(0.1, 0.2),
                         selected = 0.1)
      ),
      column(3,
             selectInput("micro_confidence_threshold_cor_input", "Select Confidence Threshold:", 
                         choices = c(0.5, 0.75, 0.9,0.95),
                         selected = 0.5)
      )
    ),
    fluidRow(
      column(12,            
             actionButton("micro_run_all_datasets_analysis", "Run Analysis", class = "btn-primary")
      )
    )
  ),
  
  # 第二行：图片输出
  bs4Card(
    title = "Plots",
    status = "info",
    solidHeader = TRUE,
    width = NULL,
    collapsible = TRUE,  
    collapsed = FALSE,    
    fluidRow(
      column(4,
             plotOutput("micro_Histogram_Count_output")  
      ),
      column(4,
             plotOutput("micro_Barplot_Confidence_output")  
      ),
      column(4,
             plotOutput("micro_Top20_Lollipop_output")  
      )
    ),
    fluidRow(
      column(12, 
             numericInput("micro_pdf_height_plots", "PDF Height", value = 8, min = 1),  # 设置PDF高度
             numericInput("micro_pdf_width_plots", "PDF Width", value = 10, min = 1),   # 设置PDF宽度
             downloadButton("micro_download_plots_pdf", "Download All Plots as PDF")    # 下载按钮
      )
    )
  ),
  
  # 第三行：Pathway Outputs
  bs4Card(
    title = "Pathway Outputs",
    status = "warning",
    solidHeader = TRUE,
    width = NULL,
    collapsible = TRUE,  
    collapsed = FALSE,   
    fluidRow(
      column(12,
             tabsetPanel(
               id = "micro_pathway_tabs",
               tabPanel("BP Output", 
                        plotOutput("micro_BP_output"),  
                        numericInput("micro_pdf_height_BP", "PDF Height", value = 6, min = 1),
                        numericInput("micro_pdf_width_BP", "PDF Width", value = 6, min = 1),
                        downloadButton("micro_download_BP", "Download BP as PDF")  
               ),
               tabPanel("KEGG Legacy Output", 
                        plotOutput("micro_KEGG_LEGACY_output"),  
                        numericInput("micro_pdf_height_KEGG", "PDF Height", value = 6, min = 1),
                        numericInput("micro_pdf_width_KEGG", "PDF Width", value = 6, min = 1),
                        downloadButton("micro_download_KEGG", "Download KEGG as PDF")  
               ),
               tabPanel("Reactome Output", 
                        plotOutput("micro_REACTOME_output"),  
                        numericInput("micro_pdf_height_REACTOME", "PDF Height", value = 6, min = 1),
                        numericInput("micro_pdf_width_REACTOME", "PDF Width", value = 6, min = 1),
                        downloadButton("micro_download_REACTOME", "Download Reactome as PDF")  
               )
             )
      )
    )
  ),
  
  # 第四行：表格输出
  bs4Card(
    title = "Correlation Table",
    status = "success",
    solidHeader = TRUE,
    width = NULL,
    collapsible = TRUE,  
    collapsed = FALSE,   
    fluidRow(
      column(12,
             DT::dataTableOutput("micro_cor_table_output"),  # 表格输出
             downloadButton("micro_download_cor_table", "Download Table as CSV")  # 表格下载按钮
      )
    )
  )
)







# GTEx-SM tab-----------------
gtex_sm_tab = bs4TabItem(
  tabName = "gtex_sm",
  
  # 第一行：输入组件
  bs4Card(
    title = "Analysis Inputs",
    status = "primary",
    solidHeader = TRUE,
    width = 12,
    collapsible = TRUE,
    collapsed = FALSE,
    fluidRow(
      column(4, 
             textInput("gtex_gene_input", "Enter Gene Name:", placeholder = "e.g., NFE2L1")
      ),
      column(4, 
             selectInput("gtex_threshold_input", "Select Expression Threshold:", 
                         choices = c(2, 1, 0.5), selected = 2)
      ),
      column(4, 
             selectInput("gtex_cor_threshold_input", "Select Correlation Threshold:", 
                         choices = c(0.5, 0.75, 0.9), selected = 0.5)
      )
    ),
    fluidRow(
      column(12, 
             actionButton("run_gtex_analysis", "Run Analysis", class = "btn-primary")
      )
    )
  ),
  
  # 第二行：显示 GTEx 趋势分析摘要表格
  bs4Card(
    title = "GTEx Trend Summary",
    status = "success",
    solidHeader = TRUE,
    width = 12,
    collapsible = TRUE,
    collapsed = FALSE,
    fluidRow(
      column(12,
             DT::dataTableOutput("gtex_trend_summary_table")  # 表格输出
      )
    )
  ),
  
  # 第二行：趋势分析结果
  bs4Card(
    title = "Trend Analysis Results",
    status = "info",
    solidHeader = TRUE,
    width = NULL,
    collapsible = TRUE,
    collapsed = FALSE,
    fluidRow(
      column(6, 
             plotOutput("gtex_male_trend_plot")),   # 男性趋势图
      column(6, 
             plotOutput("gtex_female_trend_plot"))  # 女性趋势图
    ),
    fluidRow(
      column(6, 
             textOutput("gtex_male_trend_result")),  # 男性趋势检验结果
      column(6, 
             textOutput("gtex_female_trend_result")) # 女性趋势检验结果
    ),
    fluidRow(
      column(6, 
             numericInput("gtex_pdf_height_male", "PDF Height", value = 6, min = 1),
             numericInput("gtex_pdf_width_male", "PDF Width", value = 8, min = 1),
             downloadButton("download_gtex_male_trend", "Download Male Trend Plot as PDF")
      ),
      column(6, 
             numericInput("gtex_pdf_height_female", "PDF Height", value = 6, min = 1),
             numericInput("gtex_pdf_width_female", "PDF Width", value = 8, min = 1),
             downloadButton("download_gtex_female_trend", "Download Female Trend Plot as PDF")
      )
    )
  ),
  
  # 第三行：相关分析结果
  bs4Card(
    title = "Correlation Analysis Results",
    status = "info",
    solidHeader = TRUE,
    width = NULL,
    collapsible = TRUE,
    collapsed = FALSE,
    fluidRow(
      column(6, 
             plotOutput("gtex_top20_lollipop_plot")),  # 相关分析图
      column(6, 
             DT::dataTableOutput("gtex_correlation_table"))  # 相关分析表
    ),
    fluidRow(
      column(6, 
             numericInput("gtex_pdf_height_lollipop", "PDF Height", value = 8, min = 1),
             numericInput("gtex_pdf_width_lollipop", "PDF Width", value = 10, min = 1),
             downloadButton("download_gtex_top20_lollipop", "Download Top 20 Lollipop Plot as PDF")
      ),
      column(6, 
             downloadButton("download_gtex_correlation_table", "Download Correlation Table as CSV"))
    )
  ),
  
  
  # 第四行：富集分析结果
  bs4Card(
    title = "Enrichment Analysis Results",
    status = "warning",
    solidHeader = TRUE,
    width = NULL,
    collapsible = TRUE,
    collapsed = FALSE,
    fluidRow(
      column(12,
             tabsetPanel(
               id = "gtex_enrichment_tabs",
               tabPanel("BP Output", 
                        plotOutput("gtex_BP_output"),  
                        numericInput("gtex_pdf_height_BP", "PDF Height", value = 6, min = 1),
                        numericInput("gtex_pdf_width_BP", "PDF Width", value = 8, min = 1),
                        downloadButton("gtex_download_BP", "Download BP as PDF")
               ),
               tabPanel("KEGG Legacy Output", 
                        plotOutput("gtex_KEGG_LEGACY_output"),  
                        numericInput("gtex_pdf_height_KEGG", "PDF Height", value = 6, min = 1),
                        numericInput("gtex_pdf_width_KEGG", "PDF Width", value = 8, min = 1),
                        downloadButton("gtex_download_KEGG", "Download KEGG as PDF")
               ),
               tabPanel("Reactome Output", 
                        plotOutput("gtex_REACTOME_output"),  
                        numericInput("gtex_pdf_height_REACTOME", "PDF Height", value = 6, min = 1),
                        numericInput("gtex_pdf_width_REACTOME", "PDF Width", value = 8, min = 1),
                        downloadButton("gtex_download_REACTOME", "Download Reactome as PDF")
               )
             )
      )
    )
  )
)





# atlas_tab --------------
atlas_tab <- bs4TabItem(
  tabName = "atlas",
  
  # 第一行：输入组件
  bs4Card(
    title = "Atlas Inputs",
    status = "primary",
    solidHeader = TRUE,
    width = 12,            # 占满页面宽度
    collapsible = TRUE,     # 允许收缩
    collapsed = FALSE,      # 默认展开
    fluidRow(
      column(3,             
             textInput("gene_input", "Enter Gene Name:", placeholder = "e.g., NFE2L1")
      ),
      column(3,
             selectInput("v1_input", "Select V1", 
                         choices = c('Annotation', 'Cluster'),
                         selected = 'Annotation')
      ),
      column(3,
             selectInput("v2_input", "Select V2", 
                         choices = c('All', 'Tech', 'Data_name'),
                         selected = 'All')
      )
    ),
    fluidRow(
      column(12,            
             actionButton("run_altas_analysis", "Run Analysis", class = "btn-primary")  # 运行分析按钮
      )
    )
  ),
  
  # 第二行：图片输出 (轮播图)
  bs4Card(
    title = "Atlas Plots",
    status = "info",
    solidHeader = TRUE,
    width = NULL,
    collapsible = TRUE,  
    collapsed = FALSE,    
    fluidRow(
      column(12,
             tabsetPanel(
               id = "atlas_plot_tabs",
               tabPanel("UMAP Annotation", 
                        plotOutput("atlas_UMAP_annotation_output"),  
                        numericInput("pdf_height_annotation", "PDF Height", value = 6, min = 1),
                        numericInput("pdf_width_annotation", "PDF Width", value = 6, min = 1),
                        downloadButton("download_annotation_plot", "Download UMAP as PDF")  
               ),
               tabPanel("UMAP Feature", 
                        plotOutput("atlas_UMAP_feature_output"),  
                        numericInput("pdf_height_feature", "PDF Height", value = 6, min = 1),
                        numericInput("pdf_width_feature", "PDF Width", value = 6, min = 1),
                        downloadButton("download_feature_plot", "Download Feature Plot as PDF")  
               ),
               tabPanel("CellStat Plot", 
                        plotOutput("atlas_CellStat_output"),  
                        numericInput("pdf_height_CellStat", "PDF Height", value = 6, min = 1),
                        numericInput("pdf_width_CellStat", "PDF Width", value = 6, min = 1),
                        downloadButton("download_CellStat_plot", "Download CellStat as PDF")  
               ),
               tabPanel("Feature Stat Plot", 
                        plotOutput("atlas_FeatureStat_output"),  # 新增 FeatureStat 图输出
                        numericInput("pdf_height_FeatureStat", "PDF Height", value = 6, min = 1),
                        numericInput("pdf_width_FeatureStat", "PDF Width", value = 6, min = 1),
                        downloadButton("download_FeatureStat_plot", "Download FeatureStat as PDF")
               )
             )
      )
    )
  ),
  
  # 第三行：Data_name 输出
  bs4Card(
    title = "Single Data Outputs",
    status = "warning",
    solidHeader = TRUE,
    width = NULL,
    collapsible = TRUE,  
    collapsed = FALSE,   
    fluidRow(
      column(12,
             fluidRow(
               column(6,
                      actionButton("prev_data_name", "Previous Data")  # 上一个 Data_name
               ),
               column(6,
                      actionButton("next_data_name", "Next Data")  # 下一个 Data_name
               )
             ),
             tabsetPanel(
               id = "data_name_plot_tabs",
               tabPanel("Single Data Feature Plot", 
                        plotOutput("data_name_feature_plot_output"),  
                        numericInput("pdf_height_data_name", "PDF Height", value = 6, min = 1),
                        numericInput("pdf_width_data_name", "PDF Width", value = 6, min = 1),
                        downloadButton("download_data_name_plot", "Download Data Plot as PDF")  
               )
             )
             
      )
    )
  )
)











