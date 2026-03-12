source('./tabs.R')

ui <- dashboardPage(
  # basic-------------------
  controlbar = NULL,
  footer = NULL,
  title = "IMSTAS",
  fullscreen = TRUE,
  help = TRUE,
  dark = FALSE,
  scrollToTop = TRUE,
  
  # header--------------
  header = dashboardHeader(
    title = dashboardBrand(
      title = "IMSTAS-V1",
      color = "success",
      image = "cmu_phs.png",
      href = "https://www.cmu.edu.cn/ggwsxy/"
    ),
    skin = "light",
    status = "success",
    fixed = TRUE,
    sidebarIcon = shiny::icon("bars"),
    controlbarIcon = shiny::icon("table-cells")
  ),
  
  # sidebar----------------------
  sidebar = dashboardSidebar(
    skin = "light",
    status = "success",
    elevation = 3,
    sidebarMenu(
      id = "sidebarMenuID",
      # Introduction 菜单
      menuItem("Introduction", tabName = "introduction", icon = shiny::icon("info-circle")),
      
      menuItem("RNA-seq", tabName = "rna_seq", icon = shiny::icon("dna"), 
               menuSubItem("GTEx-SM", tabName = "gtex_sm"),
               menuSubItem("DEA", tabName = "single_class"),
               menuSubItem("ICN", tabName = "all_datasets")),
      
      menuItem("Microarray", tabName = "microarray", icon = shiny::icon("microchip"), 
               menuSubItem("DEA", tabName = "micro_single_class"),
               menuSubItem("ICN", tabName = "micro_all_datasets")),
      
     
      # Sc/sn RNA-seq Atlas 菜单
      menuItem("Sc/sn RNA-seq Atlas", tabName = "atlas", icon = shiny::icon("th"))
    )
  ),
  
  
  # body----------------------
  body = dashboardBody(
    bs4TabItems(
      introduction_tab,
      
      single_class_tab,
      all_datasets_tab,
      
      micro_single_class_tab,     # 新增Microarray Single Class标签
      micro_all_datasets_tab,     # 新增Microarray All Datasets标签
      
      gtex_sm_tab,                # 新增GTEx-SM标签
      atlas_tab                   # 新增Sc/sn RNA-seq Atlas标签
    )
  )
  
  
  
  
)
