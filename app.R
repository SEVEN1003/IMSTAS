# 加载用户界面和服务器逻辑
source("global.R")
source("ui.R")
source("server.R")
# # 加载自定义绘图函数
# source("R/plot_functions.R")
# # 加载自定义的分析函数
# source("R/analysis_functions.R")


# 运行Shiny应用
shinyApp(ui = ui, server = server)



















































































