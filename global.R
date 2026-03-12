# 加载必需的库
library(shiny)
library(bs4Dash)
library(DT)
library(plotly)
library(colourpicker)
library(shinyWidgets)  
library(tidyverse)
library(readxl)
library(clinfun)
library(SCP) # plot sc and set width height
options(shiny.error = function() {
  showModal(modalDialog(
    title = "Unexpected Error",
    "An unexpected error occurred. Please try again.",
    easyClose = TRUE
  ))
})

source("R/plot_functions.R")
source('R/analysis_functions.R')

# 加载数据
load("data/All_bulk.rdata")  # 加载 All_bulk 数据
load("data/sample_info_all_final.rdata")  # 加载样本信息数据


# array 
load('data/All_array.rdata')
load('data/array_sample_info_all.rdata')

# altas

muscle_core_altas = readRDS('data/muscle_core_altas.rds')

# all genes
load('./data/all_genes.rdata')

# gtex

load('./data/gtex_muscle_clean.rdata')
gtex_gene = colnames(gtex_muscle_tpm)
# 使用read.table函数读取数据
hsa_tfs <- read.table('./data/hs_hgnc_tfs.txt', sep = "\t", stringsAsFactors = FALSE)

# load deg summary
load('./data/array_deg_summary.rdata')
load('./data/bulk_deg_summary.rdata')
load('./data/gtex_trend_summary.rdata')

