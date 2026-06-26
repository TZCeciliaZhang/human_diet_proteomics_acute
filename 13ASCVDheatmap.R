#OXIDATION#################

dirname <- "20240729ASCVDheatmap"
dir.create(dirname)


library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)


dep <- read.csv("diet_dep.csv")


calculate_diff_proteins <- function(data, diet_label) {
  sel <- data[data$diet == diet_label, ]
  diff_proteins <- sel %>%
    rename_with(~ paste0("log2_", .), starts_with("fd")) %>%
    rowwise() %>%
    mutate(
      fd_values = list(exp(c_across(starts_with("log2_fd")))),
      fd_geo_mean = exp(mean(log(unlist(fd_values)), na.rm = TRUE))
    ) %>%
    ungroup() %>%
    mutate(
      logFC = log2(ifelse(fd_geo_mean > 0, fd_geo_mean, 1e-1)),
      pvalue = min(c_across(starts_with("pval")), na.rm = TRUE)
    ) %>%
    select(SYMBOL, logFC, pvalue)
  return(diff_proteins)
}


diff_proteins_High_glucose <- calculate_diff_proteins(dep, "A")
diff_proteins_Ketogenic <- calculate_diff_proteins(dep, "B")
diff_proteins_HFHC <- calculate_diff_proteins(dep, "C")
diff_proteins_LFLC <- calculate_diff_proteins(dep, "E")


target_proteins <- c("SOD1", "TXN", "IDH1", "CAT", "PRDX1")

diff_list <- list(
  High_glucose = diff_proteins_High_glucose,
  Ketogenic = diff_proteins_Ketogenic,
  HFHC = diff_proteins_HFHC,
  LFLC = diff_proteins_LFLC
)


logFC_data <- data.frame(SYMBOL = target_proteins)
pvalue_data <- data.frame(SYMBOL = target_proteins)

for(diet in names(diff_list)) {
  diff_data <- diff_list[[diet]] %>%
    filter(SYMBOL %in% target_proteins)
  
  logFC_data <- logFC_data %>%
    left_join(diff_data %>% select(SYMBOL, logFC) %>% rename(!!diet := logFC), by = "SYMBOL")
  
  pvalue_data <- pvalue_data %>%
    left_join(diff_data %>% select(SYMBOL, pvalue) %>% rename(!!diet := pvalue), by = "SYMBOL")
}


logFC_data[is.na(logFC_data)] <- 0
pvalue_data[is.na(pvalue_data)] <- 1


rownames(logFC_data) <- logFC_data$SYMBOL
logFC_data <- logFC_data[,-1]

rownames(pvalue_data) <- pvalue_data$SYMBOL
pvalue_data <- pvalue_data[,-1]


logFC_data <- logFC_data[, c("Ketogenic", "LFLC", "HFHC", "High_glucose")]
pvalue_data <- pvalue_data[, c("Ketogenic", "LFLC", "HFHC", "High_glucose")]


pvalue_col_fun <- function(pvalues) {
  ifelse(pvalues > 0.05, "white", "black")
}

ht <- Heatmap(
  as.matrix(logFC_data), 
  name = "logFC", 
  row_names_side = "left",  
  cluster_rows = TRUE, 
  cluster_columns = FALSE,   
  show_row_dend = FALSE,  
  show_column_dend = FALSE, 
  border = TRUE,
  rect_gp = gpar(col = "black", lwd = 1), 
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (pvalue_data[i, j] <= 0.05) {
      grid.text(
        sprintf("%.2f", logFC_data[i, j]), 
        x, 
        y, 
        gp = gpar(col = pvalue_col_fun(pvalue_data[i, j]), fontsize = 10)
      )
    }
  },
  column_names_rot = 30  
)

# 保存热图到文件
tiff(paste0(dirname,"/oxidation-heatmap.tiff"), res = 300, width = 1000, height = 800)  
draw(ht, heatmap_legend_side = "right")      
dev.off()                                  




dirname <- "20240729ASCVDheatmap"
dir.create(dirname)


library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)


dep <- read.csv("diet_dep.csv")


calculate_diff_proteins <- function(data, diet_label) {
  sel <- data[data$diet == diet_label, ]
  diff_proteins <- sel %>%
    rename_with(~ paste0("log2_", .), starts_with("fd")) %>%
    rowwise() %>%
    mutate(
      fd_values = list(exp(c_across(starts_with("log2_fd")))),
      fd_geo_mean = exp(mean(log(unlist(fd_values)), na.rm = TRUE))
    ) %>%
    ungroup() %>%
    mutate(
      logFC = log2(ifelse(fd_geo_mean > 0, fd_geo_mean, 1e-1)),
      pvalue = min(c_across(starts_with("pval")), na.rm = TRUE)
    ) %>%
    select(SYMBOL, logFC, pvalue)
  return(diff_proteins)
}


diff_proteins_High_glucose <- calculate_diff_proteins(dep, "A")
diff_proteins_Ketogenic <- calculate_diff_proteins(dep, "B")
diff_proteins_HFHC <- calculate_diff_proteins(dep, "C")
diff_proteins_LFLC <- calculate_diff_proteins(dep, "E")


target_proteins <- c("CA1","CA2","FTL","FTH1","APOA5")

diff_list <- list(
  High_glucose = diff_proteins_High_glucose,
  Ketogenic = diff_proteins_Ketogenic,
  HFHC = diff_proteins_HFHC,
  LFLC = diff_proteins_LFLC
)


logFC_data <- data.frame(SYMBOL = target_proteins)
pvalue_data <- data.frame(SYMBOL = target_proteins)

for(diet in names(diff_list)) {
  diff_data <- diff_list[[diet]] %>%
    filter(SYMBOL %in% target_proteins)
  
  logFC_data <- logFC_data %>%
    left_join(diff_data %>% select(SYMBOL, logFC) %>% rename(!!diet := logFC), by = "SYMBOL")
  
  pvalue_data <- pvalue_data %>%
    left_join(diff_data %>% select(SYMBOL, pvalue) %>% rename(!!diet := pvalue), by = "SYMBOL")
}


logFC_data[is.na(logFC_data)] <- 0
pvalue_data[is.na(pvalue_data)] <- 1


rownames(logFC_data) <- logFC_data$SYMBOL
logFC_data <- logFC_data[,-1]

rownames(pvalue_data) <- pvalue_data$SYMBOL
pvalue_data <- pvalue_data[,-1]


logFC_data <- logFC_data[, c("Ketogenic", "LFLC", "HFHC", "High_glucose")]
pvalue_data <- pvalue_data[, c("Ketogenic", "LFLC", "HFHC", "High_glucose")]


pvalue_col_fun <- function(pvalues) {
  ifelse(pvalues > 0.05, "white", "black")
}

ht <- Heatmap(
  as.matrix(logFC_data), 
  name = "logFC", 
  row_names_side = "right", 
  cluster_rows = FALSE, 
  cluster_columns = FALSE,  
  show_row_dend = FALSE,   
  show_column_dend = FALSE, 
  border = TRUE,
  rect_gp = gpar(col = "black", lwd = 1),  
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")), 
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (pvalue_data[i, j] <= 0.05) {
      grid.text(
        sprintf("%.2f", logFC_data[i, j]), 
        x, 
        y, 
        gp = gpar(col = pvalue_col_fun(pvalue_data[i, j]), fontsize = 10)
      )
    }
  },
  column_names_rot = 30  
)


tiff(paste0(dirname,"/lipo-heatmap.tiff"), res = 300, width = 1000, height = 800)  
draw(ht, heatmap_legend_side = "left")     
dev.off()                                   


dirname <- "20240729ASCVDheatmap"
dir.create(dirname)


library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)


dep <- read.csv("diet_dep.csv")


calculate_diff_proteins <- function(data, diet_label) {
  sel <- data[data$diet == diet_label, ]
  diff_proteins <- sel %>%
    rename_with(~ paste0("log2_", .), starts_with("fd")) %>%
    rowwise() %>%
    mutate(
      fd_values = list(exp(c_across(starts_with("log2_fd")))),
      fd_geo_mean = exp(mean(log(unlist(fd_values)), na.rm = TRUE))
    ) %>%
    ungroup() %>%
    mutate(
      logFC = log2(ifelse(fd_geo_mean > 0, fd_geo_mean, 1e-1)),
      pvalue = min(c_across(starts_with("pval")), na.rm = TRUE)
    ) %>%
    select(SYMBOL, logFC, pvalue)
  return(diff_proteins)
}


diff_proteins_High_glucose <- calculate_diff_proteins(dep, "A")
diff_proteins_Ketogenic <- calculate_diff_proteins(dep, "B")
diff_proteins_HFHC <- calculate_diff_proteins(dep, "C")
diff_proteins_LFLC <- calculate_diff_proteins(dep, "E")


target_proteins <- c("PPBP","SRGN","THBS1","PF4V1")

diff_list <- list(
  High_glucose = diff_proteins_High_glucose,
  Ketogenic = diff_proteins_Ketogenic,
  HFHC = diff_proteins_HFHC,
  LFLC = diff_proteins_LFLC
)


logFC_data <- data.frame(SYMBOL = target_proteins)
pvalue_data <- data.frame(SYMBOL = target_proteins)

for(diet in names(diff_list)) {
  diff_data <- diff_list[[diet]] %>%
    filter(SYMBOL %in% target_proteins)
  
  logFC_data <- logFC_data %>%
    left_join(diff_data %>% select(SYMBOL, logFC) %>% rename(!!diet := logFC), by = "SYMBOL")
  
  pvalue_data <- pvalue_data %>%
    left_join(diff_data %>% select(SYMBOL, pvalue) %>% rename(!!diet := pvalue), by = "SYMBOL")
}


logFC_data[is.na(logFC_data)] <- 0
pvalue_data[is.na(pvalue_data)] <- 1


rownames(logFC_data) <- logFC_data$SYMBOL
logFC_data <- logFC_data[,-1]

rownames(pvalue_data) <- pvalue_data$SYMBOL
pvalue_data <- pvalue_data[,-1]


logFC_data <- logFC_data[, c("Ketogenic", "LFLC", "HFHC", "High_glucose")]
pvalue_data <- pvalue_data[, c("Ketogenic", "LFLC", "HFHC", "High_glucose")]


pvalue_col_fun <- function(pvalues) {
  ifelse(pvalues > 0.05, "white", "black")
}

ht <- Heatmap(
  as.matrix(logFC_data), 
  name = "logFC", 
  row_names_side = "right", 
  cluster_rows = FALSE, 
  cluster_columns = FALSE,  
  show_row_dend = FALSE,   
  show_column_dend = FALSE,
  border = TRUE,
  rect_gp = gpar(col = "black", lwd = 1),  
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (pvalue_data[i, j] <= 0.05) {
      grid.text(
        sprintf("%.2f", logFC_data[i, j]), 
        x, 
        y, 
        gp = gpar(col = pvalue_col_fun(pvalue_data[i, j]), fontsize = 10)
      )
    }
  },
  column_names_rot = 30  
)


tiff(paste0(dirname,"/platelet-heatmap.tiff"), res = 300, width = 1000, height = 800) 
draw(ht, heatmap_legend_side = "left")     
dev.off()                                


dirname <- "20240729ASCVDheatmap"
dir.create(dirname)


library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)


dep <- read.csv("diet_dep.csv")


calculate_diff_proteins <- function(data, diet_label) {
  sel <- data[data$diet == diet_label, ]
  diff_proteins <- sel %>%
    rename_with(~ paste0("log2_", .), starts_with("fd")) %>%
    rowwise() %>%
    mutate(
      fd_values = list(exp(c_across(starts_with("log2_fd")))),
      fd_geo_mean = exp(mean(log(unlist(fd_values)), na.rm = TRUE))
    ) %>%
    ungroup() %>%
    mutate(
      logFC = log2(ifelse(fd_geo_mean > 0, fd_geo_mean, 1e-1)),
      pvalue = min(c_across(starts_with("pval")), na.rm = TRUE)
    ) %>%
    select(SYMBOL, logFC, pvalue)
  return(diff_proteins)
}


diff_proteins_High_glucose <- calculate_diff_proteins(dep, "A")
diff_proteins_Ketogenic <- calculate_diff_proteins(dep, "B")
diff_proteins_HFHC <- calculate_diff_proteins(dep, "C")
diff_proteins_LFLC <- calculate_diff_proteins(dep, "E")


target_proteins <- c("CRP","ORM1","MIF","S100A4","S100A7")

diff_list <- list(
  High_glucose = diff_proteins_High_glucose,
  Ketogenic = diff_proteins_Ketogenic,
  HFHC = diff_proteins_HFHC,
  LFLC = diff_proteins_LFLC
)


logFC_data <- data.frame(SYMBOL = target_proteins)
pvalue_data <- data.frame(SYMBOL = target_proteins)

for(diet in names(diff_list)) {
  diff_data <- diff_list[[diet]] %>%
    filter(SYMBOL %in% target_proteins)
  
  logFC_data <- logFC_data %>%
    left_join(diff_data %>% select(SYMBOL, logFC) %>% rename(!!diet := logFC), by = "SYMBOL")
  
  pvalue_data <- pvalue_data %>%
    left_join(diff_data %>% select(SYMBOL, pvalue) %>% rename(!!diet := pvalue), by = "SYMBOL")
}


logFC_data[is.na(logFC_data)] <- 0
pvalue_data[is.na(pvalue_data)] <- 1


rownames(logFC_data) <- logFC_data$SYMBOL
logFC_data <- logFC_data[,-1]

rownames(pvalue_data) <- pvalue_data$SYMBOL
pvalue_data <- pvalue_data[,-1]


logFC_data <- logFC_data[, c("Ketogenic", "LFLC", "HFHC", "High_glucose")]
pvalue_data <- pvalue_data[, c("Ketogenic", "LFLC", "HFHC", "High_glucose")]


pvalue_col_fun <- function(pvalues) {
  ifelse(pvalues > 0.05, "white", "black")
}

ht <- Heatmap(
  as.matrix(logFC_data), 
  name = "logFC", 
  row_names_side = "right", 
  cluster_rows = FALSE, 
  cluster_columns = FALSE,  
  show_row_dend = FALSE,  
  show_column_dend = FALSE, 
  border = TRUE,
  rect_gp = gpar(col = "black", lwd = 1), 
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")), 
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (pvalue_data[i, j] <= 0.05) {
      grid.text(
        sprintf("%.2f", logFC_data[i, j]), 
        x, 
        y, 
        gp = gpar(col = pvalue_col_fun(pvalue_data[i, j]), fontsize = 10)
      )
    }
  },
  column_names_rot = 30  
)


tiff(paste0(dirname,"/inflam-heatmap.tiff"), res = 300, width = 1000, height = 800) 
draw(ht, heatmap_legend_side = "left")    
dev.off()                                 


dirname <- "20240729ASCVDheatmap"
dir.create(dirname)


library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)


dep <- read.csv("diet_dep.csv")


calculate_diff_proteins <- function(data, diet_label) {
  sel <- data[data$diet == diet_label, ]
  diff_proteins <- sel %>%
    rename_with(~ paste0("log2_", .), starts_with("fd")) %>%
    rowwise() %>%
    mutate(
      fd_values = list(exp(c_across(starts_with("log2_fd")))),
      fd_geo_mean = exp(mean(log(unlist(fd_values)), na.rm = TRUE))
    ) %>%
    ungroup() %>%
    mutate(
      logFC = log2(ifelse(fd_geo_mean > 0, fd_geo_mean, 1e-1)),
      pvalue = min(c_across(starts_with("pval")), na.rm = TRUE)
    ) %>%
    select(SYMBOL, logFC, pvalue)
  return(diff_proteins)
}


diff_proteins_High_glucose <- calculate_diff_proteins(dep, "A")
diff_proteins_Ketogenic <- calculate_diff_proteins(dep, "B")
diff_proteins_HFHC <- calculate_diff_proteins(dep, "C")
diff_proteins_LFLC <- calculate_diff_proteins(dep, "E")


target_proteins <- c("LMNA","ARPC2","FLNA","RAP1B")

diff_list <- list(
  High_glucose = diff_proteins_High_glucose,
  Ketogenic = diff_proteins_Ketogenic,
  HFHC = diff_proteins_HFHC,
  LFLC = diff_proteins_LFLC
)


logFC_data <- data.frame(SYMBOL = target_proteins)
pvalue_data <- data.frame(SYMBOL = target_proteins)

for(diet in names(diff_list)) {
  diff_data <- diff_list[[diet]] %>%
    filter(SYMBOL %in% target_proteins)
  
  logFC_data <- logFC_data %>%
    left_join(diff_data %>% select(SYMBOL, logFC) %>% rename(!!diet := logFC), by = "SYMBOL")
  
  pvalue_data <- pvalue_data %>%
    left_join(diff_data %>% select(SYMBOL, pvalue) %>% rename(!!diet := pvalue), by = "SYMBOL")
}


logFC_data[is.na(logFC_data)] <- 0
pvalue_data[is.na(pvalue_data)] <- 1


rownames(logFC_data) <- logFC_data$SYMBOL
logFC_data <- logFC_data[,-1]

rownames(pvalue_data) <- pvalue_data$SYMBOL
pvalue_data <- pvalue_data[,-1]


logFC_data <- logFC_data[, c("Ketogenic", "LFLC", "HFHC", "High_glucose")]
pvalue_data <- pvalue_data[, c("Ketogenic", "LFLC", "HFHC", "High_glucose")]


pvalue_col_fun <- function(pvalues) {
  ifelse(pvalues > 0.05, "white", "black")
}

ht <- Heatmap(
  as.matrix(logFC_data), 
  name = "logFC", 
  row_names_side = "left",  
  cluster_rows = FALSE, 
  cluster_columns = FALSE,  
  show_row_dend = FALSE,  
  show_column_dend = FALSE, 
  border = TRUE,
  rect_gp = gpar(col = "black", lwd = 1),  
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")), 
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (pvalue_data[i, j] <= 0.05) {
      grid.text(
        sprintf("%.2f", logFC_data[i, j]), 
        x, 
        y, 
        gp = gpar(col = pvalue_col_fun(pvalue_data[i, j]), fontsize = 10)
      )
    }
  },
  column_names_rot = 30 
)


tiff(paste0(dirname,"/cell-heatmap.tiff"), res = 300, width = 1000, height = 800) 
draw(ht, heatmap_legend_side = "right")     
dev.off()                                 

