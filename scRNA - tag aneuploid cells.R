#this script will tag cells as aneuploid or not based on the average expression of chromosomes that are expected to be aneuploid

#inputs
seuratObj <- seuratObj

#parameters 
subsample = FALSE #If true, will subsample 10% of the cells to perform the analysis. Usefull to speed up the script.
aneu_chromo <- list("B12" = c("01", "05", "10", "13", "20", "23"),
                    "E4" = c("13", "20", "23", "26"),
                    "G3" = c("20", "23", "26", "29"),
                    "G6" = c("01", "09", "13", "16", "20", "22", "23", "24", "26", "35"))

#aneu_chromo <- c("02")
#libraries
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(ggridges)
library(Seurat)
library(viridis)

#custom_functions
scale_vec <- function(x){
  x <- (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
  return(x)
}


#code
rm(genes_meta)
rm(expr_table)
rm(sum_table)
genes_meta <- seuratObj@assays$RNA@meta.features #get the metadata for the genes
expr_table <- as.matrix(seuratObj@assays$RNA@data) #get the expression matrix
if(subsample){
  cells <- c(1:ncol(expr_table))
  cells <- sample(cells, size = round(length(cells)*0.1), replace = FALSE)
  expr_table <- expr_table[,cells]
  rm(cells)
}

#add chromosome information in case it is missing
if(!"chromo" %in% colnames(seuratObj@assays$RNA@meta.features)){
  seuratObj@assays$RNA@meta.features <- seuratObj@assays$RNA@meta.features %>%
    mutate(chromo = substr(id, 7, 8))
}

#set 0s to NAs
#expr_table[expr_table == 0] <- NA

#convert to dataframe and add a column with the genes names
expr_table <- expr_table %>%
  as.data.frame()%>%
  rownames_to_column("gene")

#check in the metadata to which chromosome each gene belongs
expr_table$chromo <- genes_meta$chromo[match(expr_table$gene, genes_meta$name)]

#calculating average expression of each chromosome####
sum_table <- expr_table  %>%
  #convert to tidy dataframe
  pivot_longer(cols = c(2:(ncol(expr_table)-1)), names_to = "cell", values_to = "expression") %>%
  #remove kDNA
  filter(chromo != 37) %>%
  #for each cell and each chromosome
  group_by(cell, chromo)%>%
  #calculate the average expression of each chromosome in each cell
  summarise(mean_expr = mean(expression, na.rm = T)) %>%
  #summarise(mean_expr = sum(expression, na.rm = T)/n())%>%
  ungroup()

#ploting heatmap with calculated average expression of each chromosome####
# to_plot <- sum_table %>%
#   pivot_wider(names_from = cell, values_from = mean_expr)%>%
#   as.data.frame()
# 
# rownames(to_plot) <- to_plot$chromo
# to_plot$chromo <- NULL
# to_plot <- as.matrix(to_plot)
# 
# pheatmap(to_plot,
#          show_colnames = F, 
#          scale = "column",
#          clustering_method = "ward.D2",
#          cutree_rows = 2)

#Calculatiing for each cell the ration between the average expression of aneuploid chromosomes vs normal chromosomes####
sum_table$sample <- seuratObj@meta.data$sample[match(sum_table$cell, rownames(seuratObj@meta.data))]

cell_aneu_ratio <- sum_table %>%
  group_by(sample)%>%
  mutate(aneu_chromo = chromo %in% unlist(aneu_chromo[which(names(aneu_chromo) == sample)]))%>%
  mutate(chromo_group = ifelse(aneu_chromo, "aneuploid", "normal"))%>%
  group_by(cell, chromo_group)%>%
  summarise(mean_expr = mean(mean_expr))%>%
  pivot_wider(names_from = chromo_group, values_from = mean_expr)%>%
  mutate(aneu_ratio = aneuploid/normal)%>%
  ungroup()

#assigning the aneu_ratio to the cells metadata in the seurat seuratObj
seuratObj@meta.data$aneu_ratio <- NA
seuratObj@meta.data$aneu_ratio <- cell_aneu_ratio$aneu_ratio[match(cell_aneu_ratio$cell, rownames(seuratObj@meta.data))]


#plotting the aneu_ratio distribution
seuratObj@meta.data %>%
  ggplot(aes(x = aneu_ratio, y = sample, fill = sample))+
  geom_density_ridges(scale = 3, alpha= 0.5)+
  guides(fill = "none")+
  labs(x = "Aneuploid/disomic ratio", y = "Sample")

col_vec <- c("#f7f7f7", "#fcffa6", "yellow", "orange", "red", "dark red", "#730000", "#330000")
col_vec <- rev(viridis::magma(8))
col_vec <- c("#f7f7f7",  "#f7fabb", "#e6eb94", "#638a85", "#485463", "#2a3440", "#182029")
#plotting the aneu_ratio in a UMAP
patchwork::wrap_plots(FeaturePlot(seuratObj, 
            feature = "aneu_ratio", 
            reduction = "umap", 
            split.by = "sample", 
            keep.scale = NULL,
            combine = FALSE))&
  theme_minimal()&
  theme(panel.grid = element_blank())&
  scale_color_gradientn(colors = col_vec, limits = c(0.5, 2.1))
  #scale_color_gradientn(colors = rev(viridis::magma(8)))


#Now it will use kmeans to determine if cells are aneuploid or normal
  #NTS: for now it assumes that both normal and aneuploid cells co-exist in the data. Fix it later.
clusters <- data.frame(value = seuratObj@meta.data$aneu_ratio)
clusters$group <- kmeans(clusters$value, centers = 2)$cluster
group_name <- unique(clusters$group)
group_name[which(group_name == max(group_name))] <- "aneuploid"
group_name[which(group_name == min(group_name))] <- "normal"
clusters$group_name <- group_name[clusters$group]
clusters %>%
  ggplot(aes(x = value, fill = group_name))+
  geom_density(alpha = 0.5)

#adding the group name 
seuratObj@meta.data$aneu_group <- clusters$group_name
#adding the group name as an identity class
seuratObj <- SetIdent(seuratObj, value = "aneu_group")
#plotting the result
DimPlot(seuratObj, reduction = "umap", group.by = "aneu_group", split.by = "sample")&
  scale_color_manual(values = c("#226185", "#e6d99a"))

#finding marker genes between aneuploid and diploid cells
marker_genes <- seuratObj %>%
  SetIdent(value = "aneu_group")%>%
  FindAllMarkers(slot = "scale.data", only.pos = FALSE)
genes_to_plot <- marker_genes %>%
  group_by(cluster) %>%
  arrange(p_val_adj)%>%
  mutate(position = 1:n()) %>%
  filter(position <= 20)
DoHeatmap(seuratObj, features = genes_to_plot$gene, group.by = "aneu_group", slot = "scale.data")

