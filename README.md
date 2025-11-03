# Complex heatmap 
#Complexheatmap
setwd("C://Users/samsung/Desktop")
data1 <- read.table("Drug_rawdata.txt", header = TRUE, sep = "\t", row.names = 1)

sample_info <- data1[,15]
drug_data <- data1[,-15]


group_info <- data.frame(
  Sample = rownames(data1),
  Group = coldata$group
)

drug_order_df <- data.frame(Gene = rownames(data1), Group = sample_info)
drug_order_df <- drug_order_df %>% arrange(sample_info)
sorted_drug_data <- drug_data[drug_order_df$Gene, ]

sample_order <- order(sample_info) 
sorted_drug_data <- drug_data[, gene_info]
sorted_sample_info <- sample_info[sample_order] 

drug_order_df <- data.frame(Gene = rownames(sorted_drug_data), Group = sorted_sample_info)
drug_order_df <- drug_order_df %>% arrange(Group)
sorted_drug_data <- sorted_drug_data[drug_order_df$Gene, ]  


group_info <- data.frame(
  Sample = rownames(annotation_col), 
  Group = annotation_col$Group)
  
group_colors <- c("MYC" = "#F9766E", "PRC" = "#4DD3D6")
ha <- HeatmapAnnotation(Group = group_info$Group, col = list(Group = group_colors))
drug_category_map <- data.frame(
  Drug = drug$Drug 
  Category = drug$ProcessCategory)

heatmap_drugs <- data.frame(Drug =  colnames(sorted_drug_data))  
drug_category_df <- merge(heatmap_drugs, drug_category_map, by = "Drug", all.x = TRUE)
drug_category <- drug_category_df$Category



heatmap_drugs <- data.frame(Drug = colnames(sorted_drug_data))  # Heatmap의 약물 리스트
drug_category_df <- merge(heatmap_drugs, drug_category_map, by = "Drug", all.x = TRUE)
drug_category <- drug_category_df$Category
category_ha <- rowAnnotation(Category = drug_category, col = list(Category = category_colors))
sorted_drug_data <-as.matrix(sorted_drug_data)

Heatmap(t(sorted_drug_data), 
        name = "AUC", 
        top_annotation = ha, 
        right_annotation = category_ha,  
        col = colorRamp2(c(min(sorted_drug_data), mean(sorted_drug_data), max(sorted_drug_data)), 
                         c("blue", "white", "red")), 
        cluster_columns = FALSE,  
        cluster_rows = T,  
        row_km = 2,
        column_km = 2,
        show_row_names = TRUE, 
        show_column_names = TRUE, 
        column_names_rot = 90,  
        heatmap_legend_param = list(title = "AUC"))


