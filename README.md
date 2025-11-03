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


# Archytype Analysis
dir.create("results/signatures", showWarnings = FALSE)

# -------------------- Packages --------------------
req_cran <- c("tidyverse","matrixStats","glmnet","ggrepel")
for(p in req_cran) if(!requireNamespace(p, quietly=TRUE)) install.packages(p)
invisible(lapply(req_cran, library, character.only=TRUE))

if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
bio_pkgs <- c("edgeR","archetypes","depmap","ExperimentHub")
for (bp in bio_pkgs) if (!requireNamespace(bp, quietly=TRUE)) BiocManager::install(bp, ask=FALSE, update=FALSE)
BiocManager::install("depmap")
# -------------------- 0) Load data --------------------
grp <- readr::read_tsv("MYC_PRC_group_info.txt", col_types = "cc") %>%
  setNames(c("sample","subtype")) %>%
  mutate(subtype = toupper(subtype))

cnt <- readr::read_tsv("MYC_PRC_rawData_20250115.txt")
if (!"gene" %in% names(cnt)) {
  nm <- names(cnt); names(cnt)[1] <- "gene"
}

# 공통 샘플 정렬
samples <- intersect(grp$sample, setdiff(names(cnt), "gene"))
stopifnot(length(samples) >= 20)
grp <- grp %>% filter(sample %in% samples) %>% slice(match(samples, sample))
expr <- cnt %>% select(all_of(c("gene", samples))) %>% distinct(gene, .keep_all = TRUE)
mat <- as.matrix(expr[,-1]); rownames(mat) <- expr$gene


# CPM + log2; 
libsize <- colSums(mat)
cpm <- t(t(mat) / libsize * 1e6)
logcpm <- log2(cpm + 1)
vars <- matrixStats::rowVars(logcpm); logcpm <- logcpm[vars > 0, , drop=FALSE]

# -------------------- 2) Conformal (CV+, 5-fold, alpha=0.10) --------------------
alpha <- 0.10
y <- ifelse(grp$subtype == "MYC", 1, 0)




X <- t(logcpm[top, , drop=FALSE])   # samples x genes
X <- scale(X)                   
X <- X[, apply(X, 2, function(z) sd(z, na.rm=TRUE) > 0), drop=FALSE]
X[!is.finite(X)] <- 0


stopifnot(nrow(X) == length(y))     # 행 = 샘플 수




K <- 5
foldid <- sample(rep_len(1:K, length.out=length(y)))

p_hat <- rep(NA_real_, length(y)); nonconf <- rep(NA_real_, length(y))
for (k in 1:K) {
  tr <- foldid != k; te <- !tr
  fit <- glmnet::cv.glmnet(X[tr,], y[tr], family="binomial", alpha=1, type.measure="class")
  pr <- as.numeric(predict(fit, X[te,], s="lambda.min", type="response"))
  p_hat[te] <- pr
  nonconf[te] <- ifelse(y[te]==1, 1 - pr, pr)
}
t_MYC <- 1 - p_hat; t_PRC <- p_hat
n <- length(y)
count_ge <- function(v, thr) sum(v >= thr)
p_MYC <- sapply(t_MYC, function(th) (count_ge(nonconf, th) + 1) / (n + 1))
p_PRC <- sapply(t_PRC, function(th) (count_ge(nonconf, th) + 1) / (n + 1))

predset <- ifelse(p_MYC > alpha & p_PRC > alpha, "{MYC,PRC}",
                  ifelse(p_MYC > alpha, "{MYC}",
                         ifelse(p_PRC > alpha, "{PRC}", "{}")))



conf_tab <- conf_tab %>%
  dplyr::mutate(prediction_set = as.character(prediction_set))


conf_counts <- conf_tab %>%
  dplyr::count(prediction_set, name = "n")

ord <- c("{MYC}","{PRC}","{MYC,PRC}","{}")
conf_counts <- conf_counts %>%
  dplyr::mutate(prediction_set = factor(prediction_set, levels = ord)) %>%
  dplyr::arrange(prediction_set)


readr::write_csv(conf_counts, "박사님/conformal_set_counts.csv")





ggplot(conf_counts, aes(prediction_set, n)) +
  geom_col() + labs(title=paste0("Conformal sets (alpha=",alpha,")"), x="", y="Count") +
  theme_classic()
ggsave("results/conformal_sets_bar.png", width=4.6, height=3.3, dpi=300)


# archetypes::stepArchetypes (K=2)
library(archetypes)
coef_mat <- coef(aa)


coef_mat <- as.matrix(coef_mat)
storage.mode(coef_mat) <- "double"           # numeric 보장


grp2 <- grp[ match(rownames(coef_mat), grp$sample), , drop=FALSE ]
stopifnot(all(rownames(coef_mat) == grp2$sample))


y_bin <- as.numeric(grp2$subtype == "MYC")

r1 <- suppressWarnings(cor(coef_mat[, 1], y_bin, use = "complete.obs"))
r2 <- suppressWarnings(cor(coef_mat[, 2], y_bin, use = "complete.obs"))

myc_side <- if (r1 >= r2) 1 else 2
alpha_MYC <- coef_mat[, myc_side]

arch_df <- tibble::tibble(
  sample   = rownames(coef_mat),
  alpha_MYC = as.numeric(alpha_MYC),
  subtype   = grp2$subtype
)
readr::write_csv(arch_df, "results/archetype_alpha_k2.csv")
coef_mat <- coef(aa)  # samples x 2, 

myc_side <- if (cor(coef_mat[,1], as.numeric(grp$subtype=="MYC")) >=
                cor(coef_mat[,2], as.numeric(grp$subtype=="MYC"))) 1 else 2
alpha_MYC <- coef_mat[, myc_side]
arch_df <- tibble(sample=rownames(coef_mat), alpha_MYC=alpha_MYC,
                  subtype=grp$subtype)
readr::write_csv(arch_df, "results/archetype_alpha_k2.csv")

ggplot(arch_df, aes(subtype, alpha_MYC, fill=subtype)) +
  geom_boxplot(outlier.shape=NA, alpha=0.6) + geom_jitter(width=0.12, alpha=0.7) +
  scale_y_continuous(limits=c(0,1)) +
  labs(title="Archetypal analysis (K=2): alpha(MYC)", x="", y="alpha (MYC archetype)") +
  theme_classic() + theme(legend.position="none")
ggsave("results/archetype_alpha_box.png", width=4.6, height=3.3, dpi=300)


# TCGA 
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", "GSVA", "GSEABase"))

BiocManager::install("sva")
BiocManager::install("escape")

library(TCGAbiolinks)
library(SummarizedExperiment)
library(GSVA)
library(GSEABase)



if(!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github('sidmall/PDSclassifier')

library(TCGAbiolinks)


options <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification"
)

options

data_path <- "C:/Users/USER/Desktop/GDCdata"


files <- list.files(data_path, pattern = "\\.tsv$", recursive = TRUE, full.names = TRUE)

data_list <- lapply(files, function(file) {
  tryCatch(
    {
      
      read.table(file, header = TRUE, sep = "\t")
    },
    error = function(e) {
      cat("Error in file:", file, "\n", e$message, "\n")
      return(NULL)
    }
  )
})


data_list <- Filter(Negate(is.null), data_list)


if (length(data_list) > 0) {
  print(head(data_list[[1]]))
} else {
  cat("No valid files loaded.\n")
}


readcount_data <- do.call(cbind, data_list)



rna_df2_filtered <- cbind(Gene = rna_df2_filtered$Gene, rna_df2_filtered)
rownames(rna_df2_filtered) <- NULL  

rna_df2_filtered$Gene = NULL

head(readcount_data)
rna_df2_filtered$Gene <- rownames(rna_df2_filtered)
rna_df2_filtered <- rna_df2_filtered[, -ncol(rna_df2_filtered)]

library(PDSclassifier)
pds_calls <- PDSpredict(rna_df2_filtered, species = 'human', threshold = 0.6)
smi_data <- calculateSMI(as.matrix(testdata[,-1]), datatype = "bulk", species = "human")
