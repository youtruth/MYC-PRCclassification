# Complex heatmap 
#Complexheatmap
setwd("C://Users/samsung/Desktop")
data1 <- read.table("Drug_rawdata.txt", header = TRUE, sep = "\t", row.names = 1)
sample_info <- data1[,15]
drug_data <- data1[,-15]
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

heatmap_drugs <- data.frame(Drug = colnames(sorted_drug_data))  
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


req_cran <- c("tidyverse","matrixStats","glmnet","ggrepel")
for(p in req_cran) if(!requireNamespace(p, quietly=TRUE)) install.packages(p)
invisible(lapply(req_cran, library, character.only=TRUE))

if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
bio_pkgs <- c("edgeR","archetypes","depmap","ExperimentHub")
for (bp in bio_pkgs) if (!requireNamespace(bp, quietly=TRUE)) BiocManager::install(bp, ask=FALSE, update=FALSE)
BiocManager::install("depmap")

grp <- readr::read_tsv("MYC_PRC_group_info.txt", col_types = "cc") %>%
  setNames(c("sample","subtype")) %>%
  mutate(subtype = toupper(subtype))
cnt <- readr::read_tsv("MYC_PRC_rawData_20250115.txt")
if (!"gene" %in% names(cnt)) {
  nm <- names(cnt); names(cnt)[1] <- "gene"
}

samples <- intersect(grp$sample, setdiff(names(cnt), "gene"))
stopifnot(length(samples) >= 20)
grp <- grp %>% filter(sample %in% samples) %>% slice(match(samples, sample))
expr <- cnt %>% select(all_of(c("gene", samples))) %>% distinct(gene, .keep_all = TRUE)
mat <- as.matrix(expr[,-1]); rownames(mat) <- expr$gene
libsize <- colSums(mat)
cpm <- t(t(mat) / libsize * 1e6)
logcpm <- log2(cpm + 1)
vars <- matrixStats::rowVars(logcpm); logcpm <- logcpm[vars > 0, , drop=FALSE]

# Conformal (CV+, 5-fold, alpha=0.10)
alpha <- 0.10
y <- ifelse(grp$subtype == "MYC", 1, 0)
X <- t(logcpm[top, , drop=FALSE])   # samples x genes
X <- scale(X)                   
X <- X[, apply(X, 2, function(z) sd(z, na.rm=TRUE) > 0), drop=FALSE]
X[!is.finite(X)] <- 0
stopifnot(nrow(X) == length(y))    
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


readr::write_csv(conf_counts, "/conformal_set_counts.csv")


ggplot(conf_counts, aes(prediction_set, n)) +
  geom_col() + labs(title=paste0("Conformal sets (alpha=",alpha,")"), x="", y="Count") +
  theme_classic()
ggsave("results/conformal_sets_bar.png", width=4.6, height=3.3, dpi=300)



# archetypes::stepArchetypes (K=2)
library(archetypes)
coef_mat <- coef(aa)


coef_mat <- as.matrix(coef_mat)
storage.mode(coef_mat) <- "double"          

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

# Multi-Omics
# mixOmics, limma, ggplot2

suppressPackageStartupMessages({
  library(mixOmics)
  library(limma)
})

# 1) Load data
# path/to files (tab-delimited; rows = features, cols = samples)
expr_path   <- "rawdata.txt"  
meta_path   <- "sample_list.txt"      
X <- read.table(expr_path, header = TRUE, sep = "\t", na.strings = "NA", fill = TRUE, row.names = 1, check.names = FALSE)
meta <- read.table(meta_path, header = TRUE, sep = "\t", na.strings = "NA", fill = TRUE, row.names = 1, check.names = FALSE)

common_samples <- intersect(colnames(X), rownames(meta))
X   <- X[, common_samples, drop = FALSE]
meta <- meta[common_samples, , drop = FALSE]

group <- factor(meta$group)

# 2) PLS-DA / sPLS-DA (mixOmics)
plsda_model <- plsda(X = t(X), Y = group, ncomp = 2)  # mixOmics expects samples in rows
plotIndiv(plsda_model,
          group = group,
          title = "PLS-DA (comp 1–2)",
          legend = TRUE,
          ellipse = TRUE, ellipse.level = 0.95)

# sPLS-DA on the same matrix
keepX <- c(30, 30)
splsda_model <- splsda(X = t(X), Y = group, ncomp = 2, keepX = keepX)
plotIndiv(splsda_model,
          group = group,
          title = "sPLS-DA (comp 1–2)",
          legend = TRUE,
          ellipse = TRUE, ellipse.level = 0.95)

design <- model.matrix(~ group)
fit <- lmFit(X, design)
fit <- eBayes(fit)
res_all <- topTable(fit, coef = 2, number = Inf)          # coef=2 corresponds to group factor
res_sig <- topTable(fit, coef = 2, p.value = 0.005, number = Inf)

write.table(res_all, file = "limma_results_all.txt", sep = "\t", quote = FALSE)
write.table(res_sig, file = "limma_results_p0.005.txt", sep = "\t", quote = FALSE)

omixs_raw <- read.table("DIABLO_raw.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

if ("X5FU" %in% colnames(omixs_raw)) {
  colnames(omixs_raw)[colnames(omixs_raw) == "X5FU"] <- "5FU"
}


drug_block    <- omixs_raw[, 23:,   drop = FALSE]
rna_block     <- omixs_raw[, 54:, drop = FALSE]
pathway_block <- omixs_raw[, 70:, drop = FALSE]  # verify this range is what you intend
wes_block     <- omixs_raw[, 98:,  drop = FALSE]

common_samps2 <- intersect(rownames(drug_block), rownames(meta))
drug_block    <- drug_block[common_samps2, , drop = FALSE]
rna_block     <- rna_block[common_samps2, , drop = FALSE]
pathway_block <- pathway_block[common_samps2, , drop = FALSE]
wes_block     <- wes_block[common_samps2, , drop = FALSE]
group2        <- factor(meta[common_samps2, "group"])

data_list <- list(
  Drug       = drug_block,
  Expression = rna_block,
  Pathway    = pathway_block,
  WES        = wes_block
)

design <- matrix(0.1, nrow = length(data_list), ncol = length(data_list))
diag(design) <- 0
colnames(design) <- rownames(design) <- names(data_list)

# Fit DIABLO
diablo_model <- block.splsda(X = data_list, Y = group2, design = design, ncomp = 2)


plotIndiv(diablo_model, group = group2, legend = TRUE, title = "DIABLO: all blocks",
          ellipse = TRUE, ellipse.level = 0.95)
plotVar(diablo_model, comp = 1)  # variables contributing to comp 1 (all blocks)


keepX_list <- list(
  Expression = 15,
  Drug       = 10,
  Pathway    = 15,
  WES        = 15
)
diablo_model_top <- block.splsda(X = data_list, Y = group2, design = design, ncomp = 2, keepX = keepX_list)


set.seed(123)
perf_result <- perf(diablo_model_top, validation = "Mfold", folds = 5, nrepeat = 5)
plot(perf_result, criterion = "BER")


plotArrow(diablo_model_top, group = group2, legend = TRUE, title = "DIABLO arrow plot")
network(diablo_model_top, cutoff = 0.7)        # network of correlated features
circosPlot(diablo_model_top, comp = 1, cutoff = 0.7, showIntraLinks = TRUE)
cimDiablo(diablo_model_top, comp = 1, legend = TRUE, cluster = "both", title = "CIM for DIABLO")


load_drug <- diablo_model_top$loadings$Drug
load_rna  <- diablo_model_top$loadings$Expression
load_path <- diablo_model_top$loadings$Pathway
load_wes  <- diablo_model_top$loadings$WES

thr <- 0.1
sel_drug <- load_drug[abs(load_drug[,1]) > thr, , drop = FALSE]
sel_rna  <- load_rna [abs(load_rna [,1]) > thr, , drop = FALSE]
sel_path <- load_path[abs(load_path[,1]) > thr, , drop = FALSE]
sel_wes  <- load_wes [abs(load_wes [,1]) > thr, , drop = FALSE]

plotLoadings(diablo_model_top, block = "Drug",       comp = 1, contrib = "max")
plotLoadings(diablo_model_top, block = "Expression", comp = 1, contrib = "max")
plotLoadings(diablo_model_top, block = "Pathway",    comp = 1, contrib = "max")
plotLoadings(diablo_model_top, block = "WES",        comp = 1, contrib = "max")




if (!is.null(meta$PDS)) {
  pds <- factor(meta[common_samps2, "PDS"])
  plotIndiv(diablo_model_top, group = pds, legend = TRUE, ellipse = TRUE, ellipse.level = 0.95,
            title = "DIABLO (colored by PDS)")
}

if (!is.null(meta$CMS)) {
  cms <- factor(meta[common_samps2, "CMS"])
  plotIndiv(diablo_model_top, group = cms, legend = TRUE, ellipse = TRUE, ellipse.level = 0.95,
            title = "DIABLO (colored by CMS)")
}
