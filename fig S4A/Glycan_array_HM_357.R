library(tidyverse)
rep1 <- read.csv("rep1.csv",header=TRUE,sep=",")
rep2 <- read.csv("rep2.csv",header=TRUE,sep=",")
rep3 <- read.csv("rep3.csv",header=TRUE,sep=",")

# Function to normalize by max value from columns 126-130
normalize_by_max_5cols <- function(data) {
  data %>%
    mutate(
      # Handle NA values in columns 126-130
      col_126 = ifelse(is.na(.[[126]]), 0, .[[126]]),
      col_127 = ifelse(is.na(.[[127]]), 0, .[[127]]),
      col_128 = ifelse(is.na(.[[128]]), 0, .[[128]]),
      col_129 = ifelse(is.na(.[[129]]), 0, .[[129]]),
      col_130 = ifelse(is.na(.[[130]]), 0, .[[130]]),
      # Find max absolute value across these 5 columns for each row
      max_val = abs(pmax(col_126, col_127, col_128, col_129, col_130))
    ) %>%
    # Subtract max_val from all data columns (4 to ncol-6) for ALL rows
    mutate(across(4:(ncol(.) - 6), ~ .x - max_val)) %>%
    # Remove temporary columns
    select(-col_126, -col_127, -col_128, -col_129, -col_130, -max_val)
}

# Apply normalization to all replicates
rep1_adj_maxAz <- normalize_by_max_5cols(rep1)
rep2_adj_maxAz <- normalize_by_max_5cols(rep2)
rep3_adj_maxAz <- normalize_by_max_5cols(rep3)


# Calculate average
avg <- rep1_adj_maxAz[, 1:3]
for(i in 4:ncol(rep1_adj_maxAz)) {
  avg[, i] <- (rep1_adj_maxAz[, i] + rep2_adj_maxAz[, i] + rep3_adj_maxAz[, i]) / 3
}
colnames(avg) <- colnames(rep1_adj_maxAz)
write.csv(avg, file = "avg_adj_max_357.csv", row.names = FALSE)

# Calculate standard deviation
sd <- rep1_adj_maxAz[, 1:3]
for(i in 4:ncol(rep1_adj_maxAz)) {
  values <- cbind(rep1_adj_maxAz[, i], rep2_adj_maxAz[, i], rep3_adj_maxAz[, i])
  sd[, i] <- apply(values, 1, sd)
}
colnames(sd) <- colnames(rep1_adj_maxAz)
write.csv(sd, file = "sd_adj_max_357.csv", row.names = FALSE)

avg_clean <- avg[rowSums(!is.na(avg[, 4:130])) > 0, ]
sd_clean <- sd[rowSums(!is.na(sd[, 4:130])) > 0, ]

avg_final <- avg_clean %>%
  mutate_at(vars(4:130), ~replace_na(., 0)) %>%
  select(1:125)

sd_final <- sd_clean %>%
  mutate_at(vars(4:130), ~ifelse(is.na(.) | . == 0, 1, .)) %>%
  select(1:125)

total <- avg_final
total[, 4:ncol(total)] <- avg_final[, 4:ncol(avg_final)] - sd_final[, 4:ncol(sd_final)]

write.csv(sd_final, file = "sd_adj_max_clean_357.csv", row.names = FALSE)
write.csv(avg_final, file = "avg_adj_max_clean_357.csv", row.names = FALSE)
write.csv(total, file = "total_adj_max_clean_357.csv", row.names = FALSE)

#Set QC function
apply_qc_corrections <- function(avg_df, sd_df, cv_threshold = 0.3, outlier_threshold = 3) {
  # 1. QC with CV
  cv_matrix <- sd_df[, 4:ncol(sd_df)] / (abs(avg_df[, 4:ncol(avg_df)]) + 0.001)  # Avoid to divide with 0
  # 2. QC with Outliers
  protein_outlier_scores <- apply(avg_df[, 4:ncol(avg_df)], 1, function(row) {
    valid_values <- row[!is.na(row)]
    if(length(valid_values) < 5) return(0)  # Skip with too small
    # Scoring with MAD values
    median_val <- median(valid_values)
    mad_val <- mad(valid_values)
    if(mad_val == 0) return(0)
    # Scoring the row values
    outlier_score <- median(abs(valid_values - median_val) / mad_val)
    return(outlier_score)
  })
    # 3. QC data generation
  corrected_avg <- avg_df
  # Penalty with higher CV 
  high_cv_mask <- cv_matrix > cv_threshold
  corrected_avg[, 4:ncol(corrected_avg)][high_cv_mask] <- -9999
  # Penalty with high outliers
  outlier_proteins <- which(protein_outlier_scores > outlier_threshold)
  if(length(outlier_proteins) > 0) {
    corrected_avg[outlier_proteins, 4:ncol(corrected_avg)] <- -9999
  }
  return(list(
    corrected_data = corrected_avg,
    cv_matrix = cv_matrix,
    protein_outlier_scores = protein_outlier_scores,
    flagged_proteins = outlier_proteins
  ))
}

# Calculate QC
qc_results <- apply_qc_corrections(total, sd_final, 
                                   cv_threshold = 0.6,  
                                   outlier_threshold = 0.7)   
total_qc <- qc_results$corrected_data

# Save QC dataframe
write.csv(total_qc, file = "total_qc_357.csv", row.names = FALSE)


quick_check2 <- function(gene_name, col_name) {
  gene_row <- which(df$Protein_name == gene_name)
  if(length(gene_row) == 0) {
    cat("Gene not found!\n")
    return(NULL)
  }
  gene_row <- gene_row[1]
  cat("Gene:", gene_name, "| Column:", col_name, "\n")
  cat("Value:", df[gene_row, col_name], "\n")
}

df <- total_qc

#X104: Chitin,  X114: Chitosan
quick_check2("LYK5", "X104")
quick_check2("LYK5", "X114")
quick_check2("CERK1/LYK1", "X104")

#X115: RG-I
quick_check2("ARM1a", "X115")
quick_check2("ARM1b", "X115")
quick_check2("ARM2a/MLRK/IMK3", "X115")
quick_check2("ARM3a", "X115")
quick_check2("ARM3b", "X115")

#X115: ARAn
quick_check2("ARM1a", "X117")
quick_check2("ARM1b", "X117")
quick_check2("ARM2a/MLRK/IMK3", "X117")
quick_check2("ARM3a", "X117")
quick_check2("ARM3b", "X117")

#X122: GALn
quick_check2("ARM1a", "X122")
quick_check2("ARM1b", "X122")
quick_check2("ARM2a/MLRK/IMK3", "X122")
quick_check2("ARM3a", "X122")
quick_check2("ARM3b", "X122")

threshold <- 7430

numeric_data_qc <- df[, 4:ncol(df)]
values_above_threshold_qc <- sum(numeric_data_qc > threshold, na.rm = TRUE)
cat("total_qc - Values larger than threshold:", values_above_threshold_qc, "\n")

high_values_qc <- which(numeric_data_qc > threshold, arr.ind = TRUE)
results_qc <- data.frame(
  Gene = df$Gene[high_values_qc[, "row"]],
  Protein_name = df$Protein_name[high_values_qc[, "row"]],
  Glycan = colnames(numeric_data_qc)[high_values_qc[, "col"]],
  Value = numeric_data_qc[high_values_qc],
  stringsAsFactors = FALSE
)
write.csv(results_qc, file = "result_qc_357.csv", row.names = FALSE)

numeric_data_total <- total[, 4:ncol(total)]
values_above_threshold_total <- sum(numeric_data_total > threshold, na.rm = TRUE)
cat("total - Values larger than threshold:", values_above_threshold_total, "\n")

high_values_total <- which(numeric_data_total > threshold, arr.ind = TRUE)
results_total <- data.frame(
  Gene = total$Gene[high_values_total[, "row"]],
  Protein_name = total$Protein_name[high_values_total[, "row"]],
  Glycan = colnames(numeric_data_total)[high_values_total[, "col"]],
  Value = numeric_data_total[high_values_total],
  stringsAsFactors = FALSE
)
write.csv(results_total, file = "result_total_357.csv", row.names = FALSE)

numeric_data_avg <- avg_final[, 4:ncol(avg_final)]
values_above_threshold_avg <- sum(numeric_data_avg > threshold, na.rm = TRUE)
cat("avg_final - Values larger than threshold:", values_above_threshold_avg, "\n")

high_values_avg <- which(numeric_data_avg > threshold, arr.ind = TRUE)
results_avg <- data.frame(
  Gene = avg_final$Gene[high_values_avg[, "row"]],
  Protein_name = avg_final$Protein_name[high_values_avg[, "row"]],
  Glycan = colnames(numeric_data_avg)[high_values_avg[, "col"]],
  Value = numeric_data_avg[high_values_avg],
  stringsAsFactors = FALSE
)



#### Generate Heatmap ####
library(pheatmap)
library(RColorBrewer)

heatmap_data <- as.matrix(df[, 4:ncol(df)])
rownames(heatmap_data) <- df$Gene
colnames(heatmap_data) <- gsub("^X", "#", colnames(heatmap_data))

heatmap_data_filtered <- heatmap_data
heatmap_data_filtered[heatmap_data_filtered < 0] <- 0

breaks <- c(seq(0, 3500, length.out = 21), 
            seq(3500, 20000, length.out = 71)[-1]) 

pdf("heatmap.pdf", width = 24, height = 32) 
pheatmap(heatmap_data_filtered, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         scale = "none",
         color = colorRampPalette(brewer.pal(9, "Purples"))(100),
         breaks = breaks,
         cellwidth = 8,
         cellheight = 6,
         border_color = NA,
         fontsize = 5,
)
dev.off()

