# ============================================================
# MULTISEQ_SNPplot FINAL SCRIPT WITH CLICKABLE SUMMARY LINKS
# AND OUTPUT FOLDER NAMED WITH _Rplots
# ============================================================

library(ggplot2)
library(openxlsx)
library(dplyr)

# ==== USER SETUP ====
setwd("/Users/jerzykulski/Desktop")
input_file <- "39v47.59.72.89.92.snps.txt"
bin_width <- 1000

# ==== CREATE OUTPUT FOLDER BASED ON INPUT FILE + "_Rplots" ====
input_name <- tools::file_path_sans_ext(basename(input_file))  # remove .txt
output_dir <- file.path(getwd(), paste0(input_name, "_Rplots"))  # add "_Rplots"
if (!dir.exists(output_dir)) dir.create(output_dir)
cat("Outputs will be saved in folder:", output_dir, "\n\n")

# ==== LOAD FILE ====
snp_data <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(snp_data) <- make.names(colnames(snp_data))
snp_data <- snp_data[!grepl("N", snp_data$SNP.pattern), ]
snp_data$SNP_position <- as.numeric(snp_data$sequence_1_PosInContg)
if (!is.character(snp_data$SNP.pattern)) snp_data$SNP.pattern <- as.character(snp_data$SNP.pattern)

# ==== PARSE SNP PATTERN ====
num_seqs <- nchar(snp_data$SNP.pattern[1])
pattern_split <- do.call(rbind, strsplit(snp_data$SNP.pattern, split = ""))
colnames(pattern_split) <- paste0("seq", 1:num_seqs)
snp_data <- cbind(snp_data, pattern_split)

# ==== DETECT SEQUENCE LABELS ====
contig_cols <- grep("sequence_\\d+_Contig", colnames(snp_data), value = TRUE)
sequence_labels <- sapply(contig_cols, function(col) {
  values <- unique(snp_data[[col]])
  values <- values[values != "null"]
  if (length(values) == 0) return(NA)
  sub("\\..*$", "", values[1])
})
names(sequence_labels) <- paste0("seq", seq_along(sequence_labels))

# ==== AUTO-DETECT REFERENCE ====
first_ref <- unique(snp_data$sequence_1_Contig)
first_ref <- first_ref[first_ref != "null"][1]
ref_seq_id <- sub("\\..*$", "", first_ref)
cat("\nReference sequence detected:", ref_seq_id, "\n")

# ==== IDENTIFY REFERENCE INDEX ====
ref_index <- which(sequence_labels == ref_seq_id)
if (length(ref_index) == 0) ref_index <- 1L
reference_seq <- paste0("seq", ref_index)

# ==== LOOP THROUGH OTHER SEQUENCES + SUMMARY ====
target_seqs <- setdiff(names(sequence_labels), reference_seq)
summary_list <- list()

for (target_seq in target_seqs) {
  label <- paste0(sequence_labels[[reference_seq]], "v", sequence_labels[[target_seq]])
  cat("Processing:", label, "\n")
  
  pairwise_df <- snp_data %>%
    filter(.data[[reference_seq]] != .data[[target_seq]]) %>%
    filter(.data[[reference_seq]] != "-" & .data[[target_seq]] != "-") %>%
    select(SNP_position)
  
  if (nrow(pairwise_df) == 0 || all(is.na(pairwise_df$SNP_position))) {
    warning(paste("No SNPs found for", label, "- skipping."))
    next
  }
  
  # --- Binning SNPs ---
  range_min <- min(pairwise_df$SNP_position, na.rm = TRUE)
  range_max <- max(pairwise_df$SNP_position, na.rm = TRUE)
  if (range_min == range_max) bins <- c(range_min, range_min + bin_width)
  else {
    start <- floor(range_min / bin_width) * bin_width
    end   <- ceiling((range_max + 1) / bin_width) * bin_width
    bins  <- seq(start, end, by = bin_width)
  }
  
  binned_snps <- data.frame(Bin = bins[-length(bins)], SNPs_kb = 0)
  for (i in seq_along(bins)[-length(bins)]) {
    snps_in_bin <- sum(pairwise_df$SNP_position >= bins[i] & pairwise_df$SNP_position < bins[i + 1])
    binned_snps$SNPs_kb[i] <- snps_in_bin / (bin_width / 1000)
  }
  
  avg_snp_density <- nrow(pairwise_df) / (max(pairwise_df$SNP_position) - min(pairwise_df$SNP_position)) * 1000
  y_max <- max(binned_snps$SNPs_kb)
  
  # --- File paths ---
  plot_path <- file.path(output_dir, paste0("snp_density_", label, ".png"))
  txt_path  <- file.path(output_dir, paste0("snp_bins_", label, ".txt"))
  xlsx_path <- file.path(output_dir, paste0("snp_bins_", label, ".xlsx"))
  
  # --- Save plot ---
  snp_plot <- ggplot(pairwise_df, aes(x = SNP_position)) +
    geom_histogram(aes(y = after_stat(count) / (bin_width / 1000)),
                   binwidth = bin_width, fill = "black", alpha = 0.6, color = "black") +
    geom_hline(yintercept = avg_snp_density, color = "red", linetype = "dashed", linewidth = 1) +
    annotate("text", x = mean(pairwise_df$SNP_position), y = y_max * 0.95,
             label = sprintf("%.2f SNPs/kb (%d SNPs)", avg_snp_density, nrow(pairwise_df)),
             color = "red", size = 4, fontface = "bold") +
    labs(title = paste("SNP Density:", label),
         x = "SNP Positions (bp)", y = "SNPs per kb") +
    theme_minimal()
  ggsave(filename = plot_path, plot = snp_plot, width = 10, height = 4, dpi = 300)
  
  # --- Save TXT and XLSX ---
  write.table(binned_snps, file = txt_path, sep = "\t", row.names = FALSE)
  wb <- createWorkbook()
  addWorksheet(wb, "SNPs per kb")
  writeData(wb, "SNPs per kb", binned_snps)
  saveWorkbook(wb, xlsx_path, overwrite = TRUE)
  
  # --- Add to summary with clickable links ---
  summary_list[[label]] <- data.frame(
    Reference = sequence_labels[[reference_seq]],
    Target = sequence_labels[[target_seq]],
    Average_SNPs_per_kb = avg_snp_density,
    Total_SNPs = nrow(pairwise_df),
    Plot_Link = sprintf('HYPERLINK("%s","%s")', plot_path, basename(plot_path)),
    TXT_Link  = sprintf('HYPERLINK("%s","%s")', txt_path, basename(txt_path)),
    XLSX_Link = sprintf('HYPERLINK("%s","%s")', xlsx_path, basename(xlsx_path)),
    stringsAsFactors = FALSE
  )
  cat("Done:", label, "\n")
}

# ==== WRITE SUMMARY TABLE ====
if (length(summary_list) > 0) {
  summary_df <- do.call(rbind, summary_list)
  summary_file <- file.path(output_dir, "pairwise_snp_density_summary.xlsx")
  
  wb_sum <- createWorkbook()
  addWorksheet(wb_sum, "Summary")
  writeData(wb_sum, "Summary", summary_df)
  
  # Convert columns to Excel hyperlinks
  for (col in c("Plot_Link", "TXT_Link", "XLSX_Link")) {
    writeFormula(wb_sum, sheet = "Summary",
                 x = summary_df[[col]],
                 startCol = which(colnames(summary_df) == col),
                 startRow = 2)
  }
  
  saveWorkbook(wb_sum, summary_file, overwrite = TRUE)
  cat("\nClickable summary table saved as", summary_file, "\n")
} else {
  cat("\nNo pairwise SNPs found; summary table not created.\n")
}

cat("\nAll pairwise SNP density plots and tables completed.\n")

