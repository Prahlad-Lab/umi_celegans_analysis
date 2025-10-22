# Title: Comparative Distribution and Violin Plots of Allele Proportions
# Description: This script reads two sets of files.
#              1. per_position_results.tsv: For a histogram/density plot showing
#                 the distribution of allele proportions at each genomic site.
#              2. summary_totals.txt: For a violin plot comparing the overall
#                 ratio of alternative reads for each sample replicate.

# --- 1. Install and Load Required Libraries ---
# If you haven't installed these packages, uncomment the lines below and run them once.
# install.packages("ggplot2")
# install.packages("readr")
# install.packages("scales")
# install.packages("dplyr")
install.packages("ggpubr") # For adding statistical tests to plots

library(ggplot2)
library(readr)
library(scales)
library(dplyr)
library(ggpubr)
library(cowplot)# For stat_compare_means()
setwd("/vscratch/grp-vprahlad/umi_analysis/data/output_broad/Family_size1_allele_proportions_depth10")
# --- 2. Define Input File Paths ---
# Define the directory where the files are located.
data_directory <- "/vscratch/grp-vprahlad/umi_analysis/data/output_broad/Family_size1_allele_proportions_depth10"

# --- PART A: Files for Histogram/Density Plot ---
# IMPORTANT: Please verify that these filenames are correct.
files_n2_per_position <- file.path(data_directory, c(
  "N2.30min.HS.1.singletons.per_position_results.tsv",
  "N2.30min.HS.2.singletons.per_position_results.tsv",
  "N2.30min.HS.3.singletons.per_position_results.tsv"
))

files_prde1_per_position <- file.path(data_directory, c(
  "PRDE1.30min.HS.1.singletons.per_position_results.tsv",
  "PRDE1.30min.HS.2.singletons.per_position_results.tsv",
  "PRDE1.30min.HS.3.singletons.per_position_results.tsv"
))

# --- PART B: Files for Violin Plot ---
# IMPORTANT: Please verify that these filenames are correct.
files_n2_summary <- file.path(data_directory, c(
  "N2.30min.HS.1.singletons.summary_totals.txt",
  "N2.30min.HS.2.singletons.summary_totals.txt",
  "N2.30min.HS.3.singletons.summary_totals.txt"
))

files_prde1_summary <- file.path(data_directory, c(
  "PRDE1.30min.HS.1.singletons.summary_totals.txt",
  "PRDE1.30min.HS.2.singletons.summary_totals.txt",
  "PRDE1.30min.HS.3.singletons.summary_totals.txt"
))


# --- 3. Function to Read and Combine Data ---
read_and_label <- function(files, group_label) {
  df <- bind_rows(lapply(files, function(file_path) {
    if (!file.exists(file_path)) {
      warning(paste("Warning: File not found, skipping:", file_path))
      return(NULL)
    }
    read_tsv(file_path, col_types = cols())
  }))
  
  if (nrow(df) > 0) {
    df$SampleGroup <- group_label
  }
  return(df)
}

# --- 4. Read Data for Both Plot Types ---
# Data for Histogram/Density plot
data_n2_dist <- read_and_label(files_n2_per_position, "N2")
data_prde1_dist <- read_and_label(files_prde1_per_position, "PRDE1")
dist_plot_data <- bind_rows(data_n2_dist, data_prde1_dist)

# Data for Violin plot
data_n2_violin <- read_and_label(files_n2_summary, "N2")
data_prde1_violin <- read_and_label(files_prde1_summary, "PRDE1")
violin_plot_data <- bind_rows(data_n2_violin, data_prde1_violin)

print("Data loaded and combined successfully for both plot types.")

dist_plot_data_2 =dist_plot_data %>% filter(Proportion>0)

dist_plot_data_3 =dist_plot_data %>% filter(Proportion==0)

# 
# dist_plot_data_3 %>% group_by(SampleGroup) %>% summarise(total = n())
# dist_plot_data_2 %>% group_by(SampleGroup) %>% summarise(total = n())
# 
# dist_plot_data %>% group_by(SampleGroup) %>% summarise(total = n())

# # --- 5. Create and Save the Comparative Histogram Plot ---
# if (nrow(dist_plot_data) > 0) {
#   print(paste("Total rows for histogram plot:", nrow(dist_plot_data)))
#   distribution_plot <- ggplot(dist_plot_data, aes(x = Proportion, fill = SampleGroup)) +
#     # Use geom_histogram with default y-axis (counts)
#     # position="identity" allows the histograms to overlap.
#     geom_histogram(aes(y =count/sum(count)),binwidth = 0.05, alpha = 0.3, position = "identity") +
#     scale_fill_manual(values = c("N2" = "black", "PRDE1" = "magenta")) +
#     #scale_color_manual(values = c("N2" = "#0072B2", "PRDE1" = "#D55E00")) +
#     scale_x_continuous(name = "Proportion of Alternative Reads per Site", breaks = seq(0, 1, by = 0.1), limits = c(0, 1.01)) +
#     scale_y_continuous(name = "Number of Sites (Frequency)") + # Format y-axis for readability
#     labs(
#       fill = "Sample Group"
#     ) +
#     theme_bw(base_size   = 14) +
#     theme(
#       plot.title = element_text(hjust = 0.5, face = "bold"),
#       plot.subtitle = element_text(hjust = 0.5, color = "gray30"),
#       panel.grid.minor = element_blank(), legend.position = "bottom"
#     )
#   
#   print(distribution_plot)
#   output_filename_dist <- "comparative_histogram_plot_with_zeros.svg"
#   ggsave(output_filename_dist, plot = distribution_plot, width = 8, height = 4, dpi = 300)
#   print(paste("Histogram plot saved as:", output_filename_dist))
# } else {
#   print("Error: No data found for the histogram plot.")
# }


library(dplyr)
library(tidyverse)

# --- 6. Create and Save the Violin Plot ---
if (nrow(
  violin_plot_data) > 0) {
  print(paste("Total rows for violin plot:", nrow(violin_plot_data)))
  violin_plot <- ggplot(violin_plot_data, aes(x = SampleGroup, y = Ratio_Alt_vs_Total, fill = SampleGroup)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.3, fill = "white", outlier.shape = NA) +
    geom_jitter(shape = 16, position = position_jitter(0.1), alpha = 0.8) + # Add points for each replicate
    scale_fill_manual(values = c("N2" = "#0072B2", "PRDE1" = "#D55E00")) +
    stat_compare_means(method = "wilcox.test", comparisons = list(c("N2", "PRDE1")), label = "p.format", bracket.size = 0.5) +
    scale_y_continuous(name = "Ratio of Alternative vs. Total Reads") +
    scale_x_discrete(name = "Sample Group") +
    labs(
      fill = "Sample Group"
    ) +
    theme_cowplot(font_size  = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "gray30"),
      legend.position = "none"
    )
  
  print(violin_plot)
  
  output_filename_violin <- "comparative_violin_plot_summary_data.svg"
  ggsave(output_filename_violin, plot = violin_plot, width = 4, height = 3.5, dpi = 300)
  print(paste("Violin plot saved as:", output_filename_violin))
} else {
  print("Error: No data found for the violin plot.")
}



# 1. Generate or load your data
set.seed(123)
data_values <- dist_plot_data %>% filter(Proportion<1&Proportion>0,SampleGroup=="N2")

# dist_plot_data %>%  filter(Proportion<1&Proportion>0,SampleGroup=="N2") %>% summarise(n())
# dist_plot_data %>% filter(Proportion<1&Proportion>0,SampleGroup=="PRDE1") %>% summarise(n())


total_n2 = dist_plot_data %>%  filter(Proportion<1,SampleGroup=="N2") %>% summarise(total =n()) %>% pull(total)
total_prde1 =  dist_plot_data %>%  filter(Proportion<1,SampleGroup=="PRDE1") %>% summarise(total =n()) %>% pull(total)

# dist_plot_data %>%  filter(Proportion<1,SampleGroup=="N2") %>% summarise(total =n())
# dist_plot_data %>% filter(Proportion<1,SampleGroup=="PRDE1") %>% summarise(total= n())

# non_zero
total_n2_variant_sites = dist_plot_data %>%  filter(Proportion<1&Proportion>0,SampleGroup=="N2") %>% summarise(total = n()) %>% pull(total)
total_prde1_variant_sites =dist_plot_data %>%  filter(Proportion<1&Proportion>0,SampleGroup=="PRDE1") %>% summarise(total = n()) %>% pull(total)

variant_sites_proportion_N2 =total_n2_variant_sites/total_n2
variant_sites_proportion_PRDE1 = tolat_prde1_variant_sites/total_prde1

# zero
total_n2_nonvariant_sites  =   dist_plot_data %>%  filter(Proportion==0,SampleGroup=="N2") %>% summarise(total =n())%>% pull(total)
total_prde1_nonvariant_sites  =   dist_plot_data %>%  filter(Proportion==0,SampleGroup=="PRDE1") %>% summarise(total =n())%>% pull(total)
non_variant_sites_proportion_N2 = total_n2_nonvariant_sites/total_n2
non_variant_sites_proportion_PRDE1 = total_prde1_nonvariant_sites/total_prde1

table.number.sites.proportion.variant.sites = data.frame(SampleGroup = c("N2","PRDE1"),Total_Number_Sites = c(total_n2,total_prde1),Site_w_variants = c(total_n2_variant_sites,total_prde1_variant_sites),Proportion_Variant_Sites= c(variant_sites_proportion_N2,variant_sites_proportion_PRDE1))
write.csv(table.number.sites.proportion.variant.sites,file = "table.number.sites.proportion.variant.sites.csv")

data_values_N2 <- dist_plot_data %>% filter(Proportion<1&Proportion>0,SampleGroup=="N2")

# 2. Define the breaks for your bins
breaks <- seq(0,1,by =0.05)



# 3. Create the binned categories
binned_data <- cut(data_values_N2$Proportion, breaks = breaks, include.lowest = TRUE, right = TRUE)

# 4. Generate the frequency table
frequency_table <- as.data.frame(table(binned_data))
colnames(frequency_table) <- c("Bin", "Frequency")
frequency_table_data =frequency_table %>% mutate(ncount = Frequency/sum(Frequency)) 

print(frequency_table_data)


data_values_p <- dist_plot_data %>% filter(Proportion<1&Proportion>0,SampleGroup=="PRDE1")
binned_data_p <- cut(data_values_p$Proportion, breaks = breaks, include.lowest = TRUE, right = TRUE)
frequency_table_p <- as.data.frame(table(binned_data_p))
colnames(frequency_table_p) <- c("Bin", "Frequency")
frequency_table_p_data = frequency_table_p %>% mutate(ncount = Frequency/sum(Frequency)) 
print(frequency_table_p_data)

# plot(x = frequency_table_p_data$Bin, y =frequency_table_p_data$ncount,type ="h")
# barplot(height  = frequency_table_p_data$Bin)

# library(patchwork)
# install.packages("patchwork")
# p1 =frequency_table_p_data %>% ggplot(aes(x = Bin,y =ncount ))+geom_col(fill = "magenta" )+theme_bw()
# p2 =frequency_table_data %>% ggplot(aes(x = Bin,y =ncount ))+geom_col(fill="gray")+theme_bw()
# p1/p2+plot_layout()

# Proportion of Variant Sites



colnames(frequency_table_p_data) = c("Bins","Sites Counts PRDE1", "Percentage PRDE1")

frequency_table_PRDE1_data_variant_sites = frequency_table_p_data %>% mutate(`Percentage PRDE1` = `Percentage PRDE1`*100)

colnames(frequency_table_data) = c("Bins","Sites Counts N2", "Percentage N2")
frequency_table_N2_data_variant_sites = frequency_table_data %>% mutate(`Percentage N2` = `Percentage N2`*100)


frequency_table_data_variant_sites = cbind.data.frame(frequency_table_N2_data_variant_sites,frequency_table_PRDE1_data_variant_sites[,-1])
write.csv(frequency_table_data_variant_sites, file = "frequency_table_data_variant_sites.csv")


Total_Error_Ratios = violin_plot_data %>% mutate(Replicate = rep(1:3,2)) %>% relocate(SampleGroup,Replicate)
write.csv(violin_plot_data, file = "Total_Error_Ratios_N2_PRDE1.csv")



dist_plot_data %>%  filter(Proportion<1,SampleGroup=="N2") %>% summarise(total =n())
dist_plot_data %>% filter(Proportion<1,SampleGroup=="PRDE1") %>% summarise(total= n())
#

data_values_N2 <- dist_plot_data %>% filter(Proportion<1,SampleGroup=="N2")

# 2. Define the breaks for your bins
breaks <- seq(0,1,by =0.05)



# 3. Create the binned categories
binned_data <- cut(data_values_N2$Proportion, breaks = breaks, include.lowest = TRUE, right = TRUE)

# 4. Generate the frequency table
frequency_table <- as.data.frame(table(binned_data))
colnames(frequency_table) <- c("Bin", "Frequency")
frequency_table_data =frequency_table %>% mutate(ncount = Frequency/sum(Frequency)) 

print(frequency_table_data)


data_values_p <- dist_plot_data %>% filter(Proportion<1,SampleGroup=="PRDE1")
binned_data_p <- cut(data_values_p$Proportion, breaks = breaks, include.lowest = TRUE, right = TRUE)
frequency_table_p <- as.data.frame(table(binned_data_p))
colnames(frequency_table_p) <- c("Bin", "Frequency")
frequency_table_p_data = frequency_table_p %>% mutate(ncount = Frequency/sum(Frequency)) 
print(frequency_table_p_data)



colnames(frequency_table_p_data) = c("Bins","Sites Counts PRDE1", "Percentage PRDE1")

frequency_table_PRDE1_data_all_sites = frequency_table_p_data %>% mutate(`Percentage PRDE1` = `Percentage PRDE1`*100)

colnames(frequency_table_data) = c("Bins","Sites Counts N2", "Percentage N2")
frequency_table_N2_data_all_sites = frequency_table_data %>% mutate(`Percentage N2` = `Percentage N2`*100)


frequency_table_data_all_sites = cbind.data.frame(frequency_table_N2_data_all_sites,frequency_table_PRDE1_data_all_sites[,-1])
write.csv(frequency_table_data_all_sites, file = "frequency_table_data_all_sites.csv")


