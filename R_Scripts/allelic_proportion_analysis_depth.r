# ... existing code ...
# Define the directory where the files are located.
data_directory <- "../Family_size1_allele_proportions_depth10"
ouput_dir = "../output/Allelic_proportion_depth_10_results"

# --- Create Output Directory ---
# Create the directory if it doesn't exist
dir.create(ouput_dir, recursive = TRUE, showWarnings = FALSE)
cat("All outputs will be saved to:", ouput_dir, "\n")

# --- PART A: Files for Histogram/Density Plot ---
# ... existing code ...
#     )
#   
#   print(distribution_plot)
#   output_filename_dist <- file.path(ouput_dir, "comparative_histogram_plot_with_zeros.svg")
#   ggsave(output_filename_dist, plot = distribution_plot, width = 8, height = 4, dpi = 300)
#   print(paste("Histogram plot saved as:", output_filename_dist))
# } else {
# ... existing code ...
print(violin_plot)

output_filename_violin <- file.path(ouput_dir, "comparative_violin_plot_summary_data.svg")
ggsave(output_filename_violin, plot = violin_plot, width = 4, height = 3.5, dpi = 300)
print(paste("Violin plot saved as:", output_filename_violin))
} else {
  print("Error: No data found for the violin plot.")
}



# 1. Generate or load your data
# set.seed(123) # Uncomment if you need reproducible jitter
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
variant_sites_proportion_PRDE1 = total_prde1_variant_sites/total_prde1 # Corrected typo from 'tolat_prde1_variant_sites'

# zero
total_n2_nonvariant_sites  =   dist_plot_data %>%  filter(Proportion==0,SampleGroup=="N2") %>% summarise(total =n())%>% pull(total)
total_prde1_nonvariant_sites  =   dist_plot_data %>%  filter(Proportion==0,SampleGroup=="PRDE1") %>% summarise(total =n())%>% pull(total)
non_variant_sites_proportion_N2 = total_n2_nonvariant_sites/total_n2
non_variant_sites_proportion_PRDE1 = total_prde1_nonvariant_sites/total_prde1

table.number.sites.proportion.variant.sites = data.frame(SampleGroup = c("N2","PRDE1"),Total_Number_Sites = c(total_n2,total_prde1),Site_w_variants = c(total_n2_variant_sites,total_prde1_variant_sites),Proportion_Variant_Sites= c(variant_sites_proportion_N2,variant_sites_proportion_PRDE1))
write.csv(table.number.sites.proportion.variant.sites,file = file.path(ouput_dir, "table.number.sites.proportion.variant.sites.csv"))

data_values_N2 <- dist_plot_data %>% filter(Proportion<1&Proportion>0,SampleGroup=="N2")
# ... existing code ...
frequency_table_N2_data_variant_sites = frequency_table_data %>% mutate(`Percentage N2` = `Percentage N2`*100)


frequency_table_data_variant_sites = cbind.data.frame(frequency_table_N2_data_variant_sites,frequency_table_PRDE1_data_variant_sites[,-1])
write.csv(frequency_table_data_variant_sites, file = file.path(ouput_dir, "frequency_table_data_variant_sites.csv"))


Total_Error_Ratios = violin_plot_data %>% mutate(Replicate = rep(1:3,2)) %>% relocate(SampleGroup,Replicate)
write.csv(Total_Error_Ratios, file = file.path(ouput_dir, "Total_Error_Ratios_N2_PRDE1.csv"))



dist_plot_data %>%  filter(Proportion<1,SampleGroup=="N2") %>% summarise(total =n())
# ... existing code ...
frequency_table_N2_data_all_sites = frequency_table_data %>% mutate(`Percentage N2` = `Percentage N2`*100)


frequency_table_data_all_sites = cbind.data.frame(frequency_table_N2_data_all_sites,frequency_table_PRDE1_data_all_sites[,-1])
write.csv(frequency_table_data_all_sites, file = file.path(ouput_dir, "frequency_table_data_all_sites.csv"))

