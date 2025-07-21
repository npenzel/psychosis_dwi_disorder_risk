# ---------------------------------------------------------------------------- #
#                       SIMULATE A RANDOM DATASET                              #
# ---------------------------------------------------------------------------- #


# Simulate a random dataset with similar characteristics like the PRONIA set
# to rerun the analysis:

# Load necessary libraries
library(dplyr)
library(tidyr)

# Set seed for reproducibility
set.seed(123)

# Number of subjects in each group
n_HC <- 298
n_CHR <- 183
n_ROP <- 206
n_ROD <- 195

set.seed(123) # For reproducibility

# Number of ROIs
rois_significant <- c("ACR", "ALIC", "SLF", "CGC", "CGH", "CP", "EC", "FX", 
                      "FX_ST", "GCC", "IFO", "SS", "PCR", "PLIC", "PTR", "RLIC", 
                      "SCC", "SCP")
rois_nonsignificant <- c("ML", "BCC", "ICP", "SFO", "SCR", "UNC", "CST")
rois <- c(rois_significant, rois_nonsignificant)

# Create a data frame for the simulation
sim_data <- data.frame(
  PSN = c(sample(1000:2000, n_HC, replace = FALSE),
          sample(3000:4000, n_CHR, replace = FALSE),
          sample(5000:6000, n_ROP, replace = FALSE),
          sample(7000:8000, n_ROD, replace = FALSE)),
  Studygroup = factor(rep(c("HC", "CHR", "ROP", "ROD"), times = c(n_HC, n_CHR, n_ROP, n_ROD)))
)

# Function to simulate FA values with group-level differences
simulate_fa_values <- function(group, roi_type) {
  if (roi_type == "significant") {
    # Significant regions will have noticeable differences between groups
    if (group == "HC") return(rnorm(1, mean = 0.70, sd = 0.03))
    if (group == "CHR") return(rnorm(1, mean = 0.66, sd = 0.03))
    if (group == "ROP") return(rnorm(1, mean = 0.63, sd = 0.03))
    if (group == "ROD") return(rnorm(1, mean = 0.64, sd = 0.03))
  } else {
    # Nonsignificant regions will have similar means between groups
    return(rnorm(1, mean = 0.65, sd = 0.03))
  }
}

# Simulate FA data for each individual and each ROI
fa_data <- data.frame(
  PSN = rep(sim_data$PSN, each = length(rois)),
  Studygroup = rep(sim_data$Studygroup, each = length(rois)),
  ROI = rep(rois, times = n_HC + n_CHR + n_ROP + n_ROD),
  FA_value = NA
)

# Fill in FA data with simulated values
for (i in seq_len(nrow(fa_data))) {
  roi_type <- ifelse(fa_data$ROI[i] %in% rois_significant, "significant", "nonsignificant")
  fa_data$FA_value[i] <- simulate_fa_values(fa_data$Studygroup[i], roi_type)
}

fa_data <- fa_data %>%
  pivot_wider(names_from = ROI, values_from = FA_value)

# Simulate clinical variables (as done previously)
sim_data$age <- c(
  rnorm(n_HC, mean = 25.4, sd = 6.1),
  rnorm(n_CHR, mean = 23.9, sd = 5.7),
  rnorm(n_ROP, mean = 25.9, sd = 5.8),
  rnorm(n_ROD, mean = 25.6, sd = 5.8)
)
sim_data$sex <- c(
  rbinom(n_HC, 1, 0.60),
  rbinom(n_CHR, 1, 0.52),
  rbinom(n_ROP, 1, 0.45),
  rbinom(n_ROD, 1, 0.48)
)
sim_data$wais_v <- c(
  rnorm(n_HC, mean = 11.7, sd = 2.7),
  rnorm(n_CHR, mean = 11.1, sd = 3.3),
  rnorm(n_ROP, mean = 9.5, sd = 3.4),
  rnorm(n_ROD, mean = 10.8, sd = 3)
)
sim_data$wais_mr <- c(
  rnorm(n_HC, mean = 11.4, sd = 2.2),
  rnorm(n_CHR, mean = 11, sd = 2.4),
  rnorm(n_ROP, mean = 10, sd = 2.7),
  rnorm(n_ROD, mean = 10.9, sd = 2.5)
)
sim_data$uhr_1stdegree <- c(
  rbinom(n_HC, 1, 0.0),
  rbinom(n_CHR, 1, 0.197),
  rbinom(n_ROP, 1, 0.155),
  rbinom(n_ROD, 1, 0.036)
)
sim_data$spia_sum_t0 <- c(
  rnorm(n_HC, mean = 0.1, sd = 0.3),
  rnorm(n_CHR, mean = 6.2, sd = 5.2),
  rnorm(n_ROP, mean = 7.0, sd = 7.3),
  rnorm(n_ROD, mean = 1.7, sd = 2.3)
)
sim_data$sips_p_t0 <- c(
  rnorm(n_HC, mean = 0.3, sd = 1),
  rnorm(n_CHR, mean = 8.1, sd = 4.7),
  rnorm(n_ROP, mean = 15.7, sd = 5.2),
  rnorm(n_ROD, mean = 2.0, sd = 2.1)
)
sim_data$sips_n_t0 <- c(
  rnorm(n_HC, mean = 0.2, sd = 0.9),
  rnorm(n_CHR, mean = 10.4, sd = 6.5),
  rnorm(n_ROP, mean = 10.8, sd = 7.5),
  rnorm(n_ROD, mean = 10.5, sd = 5.6)
)
sim_data$sips_g_t0 <- c(
  rnorm(n_HC, mean = 0.3, sd = 0.8),
  rnorm(n_CHR, mean = 8.1, sd = 4.0),
  rnorm(n_ROP, mean = 7.3, sd = 4.3),
  rnorm(n_ROD, mean = 8.3, sd = 3.5)
)
sim_data$sips_d_t0 <- c(
  rnorm(n_HC, mean = 0.1, sd = 0.4),
  rnorm(n_CHR, mean = 3.5, sd = 3.0),
  rnorm(n_ROP, mean = 5.2, sd = 4.1),
  rnorm(n_ROD, mean = 2.4, sd = 2.2)
)
sim_data$GAF_S_PastMonth_Screening <- c(
  rnorm(n_HC, mean = 85.7, sd = 7.1),
  rnorm(n_CHR, mean = 51.7, sd = 11.6),
  rnorm(n_ROP, mean = 43.1, sd = 13.6),
  rnorm(n_ROD, mean = 53, sd = 12)
)
sim_data$GAF_DI_PastMonth_Screening <- c(
  rnorm(n_HC, mean = 87.1, sd = 6.7),
  rnorm(n_CHR, mean = 79.2, sd = 7.5),
  rnorm(n_ROP, mean = 78.7, sd = 9.2),
  rnorm(n_ROD, mean = 81, sd = 8)
)
sim_data$GAF_S_LifeTime_Screening <- c(
  rnorm(n_HC, mean = 87.0, sd = 7.6),
  rnorm(n_CHR, mean = 51.4, sd = 10.5),
  rnorm(n_ROP, mean = 39.2, sd = 14.7),
  rnorm(n_ROD, mean = 82, sd = 8)
)
sim_data$GAF_DI_LifeTime_Screening <- c(
  rnorm(n_HC, mean = 88.4, sd = 6.9),
  rnorm(n_CHR, mean = 79.3, sd = 8.6),
  rnorm(n_ROP, mean = 79.2, sd = 9.8),
  rnorm(n_ROD, mean = 81, sd = 8)
)
sim_data$exposome_score <- c(
  rnorm(n_HC, mean = 0, sd = 1),
  rnorm(n_CHR, mean = 6.2, sd = 5.2),
  rnorm(n_ROP, mean = 7.0, sd = 7.3),
  rnorm(n_ROD, mean = 2, sd = 1)
)
sim_data$CPZE_cum_sum <- c(
  rnorm(n_HC, mean = 0, sd = 1),
  rnorm(n_CHR, mean = 6, sd = 5),
  rnorm(n_ROP, mean = 7, sd = 6),
  rnorm(n_ROD, mean = 578, sd = 2992)
)
sim_data$bdi_t0 <- c(
  rnorm(n_HC, mean = 3.4, sd = 5.2),
  rnorm(n_CHR, mean = 29.1, sd = 11.7),
  rnorm(n_ROP, mean = 21.6, sd = 12.3),
  rnorm(n_ROD, mean = 29, sd = 13)
)

# add also missing data:
# Function to apply missing data based on presence percentages
apply_missing <- function(data, group_size, percent_present) {
  present_idx <- sample(1:group_size, size = round(group_size * percent_present), replace = FALSE)
  missing_idx <- setdiff(1:group_size, present_idx)
  data[missing_idx] <- NA
  return(data)
}

# Apply missing data patterns for CHR and ROP

# WAIS vocabulary (92% overall, 91% ROP, 94% CHR, 96 % ROD)
sim_data$wais_v[sim_data$Studygroup == "ROP"] <- apply_missing(sim_data$wais_v[sim_data$Studygroup == "ROP"], n_ROP, 0.91)
sim_data$wais_v[sim_data$Studygroup == "CHR"] <- apply_missing(sim_data$wais_v[sim_data$Studygroup == "CHR"], n_CHR, 0.94)
sim_data$wais_v[sim_data$Studygroup == "ROD"] <- apply_missing(sim_data$wais_v[sim_data$Studygroup == "ROD"], n_ROD, 0.96)

# WAIS matrices (90% overall, 91% ROP, 90% CHR, 95 % ROD)
sim_data$wais_mr[sim_data$Studygroup == "ROP"] <- apply_missing(sim_data$wais_mr[sim_data$Studygroup == "ROP"], n_ROP, 0.91)
sim_data$wais_mr[sim_data$Studygroup == "CHR"] <- apply_missing(sim_data$wais_mr[sim_data$Studygroup == "CHR"], n_CHR, 0.90)
sim_data$wais_mr[sim_data$Studygroup == "ROD"] <- apply_missing(sim_data$wais_mr[sim_data$Studygroup == "ROD"], n_ROD, 0.95)

# Familial risk (97% overall, 98% ROP, 96% CHR)
sim_data$uhr_1stdegree[sim_data$Studygroup == "ROP"] <- apply_missing(sim_data$uhr_1stdegree[sim_data$Studygroup == "ROP"], n_ROP, 0.98)
sim_data$uhr_1stdegree[sim_data$Studygroup == "CHR"] <- apply_missing(sim_data$uhr_1stdegree[sim_data$Studygroup == "CHR"], n_CHR, 0.96)
sim_data$uhr_1stdegree[sim_data$Studygroup == "ROD"] <- apply_missing(sim_data$uhr_1stdegree[sim_data$Studygroup == "ROD"], n_ROD, 0.98)

# SPI-A symptoms (84% overall, 85% ROP, 82% CHR)
sim_data$spia_sum_t0[sim_data$Studygroup == "ROP"] <- apply_missing(sim_data$spia_sum_t0[sim_data$Studygroup == "ROP"], n_ROP, 0.85)
sim_data$spia_sum_t0[sim_data$Studygroup == "CHR"] <- apply_missing(sim_data$spia_sum_t0[sim_data$Studygroup == "CHR"], n_CHR, 0.82)
sim_data$spia_sum_t0[sim_data$Studygroup == "ROD"] <- apply_missing(sim_data$spia_sum_t0[sim_data$Studygroup == "ROD"], n_ROD, 0.89)

# SIPS positive (94% overall, 95% ROP, 95% CHR)
sim_data$sips_p_t0[sim_data$Studygroup == "ROP"] <- apply_missing(sim_data$sips_p_t0[sim_data$Studygroup == "ROP"], n_ROP, 0.95)
sim_data$sips_p_t0[sim_data$Studygroup == "CHR"] <- apply_missing(sim_data$sips_p_t0[sim_data$Studygroup == "CHR"], n_CHR, 0.95)
sim_data$sips_p_t0[sim_data$Studygroup == "ROD"] <- apply_missing(sim_data$sips_p_t0[sim_data$Studygroup == "ROD"], n_ROD, 0.97)

# SIPS negative (95% overall, 96% ROP, 94% CHR)
sim_data$sips_n_t0[sim_data$Studygroup == "ROP"] <- apply_missing(sim_data$sips_n_t0[sim_data$Studygroup == "ROP"], n_ROP, 0.96)
sim_data$sips_n_t0[sim_data$Studygroup == "CHR"] <- apply_missing(sim_data$sips_n_t0[sim_data$Studygroup == "CHR"], n_CHR, 0.94)
sim_data$sips_n_t0[sim_data$Studygroup == "ROD"] <- apply_missing(sim_data$sips_n_t0[sim_data$Studygroup == "ROD"], n_ROD, 0.98)

# SIPS general (94% overall, 93% ROP, 95% CHR)
sim_data$sips_g_t0[sim_data$Studygroup == "ROP"] <- apply_missing(sim_data$sips_g_t0[sim_data$Studygroup == "ROP"], n_ROP, 0.93)
sim_data$sips_g_t0[sim_data$Studygroup == "CHR"] <- apply_missing(sim_data$sips_g_t0[sim_data$Studygroup == "CHR"], n_CHR, 0.95)
sim_data$sips_g_t0[sim_data$Studygroup == "ROD"] <- apply_missing(sim_data$sips_g_t0[sim_data$Studygroup == "ROD"], n_ROD, 0.98)

# SIPS disorganized (95% overall, 95% ROP, 95% CHR)
sim_data$sips_d_t0[sim_data$Studygroup == "ROP"] <- apply_missing(sim_data$sips_d_t0[sim_data$Studygroup == "ROP"], n_ROP, 0.95)
sim_data$sips_d_t0[sim_data$Studygroup == "CHR"] <- apply_missing(sim_data$sips_d_t0[sim_data$Studygroup == "CHR"], n_CHR, 0.95)
sim_data$sips_d_t0[sim_data$Studygroup == "ROD"] <- apply_missing(sim_data$sips_d_t0[sim_data$Studygroup == "ROD"], n_ROD, 0.98)

# GAF-Symptoms: Past Month (97% overall, 98% ROP, 97% CHR)
sim_data$GAF_S_PastMonth_Screening[sim_data$Studygroup == "ROP"] <- apply_missing(sim_data$GAF_S_PastMonth_Screening[sim_data$Studygroup == "ROP"], n_ROP, 0.98)
sim_data$GAF_S_PastMonth_Screening[sim_data$Studygroup == "CHR"] <- apply_missing(sim_data$GAF_S_PastMonth_Screening[sim_data$Studygroup == "CHR"], n_CHR, 0.97)
sim_data$GAF_S_PastMonth_Screening[sim_data$Studygroup == "ROD"] <- apply_missing(sim_data$GAF_S_PastMonth_Screening[sim_data$Studygroup == "ROD"], n_ROD, 0.98)

# GAF-DI: Past Month (97% overall, 98% ROP, 97% CHR)
sim_data$GAF_DI_PastMonth_Screening[sim_data$Studygroup == "ROP"] <- apply_missing(sim_data$GAF_DI_PastMonth_Screening[sim_data$Studygroup == "ROP"], n_ROP, 0.98)
sim_data$GAF_DI_PastMonth_Screening[sim_data$Studygroup == "CHR"] <- apply_missing(sim_data$GAF_DI_PastMonth_Screening[sim_data$Studygroup == "CHR"], n_CHR, 0.97)
sim_data$GAF_DI_PastMonth_Screening[sim_data$Studygroup == "ROD"] <- apply_missing(sim_data$GAF_DI_PastMonth_Screening[sim_data$Studygroup == "ROD"], n_ROD, 0.98)

# GAF-Symptoms: Lifetime (97% overall, 98% ROP, 97% CHR)
sim_data$GAF_S_LifeTime_Screening[sim_data$Studygroup == "ROP"] <- apply_missing(sim_data$GAF_S_LifeTime_Screening[sim_data$Studygroup == "ROP"], n_ROP, 0.98)
sim_data$GAF_S_LifeTime_Screening[sim_data$Studygroup == "CHR"] <- apply_missing(sim_data$GAF_S_LifeTime_Screening[sim_data$Studygroup == "CHR"], n_CHR, 0.97)
sim_data$GAF_S_LifeTime_Screening[sim_data$Studygroup == "ROD"] <- apply_missing(sim_data$GAF_S_LifeTime_Screening[sim_data$Studygroup == "ROD"], n_ROD, 0.98)

# GAF-DI: Lifetime (97% overall, 98% ROP, 97% CHR)
sim_data$GAF_DI_LifeTime_Screening[sim_data$Studygroup == "ROP"] <- apply_missing(sim_data$GAF_DI_LifeTime_Screening[sim_data$Studygroup == "ROP"], n_ROP, 0.98)
sim_data$GAF_DI_LifeTime_Screening[sim_data$Studygroup == "CHR"] <- apply_missing(sim_data$GAF_DI_LifeTime_Screening[sim_data$Studygroup == "CHR"], n_CHR, 0.97)
sim_data$GAF_DI_LifeTime_Screening[sim_data$Studygroup == "ROD"] <- apply_missing(sim_data$GAF_DI_LifeTime_Screening[sim_data$Studygroup == "ROD"], n_ROD, 0.98)

# Exposome score (80% overall, 79% ROP, 78% CHR)
sim_data$exposome_score[sim_data$Studygroup == "ROP"] <- apply_missing(sim_data$exposome_score[sim_data$Studygroup == "ROP"], n_ROP, 0.79)
sim_data$exposome_score[sim_data$Studygroup == "CHR"] <- apply_missing(sim_data$exposome_score[sim_data$Studygroup == "CHR"], n_CHR, 0.78)
sim_data$exposome_score[sim_data$Studygroup == "ROD"] <- apply_missing(sim_data$exposome_score[sim_data$Studygroup == "ROD"], n_ROD, 0.74)

# Define the scan site distribution for each group
scan_site_distribution <- data.frame(
  Scan_site = c("BHAM1_3.2_TRBelow8001", "LMU_4.1.3", "LMU_Other_SoftwareVersion", "Turku",
                "UBS1.1_Prisma_D13D", "UBS1.2_Prisma_E11", "UBS2_Verio", "Udine2",
                "UKK1_Achieva_b01", "UKK1_Achieva_b02", "UKK2_Ingenia", "UMUENS"),
  CHR = c(11, 53, 41, 16, 4, 1, 18, 1, 6, 6, 9, 17),
  HC = c(28, 58, 3, 43, 11, 16, 44, 14, 34, 21, 12, 14),
  ROP = c(8, 71, 36, 27, 4, 0, 18, 1, 9, 10, 11, 11),
  ROD = c(9, 68, 34, 13, 1, 0, 18, 1, 4, 14, 20, 13)
)

# Helper function to assign scan sites based on group distribution
assign_scan_sites <- function(group_data, group_name) {
  # Filter scan site distribution for the current group
  group_distribution <- scan_site_distribution[, c("Scan_site", group_name)]
  
  # Initialize empty vector for scan sites
  scan_sites <- character(nrow(group_data))
  
  # For each scan site, assign the number of individuals in the group based on distribution
  start_idx <- 1
  for (i in 1:nrow(group_distribution)) {
    num_to_assign <- group_distribution[[group_name]][i]
    if (num_to_assign > 0) {
      end_idx <- start_idx + num_to_assign - 1
      scan_sites[start_idx:end_idx] <- group_distribution$Scan_site[i]
      start_idx <- end_idx + 1
    }
  }
  
  # Shuffle the scan sites for random assignment
  scan_sites <- sample(scan_sites)
  
  return(scan_sites)
}

# Apply the function to assign scan sites to each group in sim_data
sim_data$Scan_site[sim_data$Studygroup == "CHR"] <- assign_scan_sites(sim_data[sim_data$Studygroup == "CHR", ], "CHR")
sim_data$Scan_site[sim_data$Studygroup == "HC"] <- assign_scan_sites(sim_data[sim_data$Studygroup == "HC", ], "HC")
sim_data$Scan_site[sim_data$Studygroup == "ROP"] <- assign_scan_sites(sim_data[sim_data$Studygroup == "ROP", ], "ROP")
sim_data$Scan_site[sim_data$Studygroup == "ROD"] <- assign_scan_sites(sim_data[sim_data$Studygroup == "ROD", ], "ROD")

# Check the result
table(sim_data$Scan_site, sim_data$Studygroup)

# Merge clinical and FA data
final_data <- left_join(sim_data, fa_data, by = c("PSN", "Studygroup"))

write.csv(final_data, 'simulated_data/simulated_data.csv', row.names = FALSE)


# for permutation testing of the CCA you will have to create permuted datasets.
# Since CCA is performed on ROP and CHR only subset your dataframe

simulated_data_ropchr <- final_data %>% 
  filter(Studygroup %in% c('ROP', 'CHR'))%>%
  # You want to exclude everyone who is missing all variables included in CCA
  dplyr::select(-rois, -c(age, sex))%>%
  mutate(across(everything(), ~as.character(.)))%>%
  filter(!sum(is.na(c_across())) == ncol(.) - 2)%>%
  dplyr::select(PSN)

write.csv(simulated_data_ropchr, 'simulated_data/data_for_permutation.csv', row.names = FALSE)