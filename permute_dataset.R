# Function to create the permuted datasets for PRONIA analyses.

library(tidyverse)

create_permutation <- function(df, num_perm){
  interview_num <- df %>% 
    mutate(id = row_number())%>%
    dplyr::select(id)%>%
    remove_rownames()
  set.seed(123)
  # Create shuffled row numbers for each permutation
  shuffled_row_numbers <- map(1:(num_perm-1), ~{
    #set.seed(123)  # Set seed for reproducibility within each iteration
    sample(interview_num$id)
  })
  
  # Bind original row_numbers with shuffled row numbers
  result <- bind_cols((shuffled_row_numbers))%>%
    cbind(interview_num, .)
  
  # Rename the shuffled columns
  names(result)[-1] <- paste0("shuffled_", 1:(num_perm-1))
  
  return(result)
}

# ROP-CHR: create permutation:
df_rop_chr_sim <- read.csv('simulated_data/data_for_permutation.csv')
perm_rop_chr_sim <- create_permutation(df_rop_chr_sim, 1000)
write.csv(perm_rop_chr_sim, 'simulated_data/perm_rop_chr.csv', row.names = FALSE)
