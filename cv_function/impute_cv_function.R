# ---------------------------------------------------------------------------- #
#                     Impute data within cross-validation                      #
# ---------------------------------------------------------------------------- #

# Function to perform imputation on training and test datasets within a cross-validation framework
# The function imputes missing values based on group-level statistics (mean for most variables and median for a specified variable)
impute_cv <- function(vars_train_df, vars_test_df, group_var, median_var, which_imputation){
  # Extract the names of the columns in the training data frame for later selection
  vars_to_select <- vars_train_df %>%
    as.data.frame(.)%>%
    names
  # Imputation for the training data
  if (which_imputation == 'train'){
    # Convert all variables to numeric and apply group-based imputation
    imputed_df <- vars_train_df %>%
      as.data.frame(.)%>%
      mutate(across(everything(), ~ as.numeric(.)))%>%
      group_by(get(group_var))%>%
      
      # Impute missing values for all columns except median_var using group-wise mean
      mutate(across(-(c(median_var)), 
                    ~ case_when(is.na(.) ~ mean(., na.rm = TRUE),
                                TRUE ~ .)))%>%
      
      # Impute missing values for median_var using group-wise median
      mutate(across(median_var, ~ case_when(is.na(.) ~ median(., na.rm = TRUE),
                                            TRUE ~ .)))%>%
      ungroup()%>%
      
      # Ensure columns are in the original order
      dplyr::select(vars_to_select)%>%
      as.matrix(.)
    
    # Imputation for the test data
  }else if(which_imputation == 'test'){
    # Calculate mean imputation values for all columns except median_var based on training data
    imputed_values_mean <- vars_train_df %>%
      as.data.frame(.)%>%
      mutate(across(everything(), ~as.numeric(.)))%>%
      group_by(get(group_var))%>%
      summarise(across(-c(median_var), ~mean(., na.rm = TRUE)))%>%
      dplyr::select(-`get(group_var)`)
    
    # Calculate median imputation values for median_var based on training data
    imputed_values_median <- vars_train_df %>%
      as.data.frame(.)%>%
      mutate(across(everything(), ~as.numeric(.)))%>%
      group_by(get(group_var))%>%
      summarise(across(median_var, ~median(., na.rm = TRUE)))%>%
      dplyr::rename(!!group_var := `get(group_var)`)
    
    # Apply the imputed values to the test data
    imputed_df <- vars_test_df %>%
      as.data.frame() %>%
      mutate(across(everything(), ~ as.numeric(.))) %>%
      # Left join to get the imputed mean values from the training data
      left_join(imputed_values_mean, by = group_var, suffix = c("", "_mean")) %>%
      # Left join to get the imputed median values from the training data
      left_join(imputed_values_median, by = group_var, suffix = c("", "_median")) %>%
      group_by(get(group_var)) %>%
      
      # Impute missing values in the test data using the pre-calculated mean values
      mutate(across(-(c(median_var, group_var, paste0(median_var, "_median"), ends_with('_mean'))), ~ case_when(
        is.na(.) ~ get(paste0(cur_column(), "_mean")),  # Use the pre-calculated mean
        TRUE ~ .
      ))) %>%
      
      # Impute missing values for median_var using the pre-calculated median values
      mutate(across(median_var, ~ case_when(
        is.na(.) ~ get(paste0(cur_column(), "_median")),  # Use the pre-calculated median
        TRUE ~ .
      ))) %>%
      ungroup() %>%
      
      # Ensure columns are in the original order
      dplyr::select(vars_to_select)%>%
      as.matrix(.)
  }
  
  # Return the imputed dataset either training or test data
  return(imputed_df)
}
