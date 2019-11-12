# Code

## load_imputed_data.R

This script loads the imputed data, removes selects only numeric columns, and removes all columns with 0 varaiation. It loads two dataframes into memory: 
 - `imputed_df`, the imputed metabolomics data
 - `keys_df`, the patient identifier codes