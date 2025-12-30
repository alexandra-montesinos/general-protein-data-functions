# general-protein-data-functions
Functions that clean &amp; plot protein data.

# clean_data.R
## find_regulation 
| Parameter | Type       | Description|
|-----------|------------|------------|
| dt  | data.table | finds regulation of proteins in this dt |
| dt2 | data.table | Optional. If included, will find if the regulation of the first datatable agrees with that of the second |
|ratio_cutoff| numeric | Default 0. If >0, will add a cutoff to what is considered significant up or down-regulation, marking those that don't make it as Not Significant|

## clean_dt
Used to standardize different styles of data naming styles. It is often necessary to use this function if you are calling my other datatable functions, so you'll have to use it even if you hate my particular naming style. Sorry! 
If your dt contains a gene column, use clean_gene_dt instead. It will call this function for you.
| Parameter | Type       | Description|
|-----------|------------|------------|
| dt  | data.table | dt to clean, |

## clean_gene_dt
Cleans a datatable with a Genes column by removing invalid data and standardizing names. Will also call clean_dt() on it.
WARNING: This function removes rows with NAs and NaNs in the Genes column. These missing values (especially NaNs) are sometimes caused because certain proteins in our sample don't have gene names. In these cases, dropping them like this will just remove valid data from the set. Of course, this should be a pretty rare occurence. Only ever happened to me once.
| Parameter | Type       | Description|
|-----------|------------|------------|
| dt  | data.table | dt to clean |
| remove_repeat_genes | boolean | Default FALSE. Removes repeat gnes, keeping the first instance. |
|add_regulation| boolean | Default FALSE. Adds columns containing regulation information by calling find_regulation(). |

## clean_massSpec_data
Standardizes names, drops unneeded columns, and filters out NaNs.
| Parameter | Type       | Description|
|-----------|------------|------------|
| dt  | data.table | dt to clean |
| remove_repeat_genes | boolean | Default FALSE. Removes repeat gnes, keeping the first instance. |
| keep_more_significant | boolean | Whether or not to, when facing repeated genes, keep only the most significant result by Qval. If false, the first of the repeated genes is arbitrarily removed instead. KEEP THIS TO FALSE UNLESS YOU HAVE A GREAT REASON!
|add_regulation| boolean | Default FALSE. Adds columns containing regulation information by calling find_regulation(). |
| keep_cols | vector | names of columns to keep alongside those kept by standard.
