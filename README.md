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

# plot_funcs.R
## label_top
Adds columns containing labels (names) of genes based on their regulation and qvalue. Genes without labels will have an NA in this column.
| Parameter | Type       | Description|
|-----------|------------|------------|
| dt  | data.table | dt with genes and regulation info |
| num_labels | numeric | number of genes to label in EACH GROUP of regulation and Qval. There will therefore be 2x this amount of labels total. There MAY be less total labels than this in the returned datatable if they are removed thanks to repeat labels, they have a broken label, or if you entered an uneven number of labels. |
| even_regulation_split | boolean | Default TRUE. Whether or not to enforce an even split of labels acros up- and down-regulated genes.|

## plot_volcano
Create a volcano plot for log value based on a dataset.
| Parameter | Type       | Description|
|-----------|------------|------------|
| data  | data.table | data to plot. Please call clean_dt or one of its variants on this data first. |
| title | string | Title for the plot.|
| theme | ggplot theme object | Theme to use. Default good_theme_big, as defined in this code file. |
| ratio_cutoff | numeric | Cutoff for significant log ratio. Default 0, as is standard. |
|add_labels_log| Boolean | whether to add gene labels to the plot based on significant log value. Default TRUE.|
|add_labels_q | Boolean | whether to add gene labels to the plot based on significant Qvalue. Default TRUE.|
|num_labels | numeric | number of labels to add if you are adding labels at all. Default 20. |

## plot_volcano_w_agreement
Create a volcano plot with 2 datasets, showing regulation direction agreement with color. This function is older, and therefore not quite as advanced as the basic plot_volcano() function.
| Parameter | Type       | Description|
|-----------|------------|------------|
| data1  | data.table | data to plot. Please call clean_dt or one of its variants on this data first. |
| data2 | data.table | Second data to plot. Please call clean_dt or one of its variants on this data first. |
| title | string | Title for the plot.|
| theme | ggplot theme object | Theme to use. Default good_theme_big, as defined in this code file. |
|num_labels | numeric | number of labels to add if you are adding labels at all. Default 20. |

## plot_volcano_2data
Create a volcano plot with 2 datasets, coloring them differently. 
| Parameter | Type       | Description|
|-----------|------------|------------|
| data1  | data.table | data to plot. Please call clean_dt or one of its variants on this data first. |
| data1_description | string | Short description of data1 for use in the plot key.|
| data2 | data.table | Second data to plot. Please call clean_dt or one of its variants on this data first. |
| data2_description | string | Short description of data2 for use in the plot key.|
| title | string | Title for the plot.|
| theme | ggplot theme object | Theme to use. Default good_theme_big, as defined in this code file. |
| ratio_cutoff | numeric | Cutoff for significant log ratio. Default 0, as is standard. |
|add_labels| Boolean | whether to add gene labels to the plot. Default TRUE.|
|num_labels | numeric | number of labels to add if you are adding labels at all. Default 20. |
