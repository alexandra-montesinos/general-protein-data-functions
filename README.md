# general-protein-data-functions
Functions that clean &amp; plot protein data.

# clean_data.R
## find_regulation 
| Parameter | Type       | Description|
|-----------|------------|------------|
| dt  | data.table | finds regulation of proteins in this dt |
| dt2 | data.table | Optional. If included, will find if the regulation of the first datatable agrees with that of the second |
|ratio_cutoff| numeric | Default 0. If >0, will add a cutoff to what is considered significant up or down-regulation, marking those that don't make it as Not Significant|
