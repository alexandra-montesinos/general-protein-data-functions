library(dplyr) #v1.1.4
library(stringr) #v1.5.1
library(tidyr) #v1.3.1

### Contains functions for cleaning datatables in order to create consistent naming styles. It is often required to use these functions before calling my other datatable functions, 
# especially those in plot_funcs.R

find_regulation <- function(dt, dt2 = NULL, ratio_cutoff = 0) {
  #Add a column to the data frame to specify if they are UP- or DOWN- regulated
  dt$Regulation <- "Not Significant"
  
  #If we have a second datatable, add agreement info
  if(is.null(dt2)) {
    #if log2Foldchange > 0 and pvalue < 0.05
    dt$Regulation[dt$AVG_Log2_Ratio >= 0 & dt$Qvalue < 0.05] <- "Upregulated"
    #if log2Foldchange < 0 and pvalue < 0.05
    dt$Regulation[dt$AVG_Log2_Ratio < 0 & dt$Qvalue < 0.05] <- "Downregulated"
  } else {
    dt <- semi_join(dt, dt2, by = "Genes") #Remove any genes not present in both datasets
    dt2 <- semi_join(dt2, dt, by = "Genes") 
    dt  <- dt[order(dt$Genes),] #Order both by genes
    dt2 <- dt2[order(dt2$Genes),]
    
    #if log2Foldchange > 0 and pvalue < 0.05
    dt$Regulation[dt$AVG_Log2_Ratio >= 0 & dt$Qvalue < 0.05 & dt2$AVG_Log2_Ratio >= 0] <- "Upregulated_Agrees"
    #if log2Foldchange > 0 and pvalue < 0.05
    dt$Regulation[dt$AVG_Log2_Ratio >= 0 & dt$Qvalue < 0.05 & dt2$AVG_Log2_Ratio < 0] <- "Upregulated_Disagrees"
    #if log2Foldchange > 0 and pvalue < 0.05
    dt$Regulation[dt$AVG_Log2_Ratio < 0 & dt$Qvalue < 0.05 & dt2$AVG_Log2_Ratio < 0] <- "Downregulated_Agrees"
    #if log2Foldchange > 0 and pvalue < 0.05
    dt$Regulation[dt$AVG_Log2_Ratio < 0 & dt$Qvalue < 0.05 & dt2$AVG_Log2_Ratio > 0] <- "Downregulated_Disagrees"
  }

  if(ratio_cutoff != 0) {
    dt$Regulation[dt$AVG_Log2_Ratio < ratio_cutoff & dt$AVG_Log2_Ratio > -ratio_cutoff] <- "Not Significant"
  }
  
  return(dt)
}

### FUNCTION FOR GENERAL DATATABLE CLEANING
# Used to standardize different styles of data naming styles. It is often necessary to use this function if you are calling my other datatable functions,
# so you'll have to use it even if you hate my particular naming style. Sorry! 
# If your dt contains a gene column, use clean_gene_dt instead. It will call this function for you.
# @param dt = datatable to clean
clean_dt <- function(dt) {
  #You were supposed to pass in a datatable, not dataframe >:(
  if(class(dt)[1] == "data.frame") {
    dt <- as.data.table(dt)
  }
  
  ### COLUMN NAME EDITS
  # names(dt) <- str_to_title(names(dt)) #Capitalizes each separate word
  names(dt) <- gsub(" ", "_", colnames(dt)) #Replaces spaces with underscores
  setnames(dt, gsub("[()'\"`,]", "", names(dt))) #Remove problematic characters like (), '', "", ``, and ,
  setnames(dt, old = names(dt), new = gsub("^([a-z])", "\\U\\1", names(dt), perl = TRUE)) #Capitalizes first letter
  
  #Deal with poorly named columns
  dt <- dt %>% rename_with(~ str_c("PREFIX_", .), matches("^\\d")) #Cols that start with a number
  names(dt)[names(dt) == 'NA'] <- 'PREFIX_NA' #Cols that are just called NA. WHY DO PEOPLE DO THIS??
  names(dt)[names(dt) == 'NaN'] <- 'PREFIX_NaN'
  names(dt)[names(dt) == 'NULL'] <- 'PREFIX_NULL'
  
  return(dt)
}

#Cleans a datatable with gene data, including throuhg a call of clean_dt().
# WARNING: This function removes rows with NAs and NaNs in the Genes column. These missing values (especially NaNs) are sometimes caused because
# certain proteins in our sample don't have gene names. In these cases, dropping them like this will just remove valid data from the set.
# Of course, this should be a pretty rare occurence. Only ever happened to me once!
#
# @param dt = dt to clean
# @param remove_repeat_genes = Does what it says on the tin
# @param add_regulation      = Adds regulation information through find_regulation(). The dt must contain a AVG_Log2_Ratio & Qvalue column.
clean_gene_dt <- function(dt, remove_repeat_genes = FALSE, add_regulation = FALSE) {
  dt <- clean_dt(dt) #Does a general clean
  
  if("Gene" %in% names(dt)) {
    dt <- dt %>% dplyr::rename(Genes = Gene)
  } else if("Gene_Names" %in% names(dt)) {
    dt <- dt %>% dplyr::rename(Genes = Gene_Names)
  }
  dt <- dt %>% dplyr::filter(Genes != "NaN") 
  dt <- dt %>% drop_na(Genes) #Drop NAs
  print(is.data.table(dt))
  #standardizes gene names for easy comparison across datasets
  dt[, Genes := paste0(toupper(substr(Genes, 1, 1)), tolower(substr(Genes, 2, nchar(Genes))))]
  
  if(add_regulation) {
    dt <- find_regulation(dt)
  }
  
  if(remove_repeat_genes) {
    dt <- dt %>% distinct(Genes, .keep_all = TRUE)
  }
  
  return(dt)
}

remove_less_significant_repeats <- function(dt, col, sig_col = "Qvalue") {
  setorderv(dt, cols = c(col, sig_col)) #Order by significance
  dt_filtered <- dt[!duplicated(dt[[col]])] #Keep first occurence
  return(dt_filtered)
}

### FUNCTION FOR CLEANING STANDARD MASS-SPEC DATATABLE
# @param dt == datatable to clean 
# @param bool remove_repeats        == Whether or not to remove repeated genes.
# @param bool keep_more_significant == Whether or not to, when facing repeated genes, keep only the most significant result by Qval.
#                                       If false, the first of the repeated genes is arbitrarily removed instead. KEEP THIS TO FALSE UNLESS YOU HAVE A GREAT REASON!
# @param bool add_regulation        == Whether or not to add a Regulation column containing regulation information ("Upregulated", "Downregulated", or "Not Significant").
# @param keep_cols                  == Vector of additional column names (strings) to keep.
# @return             == Cleaned datatable
clean_massSpec_data <- function(dt, remove_repeats = TRUE, keep_more_significant = FALSE, add_regulation = TRUE, keep_cols = NULL) {
  new_dt <- data.table("Genes" = dt$Genes, 
                       "UniProtIDs" = dt$UniProtIds,
                       "AVG_Log2_Ratio" = as.numeric(dt$`AVG Log2 Ratio`),
                       "Pvalue" = as.numeric(dt$Pvalue),
                       "Qvalue" = as.numeric(dt$Qvalue),
                       "Num_Unique_Peptides" = as.numeric(dt$`# Unique Total Peptides`))
  if(!is.null(keep_cols)) {
    for(col in keep_cols) {
      new_dt[[col]] <- dt[[col]]
    }
  }
  
  new_dt <- new_dt[ !(new_dt$AVG_Log2_Ratio == 'NaN'),] #Filter out NaNs
  new_dt <- new_dt[ !(new_dt$Pvalue == 'NaN'),]
  new_dt <- new_dt[ !(new_dt$Qvalue == 'NaN'),]
  new_dt <- new_dt %>% filter(Genes != "NaN") #Some proteins don't have gene names. We'll drop those.
  
  if(keep_more_significant & remove_repeats) {
    new_dt <- remove_less_significant_repeats(new_dt, "Genes", "Qvalue")
  } else if(remove_repeats) {
    new_dt <- new_dt %>% distinct(Genes, .keep_all = TRUE)
  }
  
  if(add_regulation) {
    new_dt <- find_regulation(new_dt)
  }
    
  return(new_dt)
}

