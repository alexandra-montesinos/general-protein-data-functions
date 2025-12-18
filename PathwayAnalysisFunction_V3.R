# install.packages("clusterProfiler")
# install.packages("ggplot2")
# install.packages("dplyr")
# BiocManager::install("org.Hs.eg.db", character.only = TRUE)
# BiocManager::install("org.Mm.eg.db", character.only = TRUE)

library(clusterProfiler)
library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

source("../../General Use Code/plot_funcs.R")

### Functions for running and plotting pathway analyses.

#pathway_analysis
# Runs a pathway analysis on given datatable containing gene expression data.
# Filters for significantly upregulated genes and performs Gene Ontology (GO) enrichment analysis.
# Returns a plotted linear analysis.
# @param gene_list == list of genes to perform analysis on, should be in SYMBOL form
# @param universe_gene_list == list of genes to compare against--the baseline of all genes that could've been selected, also in SYMBOL form
# @param OrgDb == OrgDb list (e.g. "org.Hs.eg.db", "org.Mm.eg.db")
# @param keyType == key type for gene list (e.g. "SYMBOL" (gene names), "UNIPROT" (uniprot ids))
# @param should_filter == bool, whether or not to use gofilter()--this should usually be TRUE
# @param filter_level == numeric, filter level to pass into gofilter(). The higher the level, the more specific the pathways, located deeper in the hierarchy.
#                        2-3 = Very broad, high-level terms, 5+ = more specific terms that represent detailed biological processes or functions.
# @param title == title string to give output graph
# @param theme == theme to use for plotting, default = good_theme_big
# @param num_categories == Number of categories to show on the plot
# @param filter_type == metric to use for filtering out insignificant pathways ("pvalue", "p.adjust", or "qvalue")
#
# @return == Named list: Plot = dotplot, Result = filtered enrichGo results.
pathway_analysis <- function(gene_list, universe_gene_list, OrgDb, keyType = "SYMBOL", should_filter = TRUE, 
                             filter_level = 5, title = "UNNAMED", theme = good_theme_big, num_categories = 15, filter_type = "p.adjust") {
  #Pass to enrichGO to automatically run FDR control and enrich column list
  # EnrichGo claims to hold all possible categories each gene might fit into.
  # Universe = all genes that could've been selected. If we had a randomly selected group,
  # what would we expect? Compares this with passed-in group of significant genes,
  # and determines which are upregulated, as scene in the produced chart. 
  enrichGO_results <- enrichGO(as.character(gene_list), OrgDb = OrgDb,
                               universe = as.character(universe_gene_list),
                               keyType = keyType,
                               ont = "BP")
  
  if(is.null(enrichGO_results)) {
    message("No enrichment results found, returning NULL")
    return(NULL)
  }
  if(should_filter) {
    enrichGO_filtered <- gofilter(enrichGO_results, filter_level)
  } else {enrichGO_filtered <- enrichGO_results}
  
  result <- as.data.table(enrichGO_filtered@result)
  result <- result %>%
    mutate(GeneRatio = sapply(GeneRatio, function(x) {
      parts <- strsplit(x, "/")[[1]]
      as.numeric(parts[1]) / as.numeric(parts[2])
    }))
  
  result_filtered <- result %>% dplyr::filter(result[[filter_type]] < 0.05)
  
  plot <- plot_pathway_analysis(result_filtered, title, theme, num_categories, filter_type)
  
  if(length(plot$data$ID) == 0) {
    message("Plotting failed! Gene list may have contained NAs")
  }
  
  return_list <- list(
    "Plot" = plot,
    "Result" = result
  )

  return(return_list)
}

#pathway_analysis_comp
#Run a pathway analysis on two gene lists, then create a plot including ONLY significantly represented by the first gene list.
# All params same as pathway_analysis, except:
# @param gene_list2 == second gene list for comparison
pathway_analysis_comp <- function(gene_list, gene_list2, universe_gene_list, OrgDb, keyType = "SYMBOL", should_filter = TRUE, 
                                  filter_level = 5, title = "UNNAMED", theme = good_theme_big, num_categories = 10, num_categories_comp = Inf) {
  result1 <- pathway_analysis(gene_list, universe_gene_list, OrgDb, keyType, should_filter, 
                              filter_level, title, theme, num_categories)
  result2 <- pathway_analysis(gene_list2, universe_gene_list, OrgDb, keyType, should_filter, 
                              filter_level, title, theme, num_categories)
  
  result1_filtered <- result1$Result %>% dplyr::filter(p.adjust < 0.05)
  result2_filtered <- result2$Result %>% dplyr::filter(p.adjust < 0.05)
  result_unique <- result1_filtered %>% dplyr::filter(!Description %in% result2_filtered$Description)
  
  plot_unique_pathways <- plot_pathway_analysis(result_unique, title, theme, num_categories_comp)
  
  return_list <- list(
    "Results1" = result1,
    "Results2" = result2,
    "ResultsComp" = list(
      "Result" = result_unique,
      "Plot" = plot_unique_pathways
    )
  )
  
  return(return_list)
}

plot_pathway_analysis <- function(dt, title, theme = good_theme, num_categories = 10, filter_type = "p.adjust") {
  filter_title <- stringr::str_to_title(gsub("\\.", " ", filter_type)) #Real quickly makes a nicer filter type title
  
  dt_subset <- dt[order(dt$p.adjust), ][1:min(num_categories, nrow(dt)), ]
  #Some of these pathway names are WAY too long, let's cut them off for displaying.
  # BUT ONLY FOR DISPLAYING--otherwise, pathways with the same name AFTER truncating 
  # will be considered the same for plotting logic, which will cause display errors.
  # description_labels <- setNames(str_trunc(dt_subset$Description, width = 50, side = "center"), dt_subset$Description)
  description_labels <- setNames(str_wrap(dt_subset$Description, width = 50), dt_subset$Description)
  
  return(ggplot(dt_subset, aes(x = GeneRatio, y = reorder(Description, GeneRatio), size = Count, color = dt_subset[[filter_type]])) +
          geom_point() +
          scale_size(range = c(4, 14)) +
          scale_y_discrete(labels = description_labels) +
          scale_color_continuous(low = "#df6664", high = "#377eb9", name = filter_title, 
                                 guide = guide_colorbar(reverse = TRUE)) +
          ggtitle(title) +
          good_theme_big + 
          theme(
            plot.title.position = "plot", #Positions the title relative to the *whole* plot, not just the panel
            plot.title = element_text(hjust = 1),
            # panel.grid.major = element_line(size = 0.7, linetype = 'solid', color = "gray90"), 
            panel.grid.major = element_blank(), 
            axis.text = element_text(size = 23),
          ) +
          labs(x = "Gene Ratio", y = NULL))
}

geneID_as_list <- function(dt, pathway) {
  genes <- dt$geneID[pathway]
  gene_list <- as.list(strsplit(genes, '\\/+'))
  return(gene_list)
}
