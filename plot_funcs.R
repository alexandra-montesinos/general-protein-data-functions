library(ggplot2) #v3.5.1
library(ggnewscale) #v0.5.0
library(ggrepel) #v0.9.6

### Functions for creating volcano plots

good_theme <- theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16), 
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 18),
    axis.title.y = element_text(face = "bold", margin = margin(0,10,0,0), size = rel(0.9), color = 'black'),
    axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(10,0,0,0), size = rel(0.9), color = 'black'),
    axis.text = element_text(size = 18),
    panel.grid.major = element_line(color = "grey96", size = 0.5),
    #panel.grid.minor = element_line(color = "grey94", size = 0.25),
    legend.title = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 16),
  )

good_theme_big <- theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 26), 
    plot.subtitle = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 24),
    axis.title.y = element_text(face = "bold", margin = margin(0,10,0,0), size = rel(1.1), color = 'black'),
    axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(10,0,0,0), size = rel(1.1), color = 'black'),
    # panel.grid.major = element_line(size = 0.7, linetype = 'solid', color = "gray94"), 
    panel.grid.major = element_blank(),
    #panel.grid.minor = element_line(size = 0.25, linetype = 'solid', color = "gray92"),
    legend.title = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 20),
    axis.line = element_line(size = 0.8, linetype = 'solid', color = "black")
  )

### FUNCTION FOR ADDING GENE LABELS TO A DATATABLE BASED ON REGULATION AND QVAL
# @param dt = dt with genes and regulation info
# @param num_labels = number of genes to label in EACH GROUP of regulation and Qval. There will therefore be 2x this amount of labels total.
#                     There MAY be less total labels than this in the returned datatable if they are removed thanks to repeat labels,
#                     they have a broken label, or if you entered an uneven number of labels.
# @param even_regulation_split = whether or not to enforce an even split of labels acros up- and down-regulated genes.
#
# @return = datatable with new columns Top_Labels_LogFold & Top_Labels_Qvalue,Top_Labels_LogFold_Up, Top_Labels_LogFold_Down,
#           Top_Labels_Qvalue_Up, and Top_Labels_Qvalue_Down containing the gene names only for the top genes in those categories.
label_top <- function(dt, num_labels = 10, even_regulation_split = TRUE) {
  if(even_regulation_split) {
    ### SPLIT SELECTED GENES EVENLY ACROSS REGULATION
    #Select top upregulated genes
    up_genes <- dt[dt$AVG_Log2_Ratio > 0 & dt$Qvalue < 0.05, ]
    up_top <- up_genes %>%
      arrange(desc(abs(AVG_Log2_Ratio))) %>%
      slice_head(n = num_labels / 2) %>%
      pull(Genes)
    #Select top downregulated genes
    down_genes <- dt[dt$AVG_Log2_Ratio < 0 & dt$Qvalue < 0.05, ]
    down_top <- down_genes %>%
      arrange(desc(abs(AVG_Log2_Ratio))) %>%
      slice_head(n = num_labels / 2) %>%
      pull(Genes)
    #Combine lists, label the corresponding column
    top <- c(up_top, down_top)
    dt$Top_Labels_LogFold <- ifelse(dt$Genes %in% top, 
                                           dt$Genes, NA)
    #Do the same again for q-value genes
    up_top <- dt %>%
      filter(AVG_Log2_Ratio > 0) %>%
      arrange(Qvalue) %>%
      slice_head(n = num_labels / 2) %>%
      pull(Genes)
    down_top <- dt %>%
      filter(AVG_Log2_Ratio < 0) %>%
      arrange(Qvalue) %>%
      slice_head(n = num_labels / 2) %>%
      pull(Genes)
    top <- c(up_top, down_top)
    
    dt$Top_Labels_Qvalue <- ifelse(dt$Genes %in% top,
                                          dt$Genes, NA)
  } else {
    dt$Top_Labels_LogFold <- ifelse(dt$Genes %in% dt$Genes[order(abs(as.numeric(dt$AVG_Log2_Ratio)), decreasing = TRUE)][1:num_labels] & dt$Qvalue < 0.05, 
                                           dt$Genes, NA)
    dt$Top_Labels_Qvalue <- ifelse(dt$Genes %in% dt$Genes[order(dt$Qvalue)][1:num_labels], 
                                          dt$Genes, NA)
  }
  #Ensure there are no genes labeled more than once across these columns
  dt <- dt[, Top_Labels_Qvalue := fifelse(Top_Labels_Qvalue %in% Top_Labels_LogFold, NA_character_, Top_Labels_Qvalue)]
  
  #Remove labels from columns with a broken label. This usually happens because there isn't a gene name given for a datapoint
  dt$Top_Labels_Qvalue <- ifelse(dt$Top_Labels_Qvalue == "NaN" | dt$Top_Labels_Qvalue == "NULL" | dt$Top_Labels_Qvalue == "NA", 
                                        NA, dt$Genes)
  dt$Top_Labels_LogFold <- ifelse(dt$Top_Labels_LogFold == "NaN" | dt$Top_Labels_LogFold == "NULL", 
                                        NA, dt$Genes)
  #Every once in a while you have a gene name that's way too long to reasonably display. Let's truncate those!
  dt[, Top_Labels_Qvalue  := str_trunc(Top_Labels_Qvalue, width = 20, side = "right", ellipsis = "...")]
  dt[, Top_Labels_LogFold := str_trunc(Top_Labels_LogFold, width = 20, side = "right", ellipsis = "...")]
  
  dt$Top_Labels_LogFold_Up   <- ifelse(dt$Regulation == "Upregulated", 
                                            dt$Top_Labels_LogFold, NA)
  dt$Top_Labels_LogFold_Down <- ifelse(dt$Regulation == "Downregulated", 
                                            dt$Top_Labels_LogFold, NA)
  dt$Top_Labels_Qvalue_Up   <- ifelse(dt$Regulation == "Upregulated", 
                                            dt$Top_Labels_Qvalue, NA)
  dt$Top_Labels_Qvalue_Down <- ifelse(dt$Regulation == "Downregulated", 
                                            dt$Top_Labels_Qvalue, NA)
  
  #Ensure there are no genes labeled more than once across these columns
  dt <- dt[, Top_Labels_Qvalue := fifelse(Top_Labels_Qvalue %in% Top_Labels_LogFold, NA_character_, Top_Labels_Qvalue)]
  
  return(dt)
}



plot_volcano <- function(data, title, theme = good_theme_big, ratio_cutoff = 0, 
                         add_labels_log = TRUE, add_labels_q = TRUE, num_labels = 20) {
  if(add_labels_log || add_labels_q) {
    data <- label_top(data, num_labels)
  }
  
  num_up   <- sum(data$AVG_Log2_Ratio > 0 & data$Qvalue < 0.05)
  num_down <- sum(data$AVG_Log2_Ratio < 0 & data$Qvalue < 0.05)
  
  plot <- ggplot(data = data, aes(x = AVG_Log2_Ratio, y = -log10(Qvalue), col = Regulation)) + 
             geom_point(size = 2.2, alpha = 0.5) +
             geom_hline(yintercept = -log10(0.05), col="#949494", linetype = "dashed", size = 0.8) +
             #coord_cartesian(ylim = c(0, 6), xlim = c(-1, 1)) + #Limits y & x axes, while keeping the nonvisible data existing
             scale_color_manual(values = c("#068ECC", "#949494", "red"), 
                                labels = c("Downregulated", "Insignificant", "Upregulated")) + #Sets labels, overwriting categories from datatable
             guides(col = guide_legend(override.aes = list(size = 8, alpha = 1))) + #Overrides geom_point point size to make them larger in legend for visibility
             scale_x_continuous(breaks = seq(floor(min(data$AVG_Log2_Ratio)), ceiling(max(data$AVG_Log2_Ratio)), 1)) + #Forces steps of 1, regardless of range
             scale_y_continuous(limits = c(0, max(ceiling(max(-log10(data$Qvalue))), 2.5)),
                                breaks = seq(0, max(ceiling(max(-log10(data$Qvalue))), 2.5), 1)) +
             theme + 
             labs(title = title, subtitle = paste0("# Upregulated: ", num_up, ", # Downregulated: ", num_down))
  
  if(add_labels_q) {
    plot <- plot + 
      geom_label_repel(data = data,
                       aes(x = AVG_Log2_Ratio, y = -log10(Qvalue), label = Top_Labels_Qvalue_Down),
                       max.overlaps = Inf, size = 9, label.padding = 0.1, 
                       box.padding = 0.2, label.r = 0.5, fontface = "bold", 
                       force = 50, alpha = 0.8, segment.size = 0.6, 
                       segment.alpha = 0.8, min.segment.length = 0, 
                       fill = "white", color = "#068ECC", label.size = NA) +
      geom_label_repel(data = data,
                       aes(x = AVG_Log2_Ratio, y = -log10(Qvalue), label = Top_Labels_Qvalue_Up),
                       max.overlaps = Inf, size = 9, label.padding = 0.1, 
                       box.padding = 0.2, label.r = 0.5, fontface = "bold", 
                       force = 50, alpha = 0.8, segment.size = 0.6, 
                       segment.alpha = 0.8, min.segment.length = 0, 
                       fill = "white", color = "red", label.size = NA)
  }
  if(add_labels_log) {
    plot <- plot + 
      geom_label_repel(data = data,
                       aes(x = AVG_Log2_Ratio, y = -log10(Qvalue), label = Top_Labels_LogFold_Down),
                       max.overlaps = Inf, size = 9, label.padding = 0.1, 
                       box.padding = 0.2, label.r = 0, fontface = "bold", 
                       force = 50, alpha = 0.8, segment.size = 0.6, 
                       segment.alpha = 0.8, min.segment.length = 0, 
                       fill = "white", color = "#068ECC", label.size = NA) +
      geom_label_repel(data = data,
                       aes(x = AVG_Log2_Ratio, y = -log10(Qvalue), label = Top_Labels_LogFold_Up),
                       max.overlaps = Inf, size = 9, label.padding = 0.1, 
                       box.padding = 0.2, label.r = 0, fontface = "bold", 
                       force = 50, alpha = 0.8, segment.size = 0.6, 
                       segment.alpha = 0.8, min.segment.length = 0, 
                       fill = "white", color = "red", label.size = NA)
  }
  
  if(ratio_cutoff != 0) {
    plot <- plot + geom_vline(xintercept = c(ratio_cutoff, -ratio_cutoff), col="#949494", linetype = "dashed", size = 0.8)
  } else {
    plot <- plot + geom_vline(xintercept = 0, col="#949494", linetype = "dashed", size = 0.8)
  }
  
  return(plot)
}

### CREATE VOLCANO PLOT WITH 2 DATASETS, SHOWING AGREEMENT WITH COLOR
# This function is older, and therefore not quite as advanced as the basic plot_volcano() function
plot_volcano_w_agreement <- function(data1, data2, title, theme = good_theme_big, num_labels = 20) {
  data <- find_regulation(data1, data2)
  data <- label_top(data, num_labels)
  View(data)
  
  data$Color <- ifelse(data$Regulation == "Downregulated_Agrees", "#068ECC", 
                       ifelse(data$Regulation == "Upregulated_Agrees", "red", 
                              ifelse(data$Regulation == "Upregulated_Disagrees", "purple", 
                                     ifelse(data$Regulation == "Downregulated_Disagrees", "#FFAC1C", "#949494"))))
  
  return(ggplot(data, aes(x = AVG_Log2_Ratio, y = -log10(Qvalue), col = Regulation)) + 
           geom_point(size = 2.2, alpha = 0.5) +
           geom_vline(xintercept = 0, col="black", linetype="dashed") + #Adds threshold lines
           geom_hline(yintercept = -log10(0.05), col="black", linetype="dashed") +
           #coord_cartesian(ylim = c(0, 6), xlim = c(-1, 1)) + #Limits y & x axes, while keeping the nonvisible data existing
           scale_color_manual(values = c("Downregulated_Agrees" = "#068ECC", 
                                         "Downregulated_Disagrees" = "purple",
                                         "Not Significant" = "#949494",
                                         "Upregulated_Agrees" = "red", 
                                         "Upregulated_Disagrees" = "#ffb000"),
                              name = "Regulation", # Legend title
                              labels = c("Downregulated (Agrees)", "Downregulated (Disagrees)", "Not Significant",
                                         "Upregulated (Agrees)", "Upregulated (Disagrees)")) + 
           #scale_x_continuous(breaks = seq(-15, 15, 1)) + # to customize the breaks in the x axis
           #scale_y_continuous(breaks = seq(-0, 30, 1)) + # to customize the breaks in the y axis
           # geom_label_repel(aes(label = Top_Labels_Qvalue),
           #                  max.overlaps = Inf, size = 9, 
           #                  label.padding = 0.1, box.padding = 0.2, label.r = 0.5,
           #                  fontface = "bold", force = 16, alpha = 0.8, 
           #                  segment.size = 0.6, segment.alpha = 0.8, min.segment.length = 0, 
           #                  fill = "white", label.size = NA) +
           geom_label_repel(aes(label = Top_Labels_LogFold),
                            max.overlaps = Inf, size = 9, 
                            label.padding = 0.1, box.padding = 0.1, label.r = 0,
                            fontface = "bold", alpha = 0.8,
                            segment.size = 0.6, segment.alpha = 0.8, min.segment.length = 0, 
                            fill = "white", label.size = NA) +
           theme +
           # theme(legend.position = "right") +
           ggtitle(title)
  )
}

### CREATES VOLCANO PLOT USING 2 DATASETS, COLORING THEM DIFFERENTLY
plot_volcano_2data <- function(data1, data1_description, data2, data2_description, title, theme = good_theme_big, 
                               add_labels = TRUE, num_labels = 20) {
  if(add_labels) {
    data1 <- label_top(data1, num_labels / 2)
    data2 <- label_top(data2, num_labels / 2)
  }
  
  plot <- ggplot() +
    geom_point(data = data1, aes(x = AVG_Log2_Ratio, y = -log10(Qvalue), col = Regulation),
               size = 2.2, alpha = 0.5, shape = 16) +
    # geom_vline(xintercept = 0, col="gray", linetype="dashed", size = 0.8, color = "red") + 
    geom_hline(yintercept = -log10(0.05), col="gray", linetype="dashed", size = 0.8, color = "red") +
    scale_x_continuous(breaks = seq(floor(min(data1$AVG_Log2_Ratio)), ceiling(max(data1$AVG_Log2_Ratio)), 1)) + #Forces steps of 1, regardless of range
    scale_y_continuous(limits = c(0, max(ceiling(max(-log10(data1$Qvalue))), 2.5)),
                       breaks = seq(0, max(ceiling(max(-log10(data1$Qvalue))), 2.5), 1)) +
    theme +
    ggtitle(title) +
    scale_color_manual(values = c("#785ef0", "#949494", "#dc267f"),
                       labels = c(paste("Downregulated", data1_description), 
                                  paste("Insignificant", data1_description),
                                  paste("Upregulated", data1_description))) +
    new_scale_color() +
    geom_point(data = data2, aes(x = AVG_Log2_Ratio, y = -log10(Qvalue), color = Regulation),
               size = 2.2, alpha = 0.5, shape = 17) +
    scale_color_manual(values = c("#ffb000", "#949494", "#fe6100"),
                       labels = c(paste("Downregulated", data2_description), 
                                  paste("Insignificant", data2_description),
                                  paste("Upregulated", data2_description))) +
    # geom_label_repel(data = data1,
    #                  aes(x = AVG_Log2_Ratio, y = -log10(Qvalue), label = Top_Labels_Qvalue_Down),
    #                  max.overlaps = Inf, size = 9, label.padding = 0.1, 
    #                  box.padding = 0.2, label.r = 0.5, fontface = "bold", 
    #                  force = 16, alpha = 0.8, segment.size = 0.6, 
    #                  segment.alpha = 0.8, min.segment.length = 0, 
    #                  fill = "white", color = "#785ef0", label.size = NA) +
    # geom_label_repel(data = data1,
    #                  aes(x = AVG_Log2_Ratio, y = -log10(Qvalue), label = Top_Labels_Qvalue_Up),
    #                  max.overlaps = Inf, size = 9, label.padding = 0.1, 
    #                  box.padding = 0.2, label.r = 0.5, fontface = "bold", 
    #                  force = 16, alpha = 0.8, segment.size = 0.6, 
    #                  segment.alpha = 0.8, min.segment.length = 0, 
    #                  fill = "white", color = "#dc267f", label.size = NA) +
    geom_label_repel(data = data1,
                     aes(x = AVG_Log2_Ratio, y = -log10(Qvalue), label = Top_Labels_LogFold_Down),
                     max.overlaps = Inf, size = 9, label.padding = 0.1, 
                     box.padding = 0.2, label.r = 0, fontface = "bold", 
                     force = 16, alpha = 0.8, segment.size = 0.6, 
                     segment.alpha = 0.8, min.segment.length = 0, 
                     fill = "white", color = "#785ef0", label.size = NA) +
    geom_label_repel(data = data1,
                     aes(x = AVG_Log2_Ratio, y = -log10(Qvalue), label = Top_Labels_LogFold_Up),
                     max.overlaps = Inf, size = 9, label.padding = 0.1, 
                     box.padding = 0.2, label.r = 0, fontface = "bold", 
                     force = 16, alpha = 0.8, segment.size = 0.6, 
                     segment.alpha = 0.8, min.segment.length = 0, 
                     fill = "white", color = "#dc267f", label.size = NA) +
    # geom_label_repel(data = data2,
                     # aes(x = AVG_Log2_Ratio, y = -log10(Qvalue), label = Top_Labels_Qvalue_Down),
                     # max.overlaps = Inf, size = 9, label.padding = 0.1, 
                     # box.padding = 0.2, label.r = 0.5, fontface = "bold", 
                     # force = 16, alpha = 0.8, segment.size = 0.6, 
                     # segment.alpha = 0.8, min.segment.length = 0, 
                     # fill = "#f6f6f6", color = "#ffb000", label.size = NA) +
    # geom_label_repel(data = data2,
    #                  aes(x = AVG_Log2_Ratio, y = -log10(Qvalue), label = Top_Labels_Qvalue_Up),
    #                  max.overlaps = Inf, size = 9, label.padding = 0.1, 
    #                  box.padding = 0.2, label.r = 0.5, fontface = "bold", 
    #                  force = 16, alpha = 0.8, segment.size = 0.6, 
    #                  segment.alpha = 0.8, min.segment.length = 0, 
    #                  fill = "#f6f6f6", color = "#fe6100", label.size = NA) +
    geom_label_repel(data = data2,
                     aes(x = AVG_Log2_Ratio, y = -log10(Qvalue), label = Top_Labels_LogFold_Down),
                     max.overlaps = Inf, size = 9, label.padding = 0.1, 
                     box.padding = 0.2, label.r = 0, fontface = "bold", 
                     force = 16, alpha = 0.8, segment.size = 0.6, 
                     segment.alpha = 0.8, min.segment.length = 0, 
                     fill = "#f6f6f6", color = "#ffb000", label.size = NA) +
    geom_label_repel(data = data2,
                     aes(x = AVG_Log2_Ratio, y = -log10(Qvalue), label = Top_Labels_LogFold_Up),
                     max.overlaps = Inf, size = 9, label.padding = 0.1, 
                     box.padding = 0.2, label.r = 0, fontface = "bold", 
                     force = 16, alpha = 0.8, segment.size = 0.6, 
                     segment.alpha = 0.8, min.segment.length = 0, 
                     fill = "#f6f6f6", color = "#fe6100", label.size = NA)
  
            
  return(plot)
}

