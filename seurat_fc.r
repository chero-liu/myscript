source("/home/chenglong.liu/RaD/myscript/clolors.r")
createDir <- function(dir){
    if(!dir.exists(dir)){
        dir.create(dir, recursive=T)
    }
}

capitalize <- function(x) {
  sapply(strsplit(x, " "), function(x) {
    paste(toupper(substring(x, 1, 1)), tolower(substring(x, 2)), sep = "", collapse = " ")
  })
}

upper <- function(x) {
  sapply(strsplit(x, " "), function(x) {
    paste(toupper(substring(x, 1, 1)), toupper(substring(x, 2)), sep = "", collapse = " ")
  })
}

write_to_file <- function(text, file) {
  cat(text, file = file, append = TRUE)
}

text_similarity <- function(str1, str2) {
  # 处理空字符串情况
  if (is.null(str1) || is.null(str2) || str1 == "" || str2 == "") {
    return(0.0)
  }
  
  # 1. 计算编辑距离相似度
  # 使用adist计算编辑距离，然后转换为相似度
  edit_dist <- adist(str1, str2, ignore.case = TRUE)
  max_len <- max(nchar(str1), nchar(str2))
  edit_sim <- ifelse(max_len == 0, 0, 1 - edit_dist / max_len)
  
  # 2. 计算Jaccard相似度
  # 拆分为字符集合（针对英文，保留大小写差异）
  chars1 <- strsplit(str1, "")[[1]]
  chars2 <- strsplit(str2, "")[[1]]
  
  set1 <- unique(chars1)
  set2 <- unique(chars2)
  
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  jaccard_sim <- ifelse(union == 0, 0, intersection / union)
  
  # 加权平均 (编辑距离权重0.7, Jaccard权重0.3)
  return(0.7 * edit_sim + 0.3 * jaccard_sim)
}

get_filepath <- function(celltype, tissue = 'default') {
  # Read whitelist index
  whitelist <- read.csv("/data/database/sc_subtype_refmarker/whitelist_index.csv", stringsAsFactors = FALSE)
  
  # Step 1: Find best matching celltype using similarity scoring
  if (nrow(whitelist) == 0) {
    stop("Whitelist is empty")
  }
  
  # Calculate similarity scores
  scores <- sapply(whitelist$celltype, function(x) text_similarity(celltype, x))
  best_match_idx <- which.max(scores)
  best_celltype <- whitelist$celltype[best_match_idx]
  
  # Only proceed if similarity score > 0.1
  if (scores[best_match_idx] <= 0.1) {
    print(paste("No sufficiently matching celltype found for:", celltype))
    return(NULL)
  }
  
  # Step 2: Filter records for the best matching celltype
  celltype_records <- whitelist[whitelist$celltype == best_celltype, ]
  
  # Step 3: Find matching tissue or default
  tissue_match <- celltype_records[celltype_records$tissue == tissue, ]
  if (nrow(tissue_match) == 0) {
    # No matching tissue, use default
    default_match <- celltype_records[celltype_records$tissue == "default", ]
    if (nrow(default_match) == 0) {
      print(paste("No default tissue found for celltype:", best_celltype))
      return(NULL)
    }
    return(default_match$filepath[1])
  } else {
    return(tissue_match$filepath[1])
  }
}


standardizeDelimitedVector <- function(input_vec, delimiter = ",") {
  processed <- as.character(input_vec)
  processed <- paste(processed, collapse = delimiter)
  processed <- unlist(strsplit(processed, split = delimiter, fixed = TRUE))
  processed <- trimws(processed)
  processed <- processed[nzchar(processed)]
  return(processed)
}


getPredicateVars <- function(predicate) {


  if(is.null(predicate)){
    return(NULL)
  }

  if(predicate == "all"){
    return(NULL)
  }

  if (!is.character(predicate)) {
    stop("Input must be a character string")
  }

  matches <- regmatches(
    predicate,
    gregexpr("\\b([[:alnum:]._]+)\\s*%in%", predicate, perl = TRUE)
  )[[1]]
  
  variables <- unique(gsub("\\s*%in%", "", matches))
  
  if (length(variables) == 0) {
    warning("No variables found matching the pattern")
    return(character(0))
  }
  
  return(variables)
}

getDownsampleRDS <- function(rds, 
                      group_by = NULL,
                      target_cells=29000, 
                      error_threshold=1000,
                      verbose = TRUE) {
  set.seed(2025)
  target_cells = as.numeric(target_cells)
  error_threshold = as.numeric(error_threshold)
  if (!inherits(rds, "Seurat")) {
    stop("Input 'rds' must be a Seurat object.")
  }
  if (target_cells <= 0) {
    stop("'target_cells' must be a positive integer.")
  }
  if (error_threshold < 0) {
    stop("'error_threshold' must be a non-negative number.")
  }

  if (is.null(group_by)) {
    group_by <- "clusters"
  }

  if (!group_by %in% colnames(rds@meta.data)) {
    stop("'group_by' must be a valid column name in the metadata.")
  }

  if (dim(rds@meta.data)[1] <= (target_cells+error_threshold)){
    return(list(optimal_subset = rds))
  }

  if (group_by == "group"){
    rds$combo <- paste(rds@meta.data$clusters,rds@meta.data[,group_by],sep = "_")
  }else{
    rds$combo <- paste(rds@meta.data$group,rds@meta.data[,group_by],sep = "_")
  }

  Idents(rds) =  'combo'

  group_counts <- table(rds@meta.data[, 'combo'])
  upper_bound <- max(group_counts)
  
  low <- 3
  high <- upper_bound
  history <- list()
  
  max_iter = target_cells/1000
  for (i in 1:max_iter) {
    if (high - low <= 2) {
      if (verbose) message("Reached minimum search interval")
      break
    }
    
    current_downsample <- floor((low + high) / 2)
    
    set.seed(1236)
    subrds <- subset(rds, downsample = current_downsample)
    actual_cells <- dim(subrds@meta.data)[1]
    error <- abs(actual_cells - target_cells)
    
    history[[i]] <- list(iteration = i,
                         downsample_numb = current_downsample,
                         actual_cells = actual_cells,
                         error = error)
    
    if (error <= error_threshold) {
      if (verbose) message(paste("Converged at iteration", i))
      break
    }
    
    if (actual_cells < target_cells) {
      low <- current_downsample
    } else {
      high <- current_downsample
    }
  }
  
  print(paste0(actual_cells," cells in downsampled RDS"))
  history_df <- do.call(rbind, lapply(history, as.data.frame))
  
  best_index <- which.min(history_df$error)
  return(list(
    optimal_subset = subset(rds, downsample = history_df$downsample_numb[best_index]),
    iteration_history = history_df,
    final_downsample = history_df$downsample_numb[best_index],
    final_error = history_df$error[best_index]
  ))
}

subsetRDS <- function(rds, sampleid = "all", group = "all", cluster = "all", new_celltype = "all", predicate = 'all') {
  library(stringr)
  if (new_celltype != "all") {
    new_celltype_values <- standardizeDelimitedVector(new_celltype,',')
    rds <- rds[, rds@meta.data$new_celltype %in% new_celltype_values]
  }
  
  if (sampleid != "all") {
    sampleid_values <- standardizeDelimitedVector(sampleid,',')
    rds <- rds[, rds@meta.data$sampleid %in% sampleid_values]
  }
  
  if (group != "all") {
    group_values <- standardizeDelimitedVector(group,',')
    rds <- rds[, rds@meta.data$group %in% group_values]
  }
  
  if (cluster != "all") {
    cluster_values <- standardizeDelimitedVector(cluster, ",")
    rds <- rds[, rds@meta.data$clusters %in% cluster_values]
  }

  if (predicate != 'all') {
      df <- rds@meta.data
      desired_cells <- subset(df, eval(parse(text = predicate)))
      rds <- rds[, rownames(desired_cells)]
  }

  use_col = unique(c("new_celltype",'sampleid','group','clusters',getPredicateVars(predicate)))
  for (Colname in use_col){
    if (Colname %in% colnames(rds@meta.data)){
      rds = standardizeLevels(rds,Colname)
    }
  }

  return(rds)
}

sort_char_numbers <- function(x) {
  numeric_x <- as.numeric(x)
  x[order(numeric_x)]
}

is_char_number <- function(x) {
  grepl("^[-+]?\\d+(\\.\\d+)?([eE][-+]?\\d+)?$", x) & is.character(x)
}

standardizeColors <- function(rds, use_col, palette = 'customecol2'){
    if (is.null(palette)){
      palette = 'customecol2'
    }
    for (Colname in use_col){
    if (Colname %in% colnames(rds@meta.data)){
      rds = standardizeLevels(rds,Colname)
      if (!paste0(Colname,"_col") %in% colnames(rds@meta.data)){
        if (Colname %in%  c('group','sampleid')){
            rds@meta.data = addColorToDataFrame(rds@meta.data,Colname,palette = 'group_color')
        }else{
            rds@meta.data = addColorToDataFrame(rds@meta.data,Colname,palette = palette)
        }
      }else if (!is2GroupEquality(rds,Colname,paste0(Colname,"_col"))){
        rds@meta.data[, paste0(Colname,"_col")] = NULL
        rds@meta.data = addColorToDataFrame(rds@meta.data,Colname,palette = palette)
      }
    }
    }
    return(rds)
}

standardizeLevels <- function(rds, Colname, Colname_levels = 'no') {
  if (!Colname %in% colnames(rds@meta.data)) {
    stop("Column '", Colname, "' not found in metadata.")
  }
  safe_factor <- function(x, levels) {
    factor(as.character(x), levels = levels)
  }
  Colname_levels = standardizeDelimitedVector(Colname_levels,',')
  if (Colname_levels == 'no') {
    current <- rds@meta.data[[Colname]]
    if (is.null(levels(current))) {
      if (is_char_number(current[1])) {
        levs <- sort_char_numbers(unique(current))
      } else {
        levs <- levels(reorderLevelsByFrequency(current))
      }
      rds@meta.data[[Colname]] <- safe_factor(current, levs)
    } else {
      levs <- intersect(levels(current), unique(current))
      rds@meta.data[[Colname]] <- safe_factor(current, levs)
    }
  } else {

    if (!all(unique(rds@meta.data[[Colname]]) %in% Colname_levels)) {
      warning("Some values not in provided levels.")
    }

    if (length(unique(rds@meta.data[[Colname]])) > length(Colname_levels)) {
      stop(paste0("列指定levels数量必须大于等于实际列中的分组数量"))
    }

    levs <- intersect(Colname_levels, unique(rds@meta.data[[Colname]]))
    rds@meta.data[[Colname]] <- safe_factor(rds@meta.data[[Colname]], levs)
  }
  return(rds)
}

mapRDS <- function(rds,anno_list){
    
    rds$new_celltype = as.character(rds$new_celltype)

    for (anno_name in names(anno_list)){
        anno = anno_list[[anno_name]]
        rownames(anno) = anno$Barcode
        rm_barcode = rds@meta.data[rds$new_celltype == anno_name,'rawbc'][!rds@meta.data[rds$new_celltype == anno_name,'rawbc'] %in% anno$Barcode]
        if (length(rm_barcode) != 0){
            rds = rds[,!rds$rawbc %in% rm_barcode]
            print(paste("remove",length(rm_barcode),"barcode from",anno_name))
        }
        anno = anno[rds@meta.data[rds$new_celltype == anno_name,'rawbc'],]
        rds@meta.data[rds$new_celltype == anno_name,'new_celltype'] = anno$new_celltype
    }
    if('new_celltype_col' %in% colnames(rds@meta.data)){
        rds@meta.data = rds@meta.data[,-which(colnames(rds@meta.data)=='new_celltype_col')]
    }
    
    return(rds)
}

#------plot functions----------------

save_plot <- function(plot, filename, file_formats = c("png","pdf"), height = 6, width = 6, dpi = 300) {
  library(ggplot2)
  if (!dir.exists(filename)) {
    createDir(dirname(filename))
  }
  for (format in file_formats) {
    ggsave(
      filename = paste0(filename, ".", format),
      plot = plot,
      height = height,
      width = width,
      dpi = dpi,
      limitsize = FALSE,
      bg="white"
    )
  }
}


plotDensity <- function(rds,reduction = 'umap') {
    library(ggplot2)
    plot_data = Seurat::Embeddings(rds, reduction = reduction)
    xx=(max(plot_data[,1])-min(plot_data[,1]))/10
    yy=(max(plot_data[,2])-min(plot_data[,2]))/10
    xlim = c(floor(min(plot_data[,1])-xx), ceiling(max(plot_data[,1])+xx))
    ylim = c(floor(min(plot_data[,2])-yy), ceiling(max(plot_data[,2])+yy))
    ptsize = abs(xlim[2]-xlim[1])*abs(ylim[2]-ylim[1])*0.4/nrow(plot_data)
    plot_data = as.data.frame(plot_data)
    p = ggplot() + 
        stat_density_2d(geom="raster",contour = F,
                        data = plot_data, alpha=1, 
                        aes_string(x = colnames(plot_data)[1], y = colnames(plot_data)[2], fill="..density..")) + 
        scale_fill_viridis_c(option = "magma") + 
        geom_point(data =plot_data, 
                    aes_string(x = colnames(plot_data)[1], y = colnames(plot_data)[2]),
                    size = ptsize,shape = 16 ,color="white",alpha=0.7) + 
        xlim(xlim) + ylim(ylim) + 
        theme(legend.position="none", plot.title = element_text(hjust = 0.5),
                axis.ticks.x=element_blank(), axis.text.x=element_blank(),
                axis.ticks.y=element_blank(), axis.text.y=element_blank(),
                panel.background = element_rect(fill = 'black', colour = 'black'), 
                panel.grid.major = element_line(colour = "black"),
                panel.grid.minor = element_line(colour = "black"))
}

plotFeaturePlot = function(rds, 
                           gene, 
                           outdir = NULL,
                           reduct = 'umap', 
                           splitby = NULL, 
                           height = 6,
                           width = 6,
                           dpi = 300,
                           pointsize = 0.5) {
    library(RColorBrewer)
    ggfeatures = Seurat::FeaturePlot(rds,reduction= reduct,
                        features = gene, split.by = splitby,
                        keep.scale = "all",
                        max.cutoff = 'q99',
                        order = T,
                        pt.size = pointsize,
                        # ncol = ifelse( is.null(splitby) ,yes = ifelse(length(gene) >1, yes=2,no=1), no =2 ) 
                        ) &
                # ggplot2::scale_colour_gradientn(colours = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))(100) ) &
                ggplot2::scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "Spectral"))) &
                ggplot2::theme(legend.position = "none",
                    legend.margin = ggplot2::margin(0,0.5,0,0, unit = "cm"),
                    plot.margin = ggplot2::margin(0.1, 0.1, 0.1, 0.1, unit = "cm"),
                    plot.title = ggplot2::element_text(hjust = 0.5) )

    ggm = ggpubr::ggarrange(ggfeatures,common.legend = T, legend = "right")

    if (!is.null(outdir)) {
        save_plot(ggm, file.path(outdir, paste0(gene, "_FeaturePlot")), file_formats = c("pdf", "png"), height = height, width = width, dpi = dpi)
    }else{
        return(ggm)
    }
}


plotVlnPlot<- function(rds, 
                      gene, 
                      outdir = NULL, 
                      pt.size = 0.1, 
                      groupby = groupby,
                      height = 6,
                      width = 6,
                      dpi = 300,
                      plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm")) {
    p<- Seurat::VlnPlot(rds, 
            features = gene, 
            pt.size = pt.size, 
            group.by = groupby, 
            cols = c("#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b",
    "#666666","#1b9e77","#d95f02","#7570b3","#d01b2a","#43acde",
    "#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff",
    "#ecb9e5","#813139","#743fd2","#434b7e","#e6908e","#214a00",
    "#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7",
    "#006874","#d2ad2c","#b5d7a5","#9e8442","#4e1737","#e482a7",
    "#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99",
    "#3fdacb","#bf5b17"))  + 
    xlab("") + ylab("") + ggtitle(gene) + 
    theme(legend.position = "none", 
        axis.text.x = element_text(angle = 90),
        plot.margin = plot.margin ) 

    if (!is.null(outdir)) {
        save_plot(p, file.path(outdir, paste0(gene, "_VlnPlot")), file_formats = c("pdf", "png"), height = height, width = width, dpi = dpi)
    }else {
        return(p)
    }
    
}


plotVlnWithPValues <- function(data, 
                               genes, 
                               x, 
                               y, 
                               group, 
                               outdir, 
                               colors, 
                               ylim = c(0, NA),
                               height = 4,
                               width = 10,
                               add_p_values = TRUE,
                               p_label = "p.signif", 
                               method = "wilcox.test", 
                               legend_position = "top", 
                               file_formats = c("pdf", "png"),
                               comparisons = NULL) {
  
  required_packages <- c("ggplot2", "ggsignif", "dplyr", "ggpubr")
  invisible(lapply(required_packages, library, character.only = TRUE))
  
  if (!all(c(x, y) %in% colnames(data))) stop("Columns specified by `x` and `y` must exist in `data`.")
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  if (is.null(colors)) stop("`colors` must be specified as a vector of color values.")

  for (gene in genes) {

    cat("Processing gene:", gene, "\n")
    
    if (add_p_values) {
      comparison_results <- data %>%
        group_by(!!sym(x)) %>%
        summarise(p_value = wilcox.test(!!sym(gene) ~ !!sym(group))$p.value, .groups = "drop")
      
      write.csv(comparison_results, file = file.path(outdir, paste0(gene, "_p_values_table.csv")), row.names = FALSE)
    }
    
    p <- ggplot(data, aes_string(x = x, y = gene, fill = group)) +
      geom_violin(trim = FALSE, scale = "width") +
      geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.shape = NA) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.position = legend_position,
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
      ) +
      scale_fill_manual(values = colors) +
      scale_x_discrete(expand = expansion(mult = c(0.05, 0.1))) +
      labs(title = paste0(" "), 
          #  x = "Cell Type", 
           y = gene) +
      coord_cartesian(ylim = ylim)
    
    if (add_p_values) {
      if (!is.null(comparisons)) {
        p <- p + stat_compare_means(comparisons = comparisons, label = p_label, method = method)
      } else {
        p <- p + stat_compare_means(aes(group = !!sym(group)), label = p_label, method = method)
      }
    }

    save_plot(p, file.path(outdir, paste0(gene, "_violin_plot")), file_formats, height = height, width = width, dpi = 300)
  }
  
  cat("All genes processed successfully.\n")
}



#------------------------------------------------
addGeneExpToMetadata <- function(rds, subColnames = "all",genes = "all", species = "no",RNA = "data", assay = 'RNA', verbose = TRUE) {
    library(Matrix)
    if (!inherits(rds, "Seurat")) {
        stop("The 'rds' object must be a Seurat object.")
    }
    if (!RNA %in% slotNames(rds[[assay]])) {
        stop(paste("The RNA slot", RNA, "is not found in the assay", assay))
    }

    if (species != "no") {
        if (species == "human") {
            genes = upper(genes)
        } else if (species == "mouse") {
            genes = capitalize(genes)
        }
    }

    mtx <- slot(rds[[assay]], RNA)
    rownames_mtx <- rownames(mtx)

    if (genes[1] == "all") {
        genes <- rownames_mtx
        # stop("All genes not supported yet")
        if (verbose) message("Using all genes in the data.")
    } else {
        # Handle regex or specific gene names
        if (any(grepl("\\*|\\[|\\(|\\.", genes))) {
            genes <- unlist(lapply(genes, function(pattern) {
                grep(gsub("\\*", ".*", pattern), rownames_mtx, value = TRUE)
            }))
            if (length(genes) == 0) stop("No genes matched the provided patterns.")
            if (verbose) message(paste("Genes matched by patterns:", paste(genes, collapse = ", ")))
        }
        duplicated_genes <- genes[duplicated(genes)]
        if (length(duplicated_genes) > 0 && verbose) {
            message(paste("The following genes are duplicated and will be removed:", paste(duplicated_genes, collapse = ", ")))
        }
        genes <- unique(genes[genes %in% rownames_mtx])
        if (length(genes) == 0) stop("No valid genes were found in the data.")
    }

    if (subColnames[1] == "all") {
      rds@meta.data <- cbind(rds@meta.data, t(as.data.frame(mtx[genes, , drop = FALSE])))
    } else {
      rds@meta.data <- cbind(rds@meta.data[,subColnames, drop = FALSE], t(as.data.frame(mtx[genes, , drop = FALSE])))
    }

    return(rds)
}



addFreqCount <- function(data, col) {
  if (!col %in% colnames(data)) {
    stop(sprintf("Column '%s' not found in the data frame.", col))
  }
  
  celltype_counts <- table(data[[col]])

  updated_labels <- sapply(as.character(data[[col]]), function(celltype) {
    count <- celltype_counts[celltype]
    paste(celltype, count, sep = "-")
  })

  if (is.factor(data[[col]])) {

    new_levels <- sapply(levels(data[[col]]), function(celltype) {
      count <- celltype_counts[celltype]
      paste(celltype, count, sep = "-")
    })
    data[[col]] <- factor(updated_labels, levels = new_levels)
  } else {

    data[[col]] <- updated_labels
  }
  
  return(data)
}


reorderLevelsByFrequency <- function(factor_var) {
  if(!is.factor(factor_var)) {
    factor_var = factor(factor_var)
  }
  freq <- table(factor_var)
  sorted_levels <- names(sort(freq, decreasing = TRUE))
  factor_var <- factor(factor_var, levels = sorted_levels)
  return(factor_var)
}

reorderClusters <- function(rds, colname) {
  rds@meta.data[[colname]] = reorderLevelsByFrequency(rds@meta.data[[colname]])
  old_clusters <- as.integer(rds@meta.data[[colname]])
  unique_clusters <- sort(unique(old_clusters))
  new_cluster_mapping <- setNames(seq_along(unique_clusters), unique_clusters)
  rds@meta.data[[colname]] <- new_cluster_mapping[as.character(old_clusters)]
  rds@meta.data[[colname]] <- factor(rds@meta.data[[colname]], levels = seq_along(unique_clusters))
  return(rds)
}

dirnameScript <- function(){
    cmd <- commandArgs(trailingOnly = FALSE)
    scriptName <- sub('--file=', "", cmd[grep('^--file=', cmd)])
    path <- system(paste0("realpath ", scriptName), intern=T)
    scriptDirname <- dirname(path)
    return(path)
}

readH5seurat <- function(input,informat = "h5seurat"){
    library(Seurat)
    library(SeuratObject)
    library(SeuratDisk)
    rds <- LoadH5Seurat(input)
    return(rds)
}


check_and_update_seurat_version <- function(rds) {
  seurat_version <- packageVersion("Seurat")
  rds_version <- rds@version
  
  seurat_major_version <- as.numeric(strsplit(as.character(seurat_version), "\\.")[[1]][1])
  rds_major_version <- as.numeric(strsplit(as.character(rds_version), "\\.")[[1]][1])
  
  if (seurat_major_version == rds_major_version) {
    message("The major version of the Seurat object matches the installed Seurat package.")
  } else if (rds_major_version > seurat_major_version) {
    warning("The Seurat object major version is greater than the installed Seurat package major version. Please update the Seurat package.")
  } else {
    message("The Seurat object major version is less than the installed Seurat package major version. Upgrading the Seurat object...")
    rds <- Seurat::UpdateSeuratObject(rds)
  }
  
  return(rds)
}

readSeurat <- function(input){
  library(Seurat)
  library(tools)
  ext <- file_ext(input)
  if (ext == 'rds'){
      rds = readRDS(input)
  }else if(ext == 'h5seurat'){
      rds = readH5seurat(input)
  }else{
    stop('informat rds or h5seurat')
  }
  
  rds = check_and_update_seurat_version(rds)

  # 解决一些Seurat对象没有Idents的问题
  if (any(is.na(Idents(rds)))){
    Idents(rds) = colnames(rds@meta.data)[1]
  }

  return(rds)
}

filterLowCountGroup <- function(rds, groupby, rm_groupby_cells_lessthan) {
  if (!is.null(rm_groupby_cells_lessthan)) {
    if (is.null(groupby)) {
      stop("rm_groupby_cells_lessthan参数需要指定groupby参数")
    }
    
    if (length(rm_groupby_cells_lessthan) == 1) {
      threshold <- as.numeric(rm_groupby_cells_lessthan)

      cross_table <- table(rds@meta.data[[groupby]])
      
      low_groups <- names(cross_table[cross_table < threshold])
      
      if (length(low_groups) > 0) {
        warning(paste0(
          "以下分组的细胞数 <", threshold, "，已删除：", paste(low_groups, collapse = ", ")
        ))
      }
      
      keep_cells <- !rds@meta.data[[groupby]] %in% low_groups
      rds <- rds[, keep_cells]
    }
    
    # 处理双参数情况（组间比较）
    else if (length(rm_groupby_cells_lessthan) == 2) {
      group_col <- rm_groupby_cells_lessthan[1]
      threshold <- as.numeric(rm_groupby_cells_lessthan[2])
      
      if (length(unique(rds@meta.data[[group_col]])) >2 ){
        warning(paste0(
          "分组变量", group_col, "中，组数大于2，不进行细胞数小于",threshold,'的过滤')
        )
        return(rds)
      }

      if (!(group_col %in% colnames(rds@meta.data))) {
        stop(paste0("元数据中不存在列：", group_col))
      }
      
      cross_table <- table(
        rds@meta.data[[group_col]],
        rds@meta.data[[groupby]]
      )
      
      cols_to_remove <- which(
        apply(cross_table, 2, function(col) any(col < threshold))
      )
      
      groupby_levels <- colnames(cross_table)
      
      if (length(cols_to_remove) > 0) {
        removed_groups <- groupby_levels[cols_to_remove]
        warning(paste0(
          "以下分组在组间比较中细胞数 <", threshold, "，已删除：", paste(removed_groups, collapse = ", ")
        ))
        groupby_to_keep <- groupby_levels[-cols_to_remove]
      }else{
        groupby_to_keep <- groupby_levels
      }
      rds <- rds[, rds@meta.data[[groupby]] %in% groupby_to_keep]
    }
    else {
      stop("rm_groupby_cells_lessthan参数必须为单值或长度为2的向量")
    }
  }
  
  return(rds)
}

verifyInputType <- function(input) {

  if (is.character(input) && length(input) == 1) {
    if (file.exists(input) && !dir.exists(input)) {
      return("path")
    } else {
      return("unknown")
    }
  } else if (inherits(input, "Seurat")) {
    return("seurat")
  } else {
    return("unknown")
  }
}

getRDS <- function(input) {
  type <- verifyInputType(input)
  if (type == "path") {
    rds <- readSeurat(input)
    if (!inherits(rds, "Seurat")) {
      stop("文件内容不是 Seurat 对象")
    }
    return(rds)
  } else if (type == "seurat") {
    return(input)
  } else {
    stop("输入类型无效")
  }
}

saveSeurat = function(rds, outdir, prefix = NULL,groupby=NULL){
    createDir(outdir)
    rds = loadRDS(rds,groupby = groupby)
    if (is.null(prefix)){
      seurat_version <- packageVersion("Seurat")
      seurat_major_version <- as.numeric(strsplit(as.character(seurat_version), "\\.")[[1]][1])
      prefix = paste0("_v",seurat_major_version)
    }
    saveRDS(rds,paste0(outdir,'/','data_ob',prefix,'.rds'))
}


loadRDS <- function(rds_filepath,
                    subnew_celltype = 'all',
                    subsampleid = 'all',
                    subgroup = 'all',
                    subcluster = 'all',
                    predicate = 'all',
                    groupby_levels = 'no',
                    downsample = NULL,
                    groupby = NULL,
                    rm_groupby_cells_lessthan = NULL,
                    palette = NULL,
                    combo = NULL) {
  # 读取RDS，验证输入，可输入路径或Seurat对象
  rds = getRDS(rds_filepath)

  # 为指定组标准化指定的levels
  if (groupby_levels != 'no'){
    rds = standardizeLevels(rds, groupby, groupby_levels)
  }

  # 生成combo分组变量
  if (!is.null(combo) && length(standardizeDelimitedVector(combo)) == 2) {
    combo = standardizeDelimitedVector(combo)
    generate_combined(rds@meta.data[, combo[1]],rds@meta.data[, combo[2]]) -> rds@meta.data[, paste0(combo[1], "_", combo[2])]
  }

  # 获取RDS的元数据列名
  use_col = unique(c("new_celltype",'sampleid','group','clusters',getPredicateVars(predicate),groupby))


  
  # 标准化颜色
  rds = standardizeColors(rds,use_col,palette = palette)

  # 根据条件对RDS进行子集化
  rds = subsetRDS(rds,
      new_celltype = subnew_celltype,
      sampleid = subsampleid,
      group = subgroup,
      cluster = subcluster,
      predicate = predicate
    )

  # 根据分组变量对RDS进行降采样
  if (!is.null(downsample)){
    rds = getDownsampleRDS(rds,groupby,downsample)$optimal_subset
  }
  
  # 解决分组变量中的元素存在细胞数小于指定数的问题
  rds = filterLowCountGroup(rds, groupby, rm_groupby_cells_lessthan)

  # 标准化常见及指定分组变量的levels
  for (Colname in use_col){
    if (Colname %in% colnames(rds@meta.data)){
      rds = standardizeLevels(rds,Colname)
    }
  }



  return(rds)
}

generate_combined <- function(vec1_fac, vec2_fac) {
  # 检查输入是否为因子类型
  if (!is.factor(vec1_fac) || !is.factor(vec2_fac)) {
    stop("输入的两个向量必须都是因子（factor）类型")
  }
  
  # 获取两个因子的levels
  lev1 <- levels(vec1_fac)
  lev2 <- levels(vec2_fac)
  
  # 生成组合因子的levels，按第一个因子的levels顺序，每个level下再按第二个因子的levels顺序
  combo_levels <- unlist(lapply(lev1, function(l1) paste0(l1, "_", lev2)))
  
  # 生成原始组合的字符向量
  combo_chars <- paste0(as.character(vec1_fac), "_", as.character(vec2_fac))
  
  # 将组合字符转换为因子，并设置正确的levels顺序
  combo_factor <- factor(combo_chars, levels = combo_levels)
  
  return(combo_factor)
}

is2GroupEquality <- function(rds, group_col1 = "group", group_col2 = "group_col") {
  if (!inherits(rds, "Seurat")) {
    stop("Input 'rds' must be a Seurat object.")
  }
  
  meta_colnames <- colnames(rds@meta.data)
  
  stopifnot(
    group_col1 %in% meta_colnames,
    group_col2 %in% meta_colnames
  )
  
  value1 <- length(unique(rds@meta.data[[group_col1]]))
  value2 <- length(unique(rds@meta.data[[group_col2]]))
  value3 <- nrow(getGroupOrder(rds, c(group_col1, group_col2)))

  if (value1 != value2) message(sprintf("%s 和 %s 唯一值数量不等", group_col1, group_col2))
  if (value2 != value3) message(sprintf("%s 和分组组合数不等", group_col2))
  if (value1 != value3) message(sprintf("%s 和分组组合数不等", group_col1))

  return(length(unique(c(value1, value2, value3))) == 1)
}

getColorOrder = function(rds,useCol,palette = "customecol2"){
    library(dplyr)
    if (!paste0(useCol,"_col") %in% colnames(rds@meta.data)) {
        rds@meta.data = addColorToDataFrame(rds@meta.data, useCol,palette)
    }
    colors_data <- rds@meta.data[!duplicated(rds@meta.data[,useCol]), c(useCol, paste0(useCol,"_col"))]
    colors_data[,useCol] = factor(colors_data[,useCol],levels = levels(rds@meta.data[,useCol]))
    colors_data = colors_data %>% arrange(!!sym(useCol))
    colors2use = colors_data[,paste0(useCol,"_col")]
    print(colors_data[,useCol])
    return(colors2use)
}

getGroupOrder = function(rds,useCol){
    return(rds@meta.data[!duplicated(rds@meta.data[,useCol]), ][,useCol])
}

getMeanMtx <- function(rds, group, genes = "all", species = "no", order_group = "no", RNA = "data", assay = 'RNA', scale = FALSE, stat = "mean") {
    library(Matrix)
    library(ComplexHeatmap)

    if (species != "no") {
        if (species == "human") {
            genes = upper(genes)
        } else if (species == "mouse") {
            genes = capitalize(genes)
        }
    }

    mtx = slot(rds[[assay]], RNA)
    group = rds@meta.data[[group]]

    if (genes[1] != "all") {
        duplicated_genes = genes[duplicated(genes)]

        if (length(duplicated_genes) > 0) {
            print(paste("The following genes are duplicated in the gene list (will be deleted):", paste(duplicated_genes, collapse = ", ")))
        }

        genes = genes[!duplicated(genes)]

        filtergenes = genes[!genes %in% rownames(mtx)]
        if (length(genes) == 0) {
            stop("No genes in the gene list are found in the data")
        }
        genes = genes[genes %in% rownames(mtx)]

        if (length(filtergenes) > 0) {
            print(paste("The following genes are not found in the data:", paste(filtergenes, collapse = ", ")))
        }
    } else {
        genes = rownames(mtx)
    }

    mtx = mtx[genes, , drop = FALSE]

    colnames(mtx) <- as.character(group)
    unique_cols <- unique(colnames(mtx))

    result_mtx <- Matrix(0, nrow = nrow(mtx), ncol = length(unique_cols), 
                         dimnames = list(rownames(mtx), unique_cols), 
                         sparse = TRUE)

    for (i in seq_along(unique_cols)) {
        cols <- which(colnames(mtx) == unique_cols[i])
        if (stat == "mean") {
            result_mtx[, i] <- rowMeans(as.matrix(mtx[, cols, drop = FALSE]), na.rm = TRUE)
        } else if (stat == "median") {
            result_mtx[, i] <- apply(as.matrix(mtx[, cols, drop = FALSE]), 1, median, na.rm = TRUE)
        } else {
            stop("Invalid value for 'stat'. Use 'mean' or 'median'.")
        }
    }

    result_df <- as.data.frame(result_mtx)

    if (scale) {
        result_df <- pheatmap:::scale_rows(log2(result_df + 0.0001))
    }

    if (order_group != "no") {
        result_df <- result_df[, order_group]
    } else {
        result_df <- result_df[, order(colnames(result_df))]
    }

    return(result_df)
}

getClusteringResult <- function(rds, outdir, reduct1 = "pca", reduct2 = "umap", saveRDS = TRUE,resolution = "None", palette = "customecol2", pointsize = 0.5){
    library(dplyr)
    library(ggplot2)
    library(Seurat)
    library(SummarizedExperiment)
    dim_outdir = file.path(outdir,paste0(reduct1, "_Dimension_Reduction"))
    if ( !dir.exists(dim_outdir) ){
        dir.create(dim_outdir, recursive = T)
    }
    reduct1_coord = FetchData(rds, vars = c("rawbc", paste0( Seurat::Key(rds)[reduct1], 1:2))) %>%
                    dplyr::rename( "Barcode" = "rawbc")
    write.table( reduct1_coord, file.path(dim_outdir, paste0(reduct1, "_Dimension_Reduction_coordination.csv")),
                sep = ",", col.names = T, row.names = F, quote = F)
    dim_outdir = file.path(outdir,paste0(reduct2, "_Dimension_Reduction"))
    if ( !dir.exists(dim_outdir) ){
        dir.create(dim_outdir, recursive = T)
    }
    reduct2_coord = FetchData(rds,
                        vars = c("rawbc", paste0( Seurat::Key(rds)[reduct2], 1:2))) %>%
                        dplyr::rename( "Barcode" = "rawbc")

    write.table( reduct2_coord, file.path(dim_outdir, paste0(reduct2, "_Dimension_Reduction_coordination.csv")),
    sep = ",", col.names = T, row.names = F, quote = F)

    rds = reorderClusters(rds, "clusters")

    Idents(rds) = rds@meta.data$clusters
    cell_id = Seurat::Idents(rds)
    new_id = as.numeric(as.vector(cell_id)) # new cluster id start from 1
    names(new_id) = names(cell_id)
    new_id = as.factor(new_id)
    cell_count_by_cluster = table( new_id )
    Seurat::Idents(rds) = new_id
    cluster_result_colname = paste(Seurat::DefaultAssay(rds), reduct2,"res", resolution, sep = ".")
    rds[[cluster_result_colname]] = Seurat::Idents(rds)
    rds[["clusters"]] = Seurat::Idents(rds)
    reordered_cell_count_by_cluster = table(Seurat::Idents(rds))
    cell_count_labels = paste(paste(names(reordered_cell_count_by_cluster),reordered_cell_count_by_cluster,sep="-")," cells")
    ggdimplot = Seurat::DimPlot(object = rds,reduction = reduct2, dims = c(1,2),
                label = T, group.by= cluster_result_colname,   pt.size = pointsize)
    custom_cols =  getColorOrder(rds,'clusters')

    ggdimplot = ggdimplot + ggplot2::labs(title = "") + theme(aspect.ratio = 1/1)+
            ggplot2::scale_colour_manual( values = custom_cols,
                                    breaks=levels(Seurat::Idents(rds)),
                                    labels = cell_count_labels)
    ggplot2::ggsave(file.path(dim_outdir,paste0(reduct2, "_groupby_cluster","_resolution", resolution,"_plot.pdf",collapse="")), 
                plot = ggdimplot, bg="white", width = ifelse(length(unique(Seurat::Idents(rds))) <= 20,7.5,9))
    ggplot2::ggsave(file.path(dim_outdir,paste0(reduct2, "_groupby_cluster","_resolution", resolution,"_plot.png",collapse="")),dpi=1000, 
                plot = ggdimplot, bg="white", width = ifelse(length(unique(Seurat::Idents(rds))) <= 20,7.5,9))

    clustering_df = rds@meta.data %>% 
    dplyr::rename(Barcode = "rawbc") %>% 
    dplyr::select(Barcode, sampleid, clusters, group)
    write.table(clustering_df, quote = F,sep =",",row.names = F,
                file.path(outdir,paste0("clustering_results.csv",collapse = "")))
                
    if (saveRDS){
        saveRDS(rds, paste0(outdir,'/data_ob.rds'))
    }
}

getManualannoResult <- function(rds,outdir,useCol="new_celltype",addFreqCount="no", saveRDS = 'yes',label = FALSE, legend = 'yes',orderbyFreq = "yes",reduct = "umap",palette="customecol2",pointsize=0.5){
    source("/home/chenglong.liu/RaD/myscript/Get_colors.R")
    library(dplyr)
    library(ggplot2)
    library(stringr)
    
    if(!file.exists(outdir)){
        dir.create(outdir,recursive = T)
    }
    
    if (addFreqCount != "no"){
      rds@meta.data$new_celltype_withFreqCount = rds@meta.data[,useCol]
      rds@meta.data = add_FreqCount(rds@meta.data,"new_celltype_withFreqCount")
      useCol = "new_celltype_withFreqCount"
    }
    if (orderbyFreq =="yes"){
      rds@meta.data[,useCol] = reorderLevelsByFrequency(rds@meta.data[,useCol])
    }else if(orderbyFreq =="no"){
      print("No reorder")
    }else{
      rds@meta.data[,useCol] = as.character(rds@meta.data[,useCol])
      orderbyFreq = unlist(str_split(orderbyFreq,","))
      rds@meta.data[,useCol] = factor(rds@meta.data[,useCol],levels = orderbyFreq)
    }
    
    if(paste0(useCol,"_col") %in% colnames(rds@meta.data)){
      if (palette != "yes"){
        rds@meta.data <- rds@meta.data[ , !(names(rds@meta.data) %in% c(paste0(useCol,"_col")))]
        get_colors_result = get_colors(rds@meta.data,useCol,palette)
        rds@meta.data[,paste0(useCol,"_col")] = get_colors_result[["object_meta"]][,paste0(useCol,"_col")]
      }
    }



    colors2use = getColorOrder(rds,useCol)

    gg = Seurat::DimPlot(object = rds, reduction = reduct,pt.size = pointsize,group.by=useCol,label=label)+
        ggplot2::theme( plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::scale_colour_manual( values = unname(colors2use)) + theme(aspect.ratio = 1/1)
    if (legend != 'yes'){
        gg = gg + NoLegend()
    }

    width_value <- max(nchar(as.vector(unique(rds@meta.data[,useCol])))) / 15 + 7

    # 构建完整文件路径
    file_path <- file.path(outdir, paste0(useCol))

    # 使用 ggsave 保存
    ggsave(
      filename = paste0(file_path,'.png'),  # 文件路径
      plot = gg,             # ggplot 对象
      width = width_value,   # 宽度（英寸）
      height = width_value * 0.75,  # 高度，通常设置为宽度的0.75倍
      dpi = 300,             # 分辨率
      limitsize = FALSE,     # 不限制图像大小
      bg = "white"           # 背景颜色
    )

    # 使用 ggsave 保存
    ggsave(
      filename = paste0(file_path,'.pdf'),  # 文件路径
      plot = gg,             # ggplot 对象
      width = width_value,   # 宽度（英寸）
      height = width_value * 0.75,  # 高度，通常设置为宽度的0.75倍
      dpi = 300,             # 分辨率
      limitsize = FALSE,     # 不限制图像大小
      bg = "white"           # 背景颜色
    )


    simplified_meta = rds@meta.data %>%
                            dplyr::rename( "Barcode" = "rawbc") %>%
                                    dplyr::select( Barcode, sampleid, clusters, group,!!useCol)
    write.table(simplified_meta, quote = F,sep =",",row.names = F,
            file.path(outdir,paste0(useCol,".metadata.csv",collapse = "")))

    if (saveRDS=='yes'){
        saveRDS(rds, paste0(outdir,'/data_ob.rds'))
    }
}


#-----------------homologene_transformed-------------------------
homologene_seurat_replace <- function(data_ob,assay2use,inTaxidgene,outTaxidgene){
    counts=data_ob@assays[[assay2use]]@counts
    counts_filted = counts[as.vector(inTaxidgene),]
    rownames(counts_filted) = as.character(outTaxidgene)
    data_ob@assays[[assay2use]]@counts=counts_filted
    expr=data_ob@assays[[assay2use]]@data
    gene=rownames(expr)
    expr_filtered = expr[as.vector(inTaxidgene),]
    rownames(expr_filtered) = as.character(outTaxidgene)
    data_ob@assays[[assay2use]]@data=expr_filtered

    print("Since the change is changed to scale only for high-variable genes, scale.data is discarded from the homologous converted RDS")
    features=data_ob@assays[[assay2use]]@meta.features
    features_filtered = features[as.vector(inTaxidgene),]
    features_filtered = as.matrix(features_filtered)
    rownames(features_filtered) = as.character(outTaxidgene)
    features_filtered = as.data.frame(features_filtered)
    data_ob@assays[[assay2use]]@meta.features=features_filtered

    return(data_ob)
}

homologene_transformed = function(data_ob,inTaxid,outTaxid,output,assay="RNA",genelist="no"){
    suppressPackageStartupMessages(library("Seurat"))
    suppressPackageStartupMessages(library("homologene"))
    assay2use = assay
    counts=data_ob@assays[[assay2use]]@counts
    if(genelist == "no"){
      gene=row.names(counts)
      in2out = homologene(gene,inTax=inTaxid,outTax=outTaxid)
      index <- duplicated(in2out[,1])
      in2out <-in2out[!index,]
      index<-duplicated(in2out[,2])
      in2out <- in2out[!index,]
      homologene <- in2out[,1:2]
      write.table(homologene,file.path(output,paste0("homologene.",inTaxid,"_2_",outTaxid,".xls",collapse=".")),sep="\t",row.names=F,quote=F)
      data_ob <- homologene_seurat_replace(data_ob,assay2use,in2out[,1],in2out[,2])
    }else{
      in2out <- read.delim(genelist,sep="\t",header=T)
      index <- duplicated(in2out[,1])
      in2out <-in2out[!index,]
      index<-duplicated(in2out[,2])
      in2out <- in2out[!index,]
      final_gene = in2out[,1][in2out[,1]  %in% row.names(counts)]
      in2out <- in2out[!is.na(in2out[,1]), ] # 解决报错missing values in 'row.names' are not allowed
      rownames(in2out)=as.vector(in2out[,1])
      in2out <- in2out[as.vector(final_gene),]
      data_ob <- homologene_seurat_replace(data_ob,assay2use,in2out[,1],in2out[,2])
    }
    return(data_ob)
}

# 定义替代函数 - 专门用于读取GTF并提取gene_id和gene_name
import_gtf_simple <- function(gtf_file, feature_type = "exon") {
  library(data.table)
  # 检查文件是否存在
  if (!file.exists(gtf_file)) {
    stop("File not found: ", gtf_file)
  }
  
  # 确定是否是压缩文件
  is_gz <- grepl("\\.gz$", gtf_file)
  
  # 构建读取命令
  if (is_gz) {
    cmd <- paste("zcat", gtf_file)
  } else {
    cmd <- paste("cat", gtf_file)
  }
  
  # 添加过滤注释行的命令
  cmd <- paste(cmd, "| grep -v '^#'")
  
  # 使用fread读取
  gtf <- fread(cmd = cmd, sep = "\t", header = FALSE, 
               quote = "", fill = TRUE, na.strings = ".")
  
  # 设置列名（标准GTF格式）
  colnames(gtf) <- c("seqname", "source", "feature", "start", "end", 
                     "score", "strand", "frame", "attributes")
  
  # 将start和end转换为数值
  gtf[, start := as.numeric(start)]
  gtf[, end := as.numeric(end)]
  
  # 过滤指定的feature类型
  if (!is.null(feature_type)) {
    gtf <- gtf[feature %in% feature_type, ]
  }
  
  # 解析attributes列 - 简化版本，只提取gene_id和gene_name
  parse_gene_info <- function(attr) {
    gene_id <- NA_character_
    gene_name <- NA_character_
    
    # 查找gene_id
    gene_id_match <- regmatches(attr, regexpr('gene_id\\s+"[^"]+"', attr))
    if (length(gene_id_match) > 0) {
      gene_id <- gsub('gene_id\\s+"|"', '', gene_id_match[1])
    }
    
    # 查找gene_name
    gene_name_match <- regmatches(attr, regexpr('gene_name\\s+"[^"]+"', attr))
    if (length(gene_name_match) > 0) {
      gene_name <- gsub('gene_name\\s+"|"', '', gene_name_match[1])
    }
    
    # 如果gene_name不存在，尝试其他可能的字段
    if (is.na(gene_name)) {
      # 尝试gene_symbol
      gene_symbol_match <- regmatches(attr, regexpr('gene_symbol\\s+"[^"]+"', attr))
      if (length(gene_symbol_match) > 0) {
        gene_name <- gsub('gene_symbol\\s+"|"', '', gene_symbol_match[1])
      }
    }
    
    return(list(gene_id = gene_id, gene_name = gene_name))
  }
  
  # 解析所有行的gene信息
  gene_info <- lapply(gtf$attributes, parse_gene_info)
  
  # 添加到数据框
  gtf[, gene_id := sapply(gene_info, function(x) x$gene_id)]
  gtf[, gene_name := sapply(gene_info, function(x) x$gene_name)]
  
  # 转换为data.frame并返回
  return(as.data.frame(gtf))
}
#-----------------mix_Species_Split-------------------------
# Helper function to read and preprocess GEM classification data
read_and_preprocess_gem_classification <- function(sample, input_dir) {
  data_path <- file.path(input_dir, sample, "outs", "analysis", "gem_classification.csv")
  if (!file.exists(data_path)) {
    stop(paste("File not found:", data_path))
  }
  data <- read.csv(data_path)
  data$cells <- paste0(sample, "-", data$barcode)
  data$cells <- gsub("-1", "", data$cells)
  return(data)
}

# Helper function to filter barcodes of interest
filter_barcodes <- function(keyData, mixSpecies) {
  call_type <- mixSpecies
  return(subset(keyData, call == call_type)$cells)
}

# Helper function to process Seurat object
process_seurat_object <- function(subrds, mixSpecies) {
  genes_of_interest <- rownames(subrds)[grep(mixSpecies, rownames(subrds))]
  
  subrds <- subset(subrds, features = genes_of_interest)

  # Clean up gene names in the assays
  rownames(subrds@assays$RNA@counts) <- gsub(paste0(mixSpecies, "-+"), "", rownames(subrds@assays$RNA@counts))

  if (!is.null(subrds@assays$RNA@data)) {
  rownames(subrds@assays$RNA@data) <- gsub(paste0(mixSpecies, "-+"), "", rownames(subrds@assays$RNA@data))
  }

  if (!is.null(subrds@assays$RNA@scale.data)) {
  rownames(subrds@assays$RNA@scale.data) <- gsub(paste0(mixSpecies, "-+"), "", rownames(subrds@assays$RNA@scale.data))
  }

  return(subrds)
}

mix_Species_Split <- function(input_dir, data_ob, mixSpecies) {
  # Main processing flow
  # mouse = "mm10", human = "GRCh38", rat = "Rat"
  samples <- list.dirs(input_dir, full.names = FALSE, recursive = FALSE)
  gem_classification_list <- lapply(samples, read_and_preprocess_gem_classification, input_dir = input_dir)
  keyData <- do.call(rbind, gem_classification_list)
  
  barcodes_of_interest <- filter_barcodes(keyData, mixSpecies)
  
  subrds <- subset(data_ob, cells = barcodes_of_interest)
  subrds <- process_seurat_object(subrds, mixSpecies)
  
  return(subrds)
}

# ----------------------
create_group_columns <- function(data, cluster_col_name) {

  unique_clusters <- sort(unique(data[[cluster_col_name]]))
  
  for (cluster in unique_clusters) {

    group_col_name <- paste0("group_", cluster)
    
    data[[group_col_name]] <- ifelse(data[[cluster_col_name]] == cluster,
                                     as.character(cluster),
                                     paste0(setdiff(unique_clusters, cluster), collapse = "_"))
  }
  
  return(data)
}
