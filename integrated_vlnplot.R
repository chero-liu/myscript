###############################
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# library(ggbeeswarm) # for ggbeeswarm mode
make_vlnplot_data <- function(seurat_obj, features, groupby = "clusters", get.color = TRUE, palette = "customecol2") {
    require(tidyr)
    # 公共数据准备部分
    vln_data = Seurat::FetchData(seurat_obj, c(groupby,features))
    data_long <- gather(vln_data, key = 'genes',value = 'expression', features )
    if(get.color){
      # 获取颜色映射
      color_use <- get_color_use(seurat_obj, groupby, palette = palette)
      user_color_pal <- color_use[["new_celltype_pal"]]
      data_long[,paste0(groupby,"_col")] = user_color_pal[data_long[, groupby, drop = T]]
    }
    return(data)
}
##############################
vlnplot_integrated  = function(features, seurat_obj = NULL, vln_data = NULL, 
        groupby =  "clusters", 
        mode = c("line_black", "light_box", "beeswarm", "black_box"), 
        color_use = NULL,
        palette = "customecol2",
        show_box = TRUE,
        return_data = FALSE,
        order = NULL,
        max_str = NULL,
        noise = FALSE # 默认不加噪点
        ){ 
  require(tidyr)
  if( ! is.null(vln_data)){
    if( all( c("genes",'expression') %in% colnames(vln_data))){
      data_long = vln_data
    } else if (any(features %in% colnames(vln_data))) {
      vln_data = subset(vln_data, select = c(groupby,features))
      data_long = gather(vln_data, key = 'genes',value = 'expression', features )
    }
    col = get_colors(data_long, groupby, palette = palette)[["new_celltype_pal"]]
  }else if( ! is.null(seurat_obj)){
    vln_data = Seurat::FetchData(seurat_obj, c(groupby,features))
    data_long <- gather(vln_data, key = 'genes',value = 'expression', features )
    col = get_colors(seurat_obj@meta.data, groupby, palette = palette)[["new_celltype_pal"]]
  }else{
    stop("Please input seurat_obj or vln_data")
  }
  if(is.null(order)){
      data_long[, groupby] = factor(data_long[,groupby], levels = names(col))
  }else{
      data_long[, groupby] = factor(data_long[,groupby], levels = order)
  }
  set.seed(seed = 42)
  if(noise){
    data_long$expression.ori = data_long$expression
    noise <- rnorm(n = length(x = data_long$expression)) / 100000
    data_long$expression <- data_long$expression.ori  + noise
  }
  # 公共颜色设置
  if (is.null(color_use)) {
      color_use <- col
  }
  if( !is.null(max_str)){
     data_long$gene_padded <- sprintf(paste0("%-", max_str, "s"), data_long$genes)
     data_long$gene_padded = factor(data_long$gene_padded, levels =  sprintf(paste0("%-", max_str, "s"), features))
  } else {
     data_long$gene_padded <- data_long$genes
  }
  # 构建图形
  p <- switch(mode,
    "line_black" = {
      ggplot(data_long, aes_string(groupby, "expression", fill = groupby, color = groupby)) +
        geom_violin(cex = 0.6, scale = "width", width = 0.8, color = "black", linewidth = 0.4) +
        geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white", size = 0.2) +
        theme_classic() +
        theme(
          panel.border = element_rect(colour = "black", fill = NA, size = 0.8),
          panel.spacing = unit(0, "lines")
        )
    },
    "light_box" = {
      ggplot(data_long, aes_string(groupby, "expression", fill = groupby, color = groupby)) +
        geom_violin(cex = 0.6, scale = "width", width = 0.8, color = "transparent", linewidth = 0.4) +
        geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white", size = 0.2) +
        theme_classic() +
        theme(
          panel.border = element_rect(colour = "black", fill = NA, size = 0.8),
          panel.spacing = unit(0, "lines")
        )
    },
    "beeswarm" = {
      require(ggbeeswarm)
      ggplot(data_long, aes_string(groupby, "expression", fill = groupby, color = groupby)) +
        geom_violin(cex = 0.6, scale = "width", width = 0.8, color = "black", fill = "white") +
        geom_quasirandom(shape = 21, size = 0.4, color = "transparent",
                          position = position_jitter(seed = 123, width = 0.1, height = 0)) +
        theme(
          strip.text.y = element_text(size = 15, color = "black",angle = 0, hjust = 0), axis.text.x = element_text(angle = 0)
        )
    },
    "black_box" = {
      ggplot(data_long, aes_string(groupby, "expression", fill = groupby, color = groupby)) +
        geom_violin(cex = 0.6, scale = "width", width = 0.8, color = "black", linewidth = 0.4, alpha = 0.6) +
        geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "black", size = 0.2) +
        stat_summary(fun = median, geom = "point", shape = 21, size = 0.3, color = "white", fill = "white") +
        theme_classic() +
        theme(
          panel.border = element_rect(colour = "black", fill = NA, size = 0.8),
          panel.spacing = unit(0, "lines")
        )
    }
  )
  if( !show_box ){
    p = remove_geoms(p,"GeomBoxplot")
    if(mode == "alpha"){p = remove_geoms(p,"GeomPoint")}
  }
  p = p + facet_grid(gene_padded ~ ., scales = "free_y") +
          scale_y_continuous(labels = function(x) ifelse(x == 0, "", sprintf("%.1f", x))) +
          scale_fill_manual(values = color_use) +
          theme(
          legend.position = 'none',
            axis.text.x = element_text(size = 28, angle = 0, color = "black"),
            axis.text.y = element_text(size = 18, angle = 0, color = "black"),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            strip.text.y = element_text(size = 28, color = "black", angle = 0, hjust = 0),
            strip.background = element_blank()
          )
  if( max(nchar(as.character(data_long[,groupby]))) > 3){
      # 获取 X 轴第一个标签文本
      label_f <- as.character(levels(p$data[,rlang::as_name(p$mapping$x)])[1])
      # 计算文本宽度（单位：pt）
      grob_width <- grid::grobWidth(grid::textGrob(label_f, gp = grid::gpar(fontsize = 28)))
      label_widths <- grid::convertWidth(grob_width, "pt", valueOnly = TRUE)

      # 获取最大宽度
      max_width <- ifelse(nchar(label_f) < 3 ,5,max(label_widths))
      p = p + theme(axis.text.x = element_text(size = 28, color = "black", angle = 45, vjust = 1, hjust = 1),
                  plot.margin = margin( t = 5,b = 5, l = max_width, unit = "pt")) 
  }

  if(return_data){
    return(list(data = data_long, plot = p))
  }else{
    return(p)
  }
}


remove_geoms <- function(x, geom_type, last_only = T) {
  # Find layers that match the requested type.
  selector <- sapply(x$layers,
                     function(y) {
                       class(y$geom)[1] == geom_type
                     })
  if(last_only) 
    selector <- max(which(selector))
  # Delete the layers.
  x$layers[selector] <- NULL
  x
}

