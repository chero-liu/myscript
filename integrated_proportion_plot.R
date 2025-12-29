make_barplot_data <- function(seurat_obj,  x_var = "clusters", fill_var = "sampleid",
     get.color = TRUE, palette = "tableau20") {
      # 计算频率数据
      freq_data <- seurat_obj@meta.data %>%
        dplyr::select(all_of(c(x_var, fill_var))) %>%
        dplyr::group_by(across(all_of(c(x_var, fill_var)))) %>%
        dplyr::summarise(cell_number = n(), .groups = "drop_last") %>%
        dplyr::mutate(freq = (cell_number / sum(cell_number)) * 100) %>%
        ungroup()
      if(get.color){
        # 获取颜色映射
        color_use <- get_color_use(seurat_obj, fill_var, palette = palette)
        user_color_pal <- color_use[["new_celltype_pal"]]
        freq_data[,paste0(fill_var,"_col")] = user_color_pal[freq_data[, fill_var, drop = T]]
      }
    return(freq_data)
}

proportion_plot_integrated <- function(
    freq_data = NULL,
    seurat_obj = NULL,
    x_var = "clusters",
    fill_var = "sampleid",
    plot_type = "default_barplot",
    palette = "tableau20",
    output_prefix = "plot",
    output_dir = ".", 
    base_height = 7,
    base_width = 4,
    angle = 90,
    text_size = 28) {
  
  # 检查必要包
  require(ggplot2)
  require(dplyr)
  if (plot_type == "alluvium_barplot") {
    if (!requireNamespace("ggalluvial", quietly = TRUE)) {
      stop("ggalluvial package is required for alluvium_barplot")
    }
    library(ggalluvial, lib.loc = c("/home/ziqingzhen/R/x86_64-conda-linux-gnu-library/4.0"))
  }
  
  if (plot_type %in% c("solid_barplot", "hollow_barplot")) {
    if (!requireNamespace("ggnewscale", quietly = TRUE)) {
      stop("ggnewscale package is required for solid_barplot and hollow_barplot")
    }
    library(ggnewscale, lib.loc = "/data/software/scrna_software_bak/miniconda3/envs/OESingleCell_3.00_visium_produce/lib/R/library")
  }  
  if( ! is.null(seurat_obj)){
      freq_data <- make_barplot_data(seurat_obj, x_var = x_var, fill_var = fill_var, get.color = TRUE, palette = palette)
  }else if( is.null(freq_data)){
      stop("Please input seurat_obj or plot_data")
  }
    # 获取颜色映射
  color_use <- get_colors(freq_data, fill_var, palette = palette)
  user_color_pal <- color_use[["new_celltype_pal"]]
  
# 计算基础参数
num_categories <- length(unique(freq_data[[x_var]]))
max_x_char <- max(nchar(as.character(unique(freq_data[[x_var]]))))
legend_items <- unique(freq_data[[fill_var]])
n_legend_items <- length(legend_items)
max_legend_char <- if (n_legend_items > 0) max(nchar(as.character(legend_items))) else 0
legend_cols <- case_when(
  n_legend_items > 30 ~ 3,
  n_legend_items > 15 ~ 2,
  TRUE ~ 1
)

# 调整基础宽度（对≤3组的特殊处理）
if (fill_var == "clusters" && legend_cols >= 2) {
  # 当fill_var是clusters时，基础宽度固定为2
  base_width_modified <- 1
} else {
  # 其他情况使用原有规则
  base_width_modified <- if (num_categories <= 3) 4 else base_width
}

# 高度和宽度调整函数
get_height_add <- function(chars) case_when(
  chars > 25 ~ 10, chars > 15 ~ 5, chars > 10 ~ 1, TRUE ~ 0
)

get_width_add <- function(chars) case_when(
  chars > 30 ~ 5, chars > 20 ~ 3, chars > 10 ~ 1, TRUE ~ 0
)

# 应用基础计算
adjusted_height <- base_height + get_height_add(max_x_char)
width_add_val <- get_width_add(max_x_char)
per_cluster_width <- 0.5

# 应用图形类型基础宽度计算
base_width_value <- if(plot_type %in% c("solid_barplot", "hollow_barplot")) 8 else base_width_modified
adjusted_width <- base_width_value + width_add_val + num_categories * per_cluster_width

# 特殊高度调整（环状图）
if(plot_type %in% c("solid_barplot", "hollow_barplot")) {
  adjusted_height <- max(base_width, base_width + num_categories*0.2)
}

# 处理coord_flip特殊规则
if (plot_type == "coord_flip_barplot") {
  legend_add <- case_when(
    legend_cols > 1 && max_legend_char > 15 ~ 7,
    legend_cols == 1 && max_legend_char > 20 ~ 6,
    legend_cols == 1 && max_legend_char > 10 ~ 3,
    TRUE ~ 0
  )
  xlabel_add <- ifelse(max_x_char > 15, 4, 0)
  cols_add <- ifelse(legend_cols > 2, 4, 0)
  adjusted_width <- 8 + legend_add + xlabel_add + cols_add
} else {
  # 其他图形类型补偿规则
  extra_width <- case_when(
    legend_cols > 1 && max_legend_char > 25 ~ 16,
    legend_cols > 1 && max_legend_char > 15 ~ 8,
    legend_cols == 1 && max_legend_char > 20 ~ 5,
    legend_cols == 1 && max_legend_char > 10 ~ 2,
    legend_cols > 1 ~ 4,
    TRUE ~ 0
  )
  
  adjusted_width <- adjusted_width + extra_width
  
  # 高度补偿
  if(legend_cols > 1 && max_legend_char <= 25) {
    adjusted_height <- adjusted_height + ifelse(max_legend_char > 15, 0, 2)
  } else if(legend_cols == 1 && max_legend_char > 20) {
    adjusted_height <- adjusted_height + 2
  }
}

#base_aes <- aes(x = !!sym(x_var), y = freq, fill = !!sym(fill_var))

  # 生成图形
if (plot_type == "default_barplot") {
  final_plot <- freq_data %>%
    ggplot(aes(x = .data[[x_var]], y = freq, fill = .data[[fill_var]])) +
    geom_col(width = 0.9) +
    scale_fill_manual(values = user_color_pal) +
    scale_x_discrete(name = "") +
    scale_y_continuous(
      expand = c(0, 0),
      name = "Proportion[%]"
    ) +
    theme_minimal() +
    theme(
      # --------------------------
      # 1. 强制图例符号宽高一致（解决高度不变问题）
      # --------------------------
      legend.key.width = unit(1.2, "lines"),    # 图例符号宽度（按需调整，例：1.2 lines）
      legend.key.height = unit(1.2, "lines"),   # 图例符号高度（与宽度完全一致，确保正方形）
      legend.key.size = NULL,                   # 显式清空该参数，避免覆盖宽高设置
      
      # --------------------------
      # 2. 增大图例与图片顶部的间距（避免覆盖）
      # --------------------------
      # 方法1：增大图例自身的顶部边缘空白（legend.margin）
      legend.margin = margin(
        t = 10,  # 顶部留白从2pt增至10pt（关键！拉开图例与图片顶部的距离）
        r = 5, 
        b = 5, 
        l = 5, 
        unit = "pt"  # 用pt更精细控制，10pt约等于3.5mm
      ),
      # 方法2：（可选）增大图片整体的顶部边距（双重保障）
      plot.margin = margin(
        t = 30,  # 图片顶部边距从20pt增至30pt
        r = 15, 
        b = 15, 
        l = 15, 
        unit = "pt"
      ),
      
      # --------------------------
      # 3. 其他原有配置（保留并微调）
      # --------------------------
      legend.spacing.y = unit(0.15, "lines"),  # 图例内部条目间距（不变）
      legend.text = element_text(
        size = 22,
        margin = margin(t = 0.2, b = 0.2)
      ),
      legend.title = element_text(size = 22),
      text = element_text(family = "sans"),
      plot.background = element_rect(fill = 'white', color = 'white'),
      panel.grid = element_blank(),
      axis.title.y = element_text(size = text_size),
      axis.text.x = element_text(
        size = text_size, color = "black",
        angle = angle, vjust = 0.5, hjust = 1
      ),
      axis.text.y = element_text(size = text_size, color = "black")
    ) +
    # --------------------------
    # 4. （可选）强制图例位置（避免默认靠上导致覆盖）
    # --------------------------
    theme(
      legend.position = "right",  # 固定图例在右侧（右侧空间更充足，减少顶部覆盖风险）
      # 若需在底部：legend.position = "bottom"，并配合调整legend.direction = "horizontal"
    ) +
    guides(
      fill = guide_legend(
        ncol = legend_cols,
        # （可选）在guides中再次确认宽高（优先级低于theme，但可双重保障）
        keywidth = unit(1.2, "lines"),
        keyheight = unit(1.2, "lines")
      )
    )    
}else if (plot_type == "alluvium_barplot") {
    # alluvium_barplot使用相同的宽度调整规则
    final_plot <- freq_data %>%
      ggplot(aes(x = .data[[x_var]], y = freq, fill = .data[[fill_var]])) +
      geom_stratum(
        aes(stratum = .data[[fill_var]]),
        color = NA,
        width = 0.8
      ) +
      geom_flow(
        aes(alluvium = .data[[fill_var]]),
        knot.pos = 0.23,
        width = 0.8,
        alpha = 0.5
      ) +
      geom_alluvium(
        aes(alluvium = .data[[fill_var]]),
        knot.pos = 0.23,
        color = "#f4e2de",
        width = 0.8,
        linewidth = 0.7,
        fill = NA,
        alpha = 1
      ) +
      scale_fill_manual(values = user_color_pal) +
      scale_x_discrete(name = "") +
      scale_y_continuous(
        expand = c(0, 0),
        name = "Proportion[%]"
      ) +
      theme_minimal() +
      theme(
         legend.spacing.y = unit(0.2, "lines"),  # 进一步减小条目间距（原1→0.2）
         legend.key.height = unit(1, "lines"),   # 缩小图例符号高度（原1.5→1）
         legend.key.size = unit(2, "lines"),     # 缩小符号整体尺寸（配合符号高度）
         legend.text = element_text(
         size = 23,
         margin = margin(t = 0.5, b = 0.5)  # 移除文本自身的额外边距（原2→0）
         ),
        legend.margin = margin(t = 2, b = 2, unit = "pt"),  # 缩小边缘空白（单位用pt更精细）
        legend.title = element_text(size = 25),
        text = element_text(family = "sans"),
        plot.background = element_rect(fill = 'white', color = 'white'),
        panel.grid = element_blank(),
        axis.title.y = element_text(size = text_size),
        axis.text.x = element_text(
          size = text_size, color = "black",
          angle = angle, vjust = 0.5, hjust = 1
        ),
        axis.text.y = element_text(size = text_size, color = "black"),
        plot.margin = margin(t = 14, r = 10, b = 10, l = 10)
      ) +
      guides(fill = guide_legend(ncol = legend_cols))
    
  }  else if (plot_type == "coord_flip_barplot") {
    adjusted_height <- max(base_width, base_width + num_categories * 0.5)
    final_plot <- freq_data %>%
      ggplot(aes(x = .data[[x_var]], y = freq, fill = .data[[fill_var]])) +
      geom_col(width = 0.9) +
      scale_fill_manual(values = user_color_pal) +
      scale_x_discrete(
        name = "",
        expand = c(-0.1, 0)
      ) +
      scale_y_continuous(
        expand = c(0, 0),
        name = "",
        labels = scales::percent_format(scale = 1)
      ) +
      theme_minimal() +
      theme(
         legend.spacing.y = unit(0.2, "lines"),  # 进一步减小条目间距（原1→0.2）
         legend.key.height = unit(1, "lines"),   # 缩小图例符号高度（原1.5→1）
         legend.key.size = unit(2, "lines"),     # 缩小符号整体尺寸（配合符号高度）
         legend.text = element_text(
         size = 23,
         margin = margin(t = 0.5, b = 0.5)  # 移除文本自身的额外边距（原2→0）
         ),
        legend.margin = margin(t = 2, b = 2, unit = "pt"),  # 缩小边缘空白（单位用pt更精细）
        legend.title = element_text(size = 25),
        text = element_text(family = "sans"),
        plot.background = element_rect(fill = 'white'),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = text_size, color = "black"),
        axis.text.y = element_text(
          size = text_size, color = "black",
          margin = margin(r = -1.5)
        ),
        plot.margin = margin(5, 5, 5, 5, "pt")
      ) +
      coord_flip() +
      guides(fill = guide_legend(ncol = legend_cols))
    
  } else if (plot_type == "solid_barplot") {
    # 为solid_barplot添加图例列数设置
    final_plot <- freq_data %>%
      ggplot(aes(x = .data[[x_var]], y = freq, fill = .data[[fill_var]])) +
      geom_stratum(
        aes(stratum = .data[[fill_var]]),
        color = NA,
        width = 0.9
      ) +
      coord_polar(theta = 'x') +
      scale_fill_manual(values = user_color_pal) +
      theme_classic() +
      theme(
        legend.spacing.y = unit(0.2, "lines"),  # 进一步减小条目间距（原1→0.2）
         legend.key.height = unit(1, "lines"),   # 缩小图例符号高度（原1.5→1）
         legend.key.size = unit(2, "lines"),     # 缩小符号整体尺寸（配合符号高度）
         legend.text = element_text(
         size = 20,
         margin = margin(t = 0.5, b = 0.5)  # 移除文本自身的额外边距（原2→0）
         ),
        legend.margin = margin(t = 2, b = 2, unit = "pt"),  # 缩小边缘空白（单位用pt更精细）
        legend.title = element_text(size = 20),
            text = element_text(family = "sans"),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = element_text(color = "black", vjust = 0.5, size = 20),
            axis.text.y = element_blank())+
            labs(x = "", y = "") +
      scale_y_continuous(expand = expansion(add = c(0.75, 0))) +
      guides(fill = guide_legend(ncol = legend_cols))  # 添加图例列数设置
  
  } else if (plot_type == "hollow_barplot") {
    # 为hollow_barplot添加图例列数设置
 final_plot <- ggplot(data = freq_data) +
      geom_col(aes(x = .data[[x_var]], y = freq, fill = .data[[fill_var]]),
               position = position_fill()) +
      scale_fill_manual(values = user_color_pal)  +
      coord_polar(theta = 'x') +
      theme_classic(base_size = 18) +
      theme(
         egend.spacing.y = unit(0.2, "lines"),  # 进一步减小条目间距（原1→0.2）
         legend.key.height = unit(1, "lines"),   # 缩小图例符号高度（原1.5→1）
         legend.key.size = unit(2, "lines"),     # 缩小符号整体尺寸（配合符号高度）
         legend.text = element_text(
         size = 20,
         margin = margin(t = 0.5, b = 0.5)  # 移除文本自身的额外边距（原2→0）
         ),
        legend.margin = margin(t = 2, b = 2, unit = "pt"),  # 缩小边缘空白（单位用pt更精细）
        legend.title = element_text(size = 20),
        text = element_text(family = "sans"),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = element_text(color = "black", vjust = 0.5, size = 20),
            axis.text.y = element_blank()) +
      labs(x = "", y = "") +
      scale_y_continuous(expand = expansion(add = c(0.75, 0))) +
      guides(fill = guide_legend(ncol = legend_cols))  # 添加图例列数设置
  
  } else {
    stop("Unsupported plot type: ", plot_type)
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  pdf_path <- file.path(output_dir, paste0(output_prefix, ".pdf"))
  png_path <- file.path(output_dir, paste0(output_prefix, ".png"))
  
  # 保存文件 (使用调整后的宽度和高度)
 ggplot2::ggsave(
    filename = pdf_path,
    plot = final_plot,
    device = "pdf",
    width = adjusted_width,
    height = adjusted_height
  )
  
  ggplot2::ggsave(
    filename = png_path,
    plot = final_plot,
    width = adjusted_width,
    height = adjusted_height,
    device = "png",
    dpi = 300,
    bg = "white"
  )
  print(adjusted_width)
  return(final_plot)
}
