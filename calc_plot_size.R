# 自定义函数：根据文本长度计算图片尺寸
calc_plot_size <- function(plot, axis = "x",unit = "in", dpi = 300, base_size = 8) { # 基础宽度/高度（单位：英寸）
  require(grid)
  # 获取绘图对象
  # grob <- ggplotGrob(plot)
  # 获取坐标轴标签文本
  if (axis == "x" | axis == "y") {
    #labels <- grob$grobs[[which(grob$layout$name == "axis-b")]]$children[[2]]$children[[1]]$label
    labels = unique(plot$data[,rlang::as_name(plot$mapping[[axis]])])

    # # 提取 x 轴变量名（可能是一个表达式，需转换为字符串）
    # x_var <- deparse(plot$mapping$x)
    # y_var <- deparse(plot$mapping$y)

    # # 示例输出
    # x_var  # 输出类似 "pred"
    # y_var  # 输出类似 "actual"

  # 计算最长文本的宽度（考虑旋转角度）
  } else {
    stop("axis 必须是 'x' 或 'y'")
  }
  # 计算最长文本的宽度（考虑旋转角度）
  rot = if (axis == "x") plot$theme$axis.text.x$angle else plot$theme$axis.text.y$angle
  fontsize = if (axis == "x") plot$theme$axis.text.x$size else plot$theme$axis.text.y$size
  max_width <- max(sapply(labels, function(label) {
    text_grob <- textGrob(
      label,
      gp = gpar(fontsize = if (is.null(fontsize)) 12 else fontsize),
      rot = if (is.null(rot)) 0 else rot
    )
    width <- convertWidth(grobWidth(text_grob), unit, valueOnly = TRUE)
    height <- convertHeight(grobHeight(text_grob), unit, valueOnly = TRUE)
    
    # 计算旋转后的投影宽度
    rad <- ifelse(is.null(rot), 0 ,rot) * pi / 180
    projected_width <- abs(width * cos(rad)) + abs(height * sin(rad))
    return(projected_width)
  }))
# 获取当前边距（默认单位：points）
current_margin <- plot$theme$plot.margin

# 将边距转换为目标单位（如 inches）
margin_in <- convertUnit(
  current_margin,
  unitTo = "inches",
  valueOnly = TRUE
)
margin_width = if (axis == "x") margin_in[2]+margin_in[4] else margin_in[1]+margin_in[3]
  # 基础图片尺寸 + 动态扩展尺寸
  dynamic_size <- base_size + max_width + margin_width  # 根据文本长度扩展
  return(dynamic_size)
}
