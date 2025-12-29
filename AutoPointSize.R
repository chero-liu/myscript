#' 智能点大小计算函数
#' 
#' @param data 数据对象
#' @param raster 是否使用栅格化
#' @param min_size 最小点大小
#' @param max_size 最大点大小
#' @param min_cells 最少细胞数（仅影响缩放比例）
#' @param max_cells 最多细胞数（细胞数大于该数字的点大小均为min_size）
#' @param method 计算方法 ("log", "linear", "sqrt")
#' @return 点大小数值
AutoPointSize <- function(data, raster = FALSE, min_size = 0.1, max_size = 2, method = "log", min_cells = 500, max_cells = 500000) {
  n_points <- get_n_points(data)
  # 如果使用栅格化，可以使用较大点
  if (raster) {
    min_size <- min_size * 2
    max_size <- max_size * 1.5
  }
  # 根据不同方法计算点大小
  if (n_points <= 0) {
    return(max_size)
  }
  # 防止log(0)
  n_points <- max(n_points, 1)
  
  # 根据不同方法计算点大小
  if (method == "log") {
    # 对数缩放
    log_min <- log10(min_cells)    # log10(500) ≈ 2.7
    log_max <- log10(max_cells)    # log10(500000) ≈ 5.7
    
    log_n <- log10(n_points)
    # 限制在合理范围内
    # clamped_log <- pmin(pmax(log_n, log_min), log_max)
    # 归一化到0-1之间（注意是反转的：细胞数越多，值越小）
    normalized <- (log_max - log_n) / (log_max - log_min)
    pt_size <- min_size + normalized * (max_size - min_size)
    
  } else if (method == "sqrt") {
    # 平方根缩放
    sqrt_min <- sqrt(min_cells)    # sqrt(500) ≈ 22.4
    sqrt_max <- sqrt(max_cells)    # sqrt(500000) ≈ 707.1
    
    sqrt_n <- sqrt(n_points)
    # 限制在合理范围内
    # clamped_sqrt <- pmin(pmax(sqrt_n, sqrt_min), sqrt_max)
    # 归一化到0-1之间
    normalized <- (sqrt_max - sqrt_n) / (sqrt_max - sqrt_min)
    pt_size <- min_size + normalized * (max_size - min_size)
    
  } else {
    # 线性缩放
    # 限制在合理范围内
    # clamped_n <- pmin(pmax(n_points, min_cells), max_cells)
    # 归一化到0-1之间
    normalized <- (max_cells - n_points) / (max_cells - min_cells)
    pt_size <- min_size + normalized * (max_size - min_size)
  }
  
  # 确保在合理范围内
  # pt_size <- pmax(min_size, pmin(max_size, pt_size))
  # 不设上限
  pt_size <- pmax(min_size, pt_size)
  
  return(round(pt_size, 2))
}
# 获取数据点数量
get_n_points <- function(x) {
        if (is.data.frame(x) || is.matrix(x)) {
          return(nrow(x))
        } else if (inherits(x, "Seurat")) {
          # 支持Seurat对象
          return(ncol(x))
        } else if (inherits(x, "SingleCellExperiment")) {
          # 支持SingleCellExperiment对象
          return(ncol(x))
        } else if (length(x) ==1 & is.numeric(x)) {
          # 支持仅传入数字，根据数值设置大小
          return(x)
        } else if (is.vector(x)) {
          return(length(x))
        } else {
          # 尝试获取长度
          return(length(x))
        }
    }

