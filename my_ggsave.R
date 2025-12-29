# 解决ggsave调用Cairo包的字体的问题
ggsave <- function(filename, plot = last_plot(), ...) {
  dots <- list(...)
#   if (!"filename" %in% names(dots)) stop("必须提供 filename 参数")
#   filename <- dots$filename
  ext <- tools::file_ext(filename)
  base_name <- tools::file_path_sans_ext(filename)
  if(ext == ""){
    ext = c("pdf","png")
  }
  # 自定义设备函数
    # 用户是否传入了 width/height？
    width <- ifelse("width" %in% names(dots), dots$width, 7)
    height <- ifelse("height" %in% names(dots), dots$height, 7)
  # 设置默认设备
  for (e in ext){
      device <- switch(
        tolower(e),
        pdf = function(...) Cairo::CairoPDF(file = paste0(base_name,".",e),width = width, height = height, bg = "white"),
        png = Cairo::CairoPNG,
        svg = function(...) Cairo::CairoSVG(file = paste0(base_name,".",e),width = width, height = height),
        jpeg = Cairo::CairoJPEG,
        jpg = Cairo::CairoJPEG,
        tiff = Cairo::CairoTIFF,
        bmp = Cairo::CairoBMP,
        stop(paste0("不支持的文件格式: ", e))
      )
      if(tolower(e) == "png"){
        dpi = ifelse("dpi" %in% names(dots), dots$dpi, 300)

        dots$width <- round(width * dpi)
        dots$height <- round(height * dpi)
        dots$dpi <- dpi
      }
      dots$limitsize <- ifelse("limitsize" %in% names(dots), dots$limitsize, FALSE)
      # 使用 do.call 调用 ggsave，避免参数冲突
      args <- c(list(paste0(base_name,".",e), plot = plot), dots, list(device = device))

      do.call(ggplot2::ggsave, args)
      message("Saved: ", paste0(base_name,".",e), " ; height:",dots$height, " | width:",dots$width)
  }
}
