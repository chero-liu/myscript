######################
make_dimplot_data <- function(seurat_obj, reduct = "umap", groupby = "clusters", splitby = NULL, facet = NULL, get.color = TRUE, palette = "customecol2") {
    require("dplyr")
    # 公共数据准备部分
    data <- as.data.frame(seurat_obj@reductions[[reduct]]@cell.embeddings) %>% 
                tibble::rownames_to_column(var = "barcodes")
    if(get.color){
        object_meta <- get_color_use(seurat_obj, groupby, palette = palette)[["object_meta"]]
        anno_col = c(groupby,splitby,facet, paste0(groupby,"_col"))
    } else {
        object_meta <- seurat_obj@meta.data
        anno_col = c(groupby,splitby,facet)
    }
    anno <- object_meta[, anno_col, drop = FALSE] %>% 
                tibble::rownames_to_column(var = "barcodes")
    data <- merge(data, anno, by = "barcodes")
    if(class(seurat_obj@meta.data[,groupby]) == "factor"){
        data[,groupby] = factor(data[,groupby], levels = levels(seurat_obj@meta.data[,groupby]))
    }
    return(data)
}
dimplot_integrated <- function(
    data = NULL, 
    seurat_obj = NULL,
    mode = c("level", "contour", "polygon", "none"),  # 新增模式选择参数(对应：等高线、轮廓线、背景阴影,无任何修饰)
    reduct = "umap", 
    pt.size = 0.5,
    groupby = "clusters", 
    splitby = NULL,
    ncol = 2,
    color_use = NULL, 
    palette = "customecol2",
    special_cells = NULL, # 指定突出显示的细胞类型（仅限轮廓线、背景阴影图形）
    combine = FALSE, # 是否合并special_cells出轮廓图
    total_cells = FALSE, # 是否显示总细胞数量
    cell_numbers = FALSE, # 是否显示细胞数量
    cell_tag = "cells", # 是否显示细胞数量
    label = TRUE,  # 是否显示标签
    label_color = FALSE, # 标签颜色是否跟随细胞颜色
    panel_border = FALSE, # 是否显示边框
    minDensity = 0.1, # 最小密度阈值，数值越大，密度范围越小,轮廓更贴合
    smoothSigma = 0.02, # 数值越小，轮廓越清晰(越平滑)
    masktype = "partition" # polygon模式: partition 互不重叠; independent 允许重叠
) {
    require("ggplot2")
    if( ! is.null(seurat_obj)){
        data <- make_dimplot_data(seurat_obj, reduct = reduct, groupby = groupby, splitby = splitby, facet = NULL, get.color = TRUE, palette = palette)
    }else if( is.null(data)){
        stop("Please input seurat_obj or plot_data")
    }
    # 公共数据准备部分
    reduct1 <- colnames(data)[2]
    reduct2 <- colnames(data)[3]

    # 公共颜色设置
    if (is.null(color_use)) {
        color_use <- get_colors(data, groupby, palette = palette)[["new_celltype_pal"]]
    }
    counts <- table(data[, groupby])
     if (cell_numbers) {
        cell_count_labels <- paste(names(counts), counts, sep = "-") %>% paste(., cell_tag)
    } else {
        if( is.null(names(color_use))){
            cell_count_labels <- names(counts)
        }else{
            cell_count_labels <-  names(color_use) 
        }
        if( groupby == "clusters" ){
            cell_count_labels = paste0("cluster ",cell_count_labels)
        }        
    }
    
    # # 公共坐标范围计算 (长宽比例尺一致版本)
    # max_diff = ceiling(max(diff(range(data[, reduct1])), diff(range(data[, reduct2])))) * 1.1
    # diff_x = (max_diff - diff(range(data[, reduct1])))/2 
    # diff_y = (max_diff - diff(range(data[, reduct2])))/2 
    
    # 输入：原始数据坐标范围
    x_range <- range(data[, reduct1])  # c(min_x, max_x)
    y_range <- range(data[, reduct2])  # c(min_y, max_y)
    # 将坐标轴的扩大1.3倍
    expand_factor <- 0.15  # (1.1 - 1) / 2
    x_limits <- c(x_range[1] - diff(x_range) * expand_factor, 
                x_range[2] + diff(x_range) * expand_factor)
    y_limits <- c(y_range[1] - diff(y_range) * expand_factor, 
                y_range[2] + diff(y_range) * expand_factor)

    # 确定绘图长宽比例
    # ratio = diff(range(y_limits))/diff(range(x_limits))

    # 使用switch控制不同模式
    gg_tsne <- switch(
        match.arg(mode),
        "level" = {
            # level等高线模式专属逻辑
            bwidth <- 0.05 * c(sum(abs(range(data[, reduct1]))), 
                              sum(abs(range(data[, reduct2]))))
            ggplot(data, aes_string(reduct1, reduct2, colour = groupby)) +
                stat_density2d(
                    aes_string(x = reduct1, y = reduct2, colour = groupby),
                    alpha = 0.5, h = bwidth, contour_var = "ndensity",
                    show.legend = FALSE
                ) +
                geom_point(size = pt.size, stroke = 0.1)
        },
        "contour" = {
            # contour轮廓线模式专属逻辑
            if( combine & !is.null(special_cells)){
                special_list = unlist(strsplit(special_cells,":", perl =T)) %>% strsplit(split = ",",perl =T)
                data_cluster = NULL
                for(i in 1:length(special_list)){
                    special = special_list[[i]]
                    data_sub = data[which(data[, groupby] %in% special ),]
                    data_sub[,groupby] <- paste0("combine",i)
                    data_cluster = rbind(data_cluster,data_sub)
                }
                maskTable <- generateMask(
                        dims = as.data.frame(data[,c(reduct1,reduct2)])[rownames(data_cluster),],
                        clusters = data_cluster[, groupby],
                        gridSize = 500, # 轮廓的分辨率
                        minDensity = minDensity, # 最小密度阈值，数值越大，密度范围越小,轮廓更贴合
                        smoothSigma = smoothSigma, # 数值越小，轮廓越清晰(越平滑)
                        type="partition" # partition 互不重叠; independent 允许重叠
                    ) %>% as.data.frame() %>% dplyr::rename(!!groupby := cluster) 
            } else {
                maskTable <- generateMask(
                        dims = as.data.frame(data[,c(reduct1,reduct2)]),
                        clusters = data[, groupby],
                        gridSize = 500, # 轮廓的分辨率
                        minDensity = minDensity, 
                        smoothSigma = smoothSigma,
                        type = masktype
                    ) %>% as.data.frame() %>% dplyr::rename(!!groupby := cluster)
                if(!is.null(special_cells)){
                    special_cells = unlist(strsplit(special_cells,":", perl =T)) %>% strsplit(split = ",",perl =T) %>% unlist()
                    maskTable <- maskTable[which(maskTable[,groupby] %in% special_cells ),]
                }
            }
            ggplot(data, aes_string(reduct1, reduct2, colour = groupby)) +
                geom_point(size = pt.size, stroke = 0.1) +
                geom_path(
                    data = maskTable,
                    aes_string(reduct1, reduct2, colour = NULL, group = "part"),
                    linewidth = 0.5, linetype = "longdash",
                    show.legend = FALSE, color = "black"
                )
        },
        "polygon" = {
            # polygon阴影背景模式专属逻辑
            maskTable <- generateMask(
                dims = as.data.frame(data[,c(reduct1,reduct2)]),
                clusters = data[, groupby],
                minDensity = minDensity,
                smoothSigma = smoothSigma
            ) %>% as.data.frame() %>% dplyr::rename(!!groupby := cluster)
            if(!is.null(special_cells)){
                special_cells = unlist(strsplit(special_cells,":", perl =T)) %>% strsplit(split = ",",perl =T) %>% unlist()
                maskTable <- maskTable[which(maskTable[,groupby] %in% special_cells ),]
            }
            part = unique(maskTable[,c("part",groupby)])
            color_use_fill = color_use[match(part[,groupby],names(color_use))]
            names(color_use_fill) = part$part
            ggplot(data, aes_string(reduct1, reduct2, colour = groupby)) +
                geom_polygon(
                    data = maskTable,
                    aes_string(reduct1, reduct2, colour = NULL, fill = "part"),
                    alpha = 0.3,
                    show.legend = FALSE
                ) +
                geom_point(size = pt.size, stroke = 0.1) +
                scale_fill_manual(values = color_use_fill)
        },
        "none" = {
            # 无任何背景修饰模式专属逻辑
            ggplot(data, aes_string(reduct1, reduct2, colour = groupby)) +
                geom_point(size = pt.size, stroke = 0.1) 
        }
    )
    # min_x = ifelse( diff_x/3 + max_diff/30 > diff_x, diff_x/3 + max_diff/25, diff_x)
    # min_y = ifelse( diff_y/3 + max_diff/30 > diff_y, diff_y/3 + max_diff/25, diff_y)
    # 公共图形设置部分
    gg_tsne <- gg_tsne +
        theme_test() +
        scale_x_continuous(
            limits = x_limits
        ) +
        scale_y_continuous(
            limits = y_limits
        ) +
        scale_colour_manual(
            values = color_use,
            breaks = names(color_use),
            labels = cell_count_labels
        ) +
        guides(fill = FALSE) +
        Seurat::NoAxes() +
        # 公共坐标轴箭头设置
        geom_segment(
            aes(x = (x_range[1] + x_limits[1]) / 2, 
                y = (y_range[1] + y_limits[1]) / 2,
                xend = (x_range[1] + x_limits[1]) / 2 + diff(x_limits)/10, 
                yend = (y_range[1] + y_limits[1]) / 2 ),
            colour = "black", size = 0.6,
            arrow = arrow(length = unit(0.3, "cm"))
        ) +
        geom_segment(
            aes(x = (x_range[1] + x_limits[1]) / 2, 
                y = (y_range[1] + y_limits[1]) / 2,
                xend = (x_range[1] + x_limits[1]) / 2, 
                yend = (y_range[1] + y_limits[1]) / 2 + diff(y_limits)/10),
            colour = "black", size = 0.6,
            arrow = arrow(length = unit(0.3, "cm"))
        ) +
        # 公共文本标注
        annotate(
            "text", 
            x = (x_range[1] + x_limits[1]) / 2 + diff(x_limits)/20,
            y = y_range[1] - (y_range[1] - y_limits[1]) * 0.75,
            label = gsub("_"," ",reduct1), color = "black", size = 5
        ) +
        annotate(
            "text", 
            x = x_range[1] - (x_range[1] - x_limits[1]) * 0.75,
            y = (y_range[1] + y_limits[1]) / 2 + diff(y_limits)/20,
            label = gsub("_"," ",reduct2), color = "black", size = 5, angle = 90
        ) +
        theme(
            legend.title = element_blank(),
            legend.key = element_rect(fill = "white"),
            legend.text = element_text(size = 20),
            aspect.ratio = 1/1
        ) +
        guides(color = guide_legend(override.aes = list(size = 5)))
    
    # 公共边框设置
    if (total_cells) {
        gg_tsne <- gg_tsne + 
                annotate(
                "text",
                x = x_range[1] - (x_range[1] - x_limits[1]) * 0.75,
                y = (y_range[2] + y_limits[2]) / 2,
                label = paste0("Total ", format(dim(data)[1], big.mark = ","), " ", cell_tag),
                color = "black", size = 8, hjust = 0, vjust = 0
            ) 
    }
    # 公共边框设置
    if (panel_border) {
        gg_tsne <- gg_tsne + 
            theme(panel.border = element_rect(size = 1.5, color = "black"))
    } else {
        gg_tsne <- gg_tsne + 
            theme(panel.border = element_blank())
    }
    
    # 公共标签设置
    if (label) {
        cell_type_med <- data %>%
            group_by(!!sym(groupby)) %>%
            summarise(
                reduct_1 = median(!!sym(reduct1)),
                reduct_2 = median(!!sym(reduct2))
            )
        
        text_color <- if (label_color) {
            adjustcolor(
                color_use,
                alpha.f = 1, 
                blue.f = 0.6, 
                green.f = 0.6, 
                red.f = 0.6
            )
        } else {
            "black"
        }
        
        gg_tsne <- gg_tsne +
            ggrepel::geom_text_repel(
                aes_string("reduct_1", "reduct_2", label = groupby),
                data = cell_type_med,
                point.padding = unit(0, "lines"),
                colour = text_color,
                size = 7, seed = 2025
            )
    }
    if(!is.null(splitby)){
        gg_tsne <- gg_tsne + facet_wrap(facets = vars(!!sym(x = splitby)), ncol = ncol) + theme(strip.background = element_blank(), strip.text = element_text(size = 20))
    }
    
    return(gg_tsne)
}
###################################################################
#### 计算分群轮廓 改编来源mascarade包   #############################
#### https://github.com/alserglab/mascarade  ######################
expandedRange2d <- function(x, y, fraction=0.05, fixAspectRatio=TRUE) {
    xRange <- range(x)
    xWidth <- (xRange[2] - xRange[1]) * (1 + fraction)

    yRange <- range(y)
    yWidth <- (yRange[2] - yRange[1]) * (1 + fraction)

    if (fixAspectRatio) {
        xWidth <- yWidth <- max(xWidth, yWidth)
    }

    xCenter <- mean(xRange)
    yCenter <- mean(yRange)

    return(
        c(xCenter - xWidth/2, xCenter + xWidth/2,
          yCenter - yWidth/2, yCenter + yWidth/2))
}

#' @importFrom spatstat.geom tiles tess connected as.polygonal
#' @importFrom  data.table rbindlist
borderTableFromMask <- function(curMask) {
    require(spatstat.geom)
    require(data.table)
    parts <- spatstat.geom::tiles(spatstat.geom::tess(image=spatstat.geom::connected(curMask > 0)))

    curBorderTable <- list()

    for (partIdx in seq_along(parts)) {
        part <- parts[[partIdx]]

        partBoundary <- spatstat.geom::as.polygonal(part)
        lines <- partBoundary$bdry

        curBorderTable <- c(curBorderTable, lapply(seq_along(lines), function(lineIdx) {
            curLine <- lines[[lineIdx]]
            xs <- curLine$x
            ys <- curLine$y

            # make lines closed
            xs <- c(xs, xs[1])
            ys <- c(ys, ys[1])

            # remove steps
            xs <- (head(xs, -1) + tail(xs, -1)) / 2
            ys <- (head(ys, -1) + tail(ys, -1)) / 2

            # make lines closed again
            xs <- c(xs, xs[1])
            ys <- c(ys, ys[1])

            res <- data.table::data.table(x=xs, y=ys)
            res[, part := as.character(partIdx)] # 避免因为partIdx为数值而导致数据格式问题
            res[, linegroup := as.character(lineIdx)]
            res[]
        }))
    }
    data.table::rbindlist(curBorderTable)
}

#' Generate mask for clusters on 2D dimensional reduction plots
#'
#' Internally the function rasterizes and smoothes the density plots.
#' @param dims matrix of point coordinates.
#'      Rows are points, columns are dimensions. Only the first two columns are used.
#' @param clusters vector of cluster annotations.
#'      Should be the same length as the number of rows in `dims`.
#' @param gridSize width and height of the raster used internally
#' @param minDensity minimal required density for the grid cells to be included in the mask.
#'      Decreasing this parameter will expand masks.
#' @param smoothSigma sigma used in Gaussian smoothing represented as a fraction of plot width.
#'      Increasing this parameter can help dealing with sparse regions.
#' @param type controls the behavior of the method.
#'      When set to "partition" (default) generated masks are mutually exclusive.
#'      When set to "independent" masks can overlap.
#' @returns data.table with points representing the mask borders.
#'      Each individual border line corresponds to a single level of `group` column.
#'      Cluster assignment is in `cluster` column.
#' @importFrom data.table rbindlist data.table setnames
#' @importFrom utils head tail
#' @export
#' @examples
#' data("exampleMascarade")
#' res <- generateMask(dims=exampleMascarade$dims,
#'                     clusters=exampleMascarade$clusters)
#' \dontrun{
#' data <- data.table(exampleMascarade$dims,
#'                    cluster=exampleMascarade$clusters,
#'                    exampleMascarade$features)
#' ggplot(data, aes(x=UMAP_1, y=UMAP_2)) +
#'     geom_point(aes(color=cluster)) +
#'     geom_path(data=maskTable, aes(group=group)) +
#'     coord_fixed() +
#'     theme_classic()
#' }
generateMask <- function(dims, clusters,
                         gridSize=200,
                         minDensity=0.1,
                         smoothSigma=0.025,
                         type=c("partition", "independent")) {
    type <- match.arg(type)
    require(spatstat.geom)
    require(data.table)

    clusterLevels <- unique(clusters)

    gridRange <- expandedRange2d(dims[, 1], dims[ ,2])
    windowWidth <- gridRange[2] - gridRange[1]

    smoothSigma <- smoothSigma * windowWidth

    window <- spatstat.geom::as.mask(spatstat.geom::owin(xrange = gridRange[1:2],
                                                         yrange = gridRange[3:4]),
                                     dimyx=gridSize)

    dims <- dims[, 1:2]

    if (is.null(colnames(dims))) {
        colnames(dims) <- c("x", "y")
    }

    points <- spatstat.geom::ppp(dims[, 1], dims[, 2], window=window)


    allDensities <- lapply(clusterLevels, function(cluster) {
        res <- spatstat.geom::pixellate(points[clusters == cluster], dimyx=gridSize)
    })

    allDensitiesSmoothed <- lapply(allDensities, spatstat.explore::blur, sigma = smoothSigma)

    densityThresholds <- pmin(vapply(allDensitiesSmoothed, max, numeric(1)) / 2,
                              minDensity)

    if (type == "partition") {
        # backgroundDensity <- spatstat.geom::as.im(window) * minDensity

        allDensitiesMax <- spatstat.geom::im.apply(allDensitiesSmoothed, which.max)
        allDensitiesSmoothed <- lapply(seq_along(clusterLevels), function(i) {
            allDensitiesSmoothed[[i]] * (allDensitiesMax == i)
        })
    }

    borderTable <- data.table::rbindlist(lapply(seq_along(clusterLevels), function(i) {
        curMask <- (allDensitiesSmoothed[[i]] >= densityThresholds[i]) * allDensitiesSmoothed[[i]]
        if (sum(curMask) == 0) {
            warning(sprintf("Mask is empty for cluster %s", clusterLevels[i]))
            return(NULL)
        }
        # 使用 connected 函数确保每个离散部分被视为独立区域
        curMaskConnected <- spatstat.geom::connected(curMask > 0)
        curTableList <- lapply(levels(curMaskConnected ), function(j) {
            curPartMask <- curMask * (curMaskConnected == j)
            if (sum(curPartMask) == 0) return(NULL)
            curTable <- borderTableFromMask(curPartMask)
            curTable[, cluster := clusterLevels[i]]
            curTable[, part := paste0(cluster, "#", j)]
            curTable[, linegroup := paste0(part, "#", linegroup)]
            curTable[]
        })
        curTableList <- curTableList[!sapply(curTableList, is.null)]
        data.table::rbindlist(curTableList)
    }))

    data.table::setnames(borderTable, c("x", "y"), colnames(dims))

    return(borderTable)
}

