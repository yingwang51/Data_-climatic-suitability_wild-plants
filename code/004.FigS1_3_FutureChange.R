# ============================================================================
# 004.FigS1_3_FutureChange.R   (独立出图脚本)
# FigS1-3 — 未来适生区变化图 (每 SSP 一张大图; 逐 GCM 展示, 不做 ensemble)
#
# 布局 (每张图, 列=物种 6, 行=GCM 3):
#   A: 3行(GCM) × 6列(物种) 小地图, 三色展示 Lost / Stable / Gained
#      —— 每张图为单个 GCM 的当前→未来变化 (非 ensemble)
#   B: 6列柱形图, 每列一个物种内含 3 个 GCM 的 diverging bar
#      (丧失向下 / 保留+扩张向上)
#
# 优点: 每个 GCM 地图的色块面积 == 它对应的柱子, 地图与柱子口径完全一致;
#       且与正文 Fig3/Fig4 (GCM 为重复) 的原始数据同源 (area_raw_data.csv)。
#
# 物种按当前面积从小到大排序 (与 Fig6 一致)
# ============================================================================

rm(list = ls()); gc()
suppressMessages({
  library(terra); library(sf); library(tidyverse)
  library(patchwork)
  library(scales)
})

# --- 路径 ---
work_dir      <- "K:/周建maxent"
output_dir    <- file.path(work_dir, "MaxEnt_Fixed_Results")
res_dir       <- file.path(output_dir, "gcm_replicate_analysis")
fig_dir       <- file.path(res_dir, "figures")
tbl_dir       <- file.path(res_dir, "tables")
boundary_file <- file.path(work_dir, "1矢量边界/大兴安岭范围准确.shp")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# --- 物种 (与 004b.Fig1 完全一致：面积升序 + 拉丁名缩写) ---
species_cn_map <- c("偃松" = "yangsong", "杜香" = "duxiang", "越桔" = "yuejv",
                     "笃斯越桔" = "dusiyuejv", "小黄花菜" = "xiaohuanghuacai",
                     "短瓣金莲花" = "duanbaijinlianhua")
species_latin <- c(yangsong = "Pinus pumila", duxiang = "Rhododendron tomentosum",
                   yuejv = "Vaccinium vitis-idaea", dusiyuejv = "Vaccinium uliginosum",
                   xiaohuanghuacai = "Hemerocallis minor",
                   duanbaijinlianhua = "Trollius ledebourii")
species_latin_abbr <- sub("^([A-Z])[a-z]+ ", "\\1. ", species_latin)

# 面积升序 (与 004b 相同)
area_raw   <- read.csv(file.path(tbl_dir, "area_raw_data.csv"), fileEncoding = "UTF-8")
cur_area   <- area_raw %>% distinct(species_code, current_area_km2) %>% arrange(current_area_km2)
species_order <- cur_area$species_code   # c("yangsong", "duxiang", "yuejv", ...)

# 统一 y 轴范围 (三张 SSP 图同尺度, 便于比较); 基于各 GCM 原始值, 留余量取整到万
y_hi <- ceiling(max(area_raw$stable_km2 + area_raw$gained_km2) * 1.05 / 10000) * 10000
y_lo <- -ceiling(max(area_raw$lost_km2) * 1.05 / 10000) * 10000
y_brks <- seq(y_lo, y_hi, by = 20000)

# --- GCM 与 SSP ---
gcms      <- c("BCC-CSM2-MR", "MIROC6", "IPSL-CM6A-LR")
gcm_abbr  <- c("BCC-CSM2-MR" = "BCC", "MIROC6" = "MIROC6", "IPSL-CM6A-LR" = "IPSL")
ssps      <- c("ssp126", "ssp245", "ssp585")
ssp_labels <- c(ssp126 = "SSP1-2.6", ssp245 = "SSP2-4.5", ssp585 = "SSP5-8.5")

# --- 颜色 ---
COL_LOST   <- "#D32F2F"
COL_STABLE <- "#2E7D32"
COL_GAINED <- "#1976D2"

# --- 边界 ---
boundary_wgs <- st_transform(st_read(boundary_file, quiet = TRUE), 4326)

# ============================================================================
# 辅助函数
# ============================================================================
read_binary <- function(sp_code, fname) {
  mask(crop(rast(file.path(output_dir, sp_code, fname)), boundary_wgs), boundary_wgs)
}

gcm_ensemble <- function(sp_code, ssp) {
  rlist <- lapply(gcms, function(gcm) {
    read_binary(sp_code, sprintf("binary_%s_%s_tss.tif", gcm, ssp))
  })
  stk     <- rast(rlist)
  n_valid <- sum(!is.na(stk))
  gcm_sum <- sum(stk, na.rm = TRUE)
  ens     <- gcm_sum >= ceiling(n_valid / 2)
  ifel(n_valid < 2, NA, ens)
}

compute_change <- function(cur_rast, fut_rast) {
  cur_v <- values(cur_rast)
  fut_v <- values(fut_rast)
  change <- rast(cur_rast)
  values(change) <- 0L
  change[cur_v == 1 & fut_v == 0] <- 1L
  change[cur_v == 1 & fut_v == 1] <- 2L
  change[cur_v == 0 & fut_v == 1] <- 3L
  change
}

extract_areas <- function(change_rast) {
  cell_area_km2 <- cellSize(change_rast, unit = "km")
  v  <- values(change_rast)
  ca <- values(cell_area_km2)
  c(lost   = sum(ca[v == 1], na.rm = TRUE),
    stable = sum(ca[v == 2], na.rm = TRUE),
    gained = sum(ca[v == 3], na.rm = TRUE))
}

# 三色变化地图 (data.frame 方式, 最稳定; 不设独立 legend, 由外部统一收集)
mk_change_map <- function(change_rast, title) {
  df <- as.data.frame(change_rast, xy = TRUE, na.rm = TRUE)
  names(df)[3] <- "value"
  df <- df[df$value > 0, ]
  df$category <- factor(df$value, levels = 1:3,
                        labels = c("Lost", "Stable", "Gained"))

  ggplot() +
    geom_raster(data = df, aes(x, y, fill = category)) +
    scale_fill_manual(
      values = c("Lost" = COL_LOST, "Stable" = COL_STABLE, "Gained" = COL_GAINED),
      name = NULL, drop = FALSE, guide = "none"
    ) +
    geom_sf(data = boundary_wgs, fill = NA, color = "grey30", linewidth = 0.15) +
    coord_sf(expand = FALSE) +
    labs(title = title) +
    theme_void(base_size = 10) +
    theme(plot.title = element_text(face = "italic", hjust = 0.5, size = 10.5))
}

# ============================================================================
# 逐 SSP 出大图
# ============================================================================

fig_idx <- 0
for (ssp in ssps) {
  fig_idx <- fig_idx + 1
  cat(sprintf("\n=== %s ===\n", ssp_labels[ssp]))

  # --- 逐物种构建一列: 3 个 GCM 地图(上) + 3-GCM diverging 柱子(下) ---
  combined_panels <- lapply(seq_along(species_order), function(i) {
    code     <- species_order[i]
    sp_latin <- species_latin_abbr[[code]]
    is_first <- (i == 1)

    cur_rast <- read_binary(code, "binary_current_tss.tif")

    # 3 个 GCM 各一张变化地图 (首行附物种名标题)
    maps <- lapply(seq_along(gcms), function(g) {
      cat(sprintf("  %s/%s...", code, gcms[g]))
      gcm_rast <- read_binary(code, sprintf("binary_%s_%s_tss.tif", gcms[g], ssp))
      chg <- compute_change(cur_rast, gcm_rast)
      mk_change_map(chg, if (g == 1) sp_latin else "")
    })

    # 3-GCM diverging 柱子: x = 1,2,3 对应三个 GCM
    ar <- area_raw[area_raw$species_code == code & area_raw$scenario == ssp, ]
    ar <- ar[match(gcms, ar$gcm), ]
    bw <- 0.38
    pos_df <- do.call(rbind, lapply(seq_len(3), function(g) data.frame(
      x = g,
      ymin = c(0, ar$stable_km2[g]),
      ymax = c(ar$stable_km2[g], ar$stable_km2[g] + ar$gained_km2[g]),
      category = factor(c("Stable", "Gained"), levels = c("Stable", "Gained"))
    )))
    neg_df <- data.frame(x = 1:3, ymin = -ar$lost_km2, ymax = 0,
                         category = factor("Lost"))
    y_lab <- if (is_first) expression("Area (km"^2*")") else NULL

    bar_p <- ggplot() +
      geom_rect(data = pos_df,
                aes(xmin = x - bw, xmax = x + bw, ymin = ymin, ymax = ymax, fill = category)) +
      geom_rect(data = neg_df,
                aes(xmin = x - bw, xmax = x + bw, ymin = ymin, ymax = ymax, fill = category)) +
      geom_hline(yintercept = 0, linewidth = 0.3) +
      scale_fill_manual(
        values = c("Lost" = COL_LOST, "Stable" = COL_STABLE, "Gained" = COL_GAINED),
        guide = "none"
      ) +
      scale_y_continuous(limits = c(y_lo, y_hi), expand = c(0, 0),
                         breaks = y_brks,
                         labels = function(x) label_comma()(abs(x))) +
      scale_x_continuous(breaks = 1:3, labels = gcm_abbr[gcms], limits = c(0.4, 3.6)) +
      labs(x = NULL, y = y_lab) +
      theme_classic(base_size = 11) +
      theme(axis.text.x     = element_text(size = 7, angle = 45, hjust = 1),
            axis.text.y     = if (is_first) element_text(size = 9) else element_blank(),
            axis.ticks.y    = if (is_first) element_line() else element_blank(),
            axis.line.y     = if (is_first) element_line() else element_blank(),
            axis.ticks.length = unit(0.12, "cm"))

    # 一列纵向拼接: 3 地图 + 柱子
    wrap_plots(c(maps, list(bar_p)), ncol = 1, heights = c(2, 2, 2, 2.4))
  })
  cat("\n")

  # --- 6 列并排; 顶部加 spacer 空行, 留白给图例 (避免压住物种名标题) ---
  fig <- wrap_plots(combined_panels, nrow = 1)
  fig <- wrap_plots(plot_spacer(), fig, ncol = 1, heights = c(0.045, 0.955))

  out_file <- file.path(fig_dir, sprintf("FigS%d_FutureChange_%s.png", fig_idx, ssp))
  png(out_file, width = 14, height = 13, units = "in", res = 300)
  grid::grid.newpage()
  grid::grid.draw(fig)

  # A / B 面板标签 (主体占下方 95.5%, 故乘以缩放)
  grid::grid.text("A", x = unit(0.015, "npc"), y = unit(0.925, "npc"),
                  just = c("left", "top"), gp = grid::gpar(fontsize = 20))
  grid::grid.text("B", x = unit(0.015, "npc"), y = unit(0.225, "npc"),
                  just = c("left", "top"), gp = grid::gpar(fontsize = 20))

  # 左侧 GCM 行标签 (竖排, 对齐 3 行地图中心)
  gcm_y <- c(0.775, 0.555, 0.335)
  for (g in seq_along(gcms)) {
    grid::grid.text(gcm_abbr[gcms][g], x = unit(0.012, "npc"), y = unit(gcm_y[g], "npc"),
                    rot = 90, gp = grid::gpar(fontsize = 13))
  }

  # 三色图例 (Lost / Stable / Gained), 顶部空行居中横排
  leg <- data.frame(
    lab = c("Lost", "Stable", "Gained"),
    col = c(COL_LOST, COL_STABLE, COL_GAINED),
    x   = c(0.40, 0.49, 0.60)
  )
  for (k in seq_len(nrow(leg))) {
    grid::grid.rect(x = unit(leg$x[k], "npc"), y = unit(0.978, "npc"),
                    width = unit(0.012, "npc"), height = unit(0.014, "npc"),
                    just = c("left", "center"),
                    gp = grid::gpar(fill = leg$col[k], col = NA))
    grid::grid.text(leg$lab[k], x = unit(leg$x[k] + 0.015, "npc"),
                    y = unit(0.978, "npc"), just = c("left", "center"),
                    gp = grid::gpar(fontsize = 13))
  }
  dev.off()
  cat(sprintf("  -> %s\n", basename(out_file)))
}

cat("\n=== All done ===\n")
