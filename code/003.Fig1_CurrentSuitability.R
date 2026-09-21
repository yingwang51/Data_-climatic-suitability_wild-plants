# ============================================================================
# 003_CurrentSuitability.R   (独立出图脚本)
# Fig.1 — A: 6物种当前适生二值图(单排)  +  B: 当前适生面积柱形图
#
# 独立运行, 不重跑面积提取:
#   - 地图直接读各物种 binary_current_tss.tif
#   - 面积读 002 已产出的 area_raw_data.csv
#   便于快速调整 Fig1 而无需重跑 004 全流程。
# 全图统一颜色 = 适生绿 (GREEN, 一处可改); 物种按当前面积从小到大排序。
# ============================================================================

rm(list = ls()); gc()

library(terra)
library(sf)
library(tidyverse)
library(tidyterra)   # geom_spatraster
library(patchwork)   # A/B 拼图

# --- 路径 ---
work_dir      <- "K:/周建maxent"
output_dir    <- file.path(work_dir, "MaxEnt_Fixed_Results")
res_dir       <- file.path(output_dir, "gcm_replicate_analysis")
fig_dir       <- file.path(res_dir, "figures")
tbl_dir       <- file.path(res_dir, "tables")
boundary_file <- file.path(work_dir, "1矢量边界/大兴安岭范围准确.shp")

# --- 物种与拉丁名 (与 004 一致) ---
species_map <- c("偃松" = "yangsong", "杜香" = "duxiang", "越桔" = "yuejv",
                 "笃斯越桔" = "dusiyuejv", "小黄花菜" = "xiaohuanghuacai",
                 "短瓣金莲花" = "duanbaijinlianhua")
species_latin <- c(yangsong = "Pinus pumila", duxiang = "Rhododendron tomentosum",
                   yuejv = "Vaccinium vitis-idaea", dusiyuejv = "Vaccinium uliginosum",
                   xiaohuanghuacai = "Hemerocallis minor", duanbaijinlianhua = "Trollius ledebourii")
species_latin_abbr <- sub("^([A-Z])[a-z]+ ", "\\1. ", species_latin)   # 属名缩写

GREEN <- "#2E7D32"   # 全图统一颜色 (适生绿); 改这里即可

# --- 数据: 边界 + 当前面积(面积升序) ---
boundary_wgs  <- st_transform(st_read(boundary_file, quiet = TRUE), 4326)
area_raw      <- read.csv(file.path(tbl_dir, "area_raw_data.csv"), fileEncoding = "UTF-8")
cur_area      <- area_raw %>% distinct(species_code, current_area_km2) %>% arrange(current_area_km2)
species_order <- cur_area$species_code   # 面积升序

# ----------------------------------------------------------------------------
# Fig6: 2行×6列网格 — 每列一个物种，地图在上、柱子在下，6列并排
# ----------------------------------------------------------------------------
mk_cur_map <- function(sp_code, ttl) {
  r <- as.factor(mask(crop(rast(file.path(output_dir, sp_code, "binary_current_tss.tif")),
                           boundary_wgs), boundary_wgs))
  ggplot() +
    geom_spatraster(data = r) +
    geom_sf(data = boundary_wgs, fill = NA, color = "grey30", linewidth = 0.15) +
    scale_fill_manual(values = c("0" = "grey90", "1" = GREEN), na.value = "transparent",
                      na.translate = FALSE, guide = "none") +
    coord_sf(expand = FALSE) +
    labs(title = ttl) +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(face = "italic", hjust = 0.5, size = 10.5),
          legend.text = element_text(size = 10))
}

# 统一柱形图 y 轴上限
max_area <- max(cur_area$current_area_km2)

# --- Fig6: 每列内部「地图上 + 柱子下」拼接，保证对齐 ---
# y 轴上限取整到万
max_y  <- ceiling(max_area * 1.12 / 10000) * 10000

combined_panels <- purrr::map(species_order, function(code) {
  sp_latin  <- species_latin_abbr[[code]]
  area_val  <- cur_area$current_area_km2[cur_area$species_code == code]
  is_first  <- code == species_order[1]
  y_lab     <- if (is_first) expression("Current suitable area (km"^2*")") else NULL
  bar_df    <- data.frame(x = 1, sp = sp_latin, area = area_val)

  # 地图 (无标题, 无图例)
  map_p <- mk_cur_map(code, "")

  # 柱子 (连续 x 轴, 与 Fig7 等宽; y轴刻度 0-40000; y轴仅首列显示; x轴拉丁名斜体)
  bar_p <- ggplot(bar_df, aes(x, area)) +
    geom_col(width = 0.7, fill = GREEN) +
    geom_text(aes(label = sprintf("%.0f", area)), vjust = -0.5, size = 5) +
    scale_y_continuous(limits = c(0, max_y), expand = c(0, 0),
                       breaks = seq(0, max_y, by = 10000),
                       labels = scales::label_comma()) +
    scale_x_continuous(breaks = 1, labels = sp_latin) +
    labs(x = NULL, y = y_lab) +
    theme_classic(base_size = 12) +
    theme(axis.text.x     = element_text(face = "italic", hjust = 0.5, size = 12),
          axis.text.y     = if (is_first) element_text(size = 12) else element_blank(),
          axis.ticks.y    = if (is_first) element_line() else element_blank(),
          axis.line.y     = if (is_first) element_line() else element_blank(),
          axis.ticks.length = unit(0.15, "cm"))

  # 纵向拼接：地图在上，柱子在下
  map_p / bar_p + plot_layout(heights = c(2, 1))
})

# --- 6 列并排 ---
fig6 <- wrap_plots(combined_panels, nrow = 1)

# --- 输出 PNG，用 grid.text 添加粗体 A / B 标签 ---
png(file.path(fig_dir, "Fig1_CurrentSuitability.png"),
    width = 14, height = 7.5, units = "in", res = 300)
grid::grid.newpage()
grid::grid.draw(fig6)
grid::grid.text("A", x = unit(0.02, "npc"), y = unit(0.97, "npc"),
                just = c("left", "top"), gp = grid::gpar(fontsize = 18))
grid::grid.text("B", x = unit(0.02, "npc"), y = unit(0.46, "npc"),
                just = c("left", "top"), gp = grid::gpar(fontsize = 18))
dev.off()

cat("Fig6 已生成:", file.path(fig_dir, "Fig6_CurrentSuitability.png"), "\n")
