# =============================================================================
# 大兴安岭野生经济植物 MaxEnt 建模 — 发表级工作流 (001b, v5)
#
# 核心方法引用：
#  [ENMeval]  Kass et al. (2021) Methods Ecol Evol 12:1602–1608
#  [maxnet]   Phillips et al. (2017) Ecography 40:887–896
#  [MaxEnt]   Phillips et al. (2006) Ecol Model 190:231–259
#  [BlockCV]  Roberts et al. (2017) Ecography 40:913–929
#  [AICc]     Warren & Seifert (2011) Ecol Appl 21:335–342
#  [OR10]     Radosavljevic & Anderson (2014) J Biogeogr 41:629–643
#  [VIF]      Naimi et al. (2014) Ecography 37:1192–1195 (usdm)
#  [Collin]   Dormann et al. (2013) Ecography 36:27–46
#  [Ref]      Gu et al. (2023) Nat Ecol Evol 7:1652–1665
#
# == v5 改进 ==
# - 添加缓存机制：occurrence_thinned.csv / vif_selection.csv / env_final.rds
# - 预处理步骤（§2-§4）若已有缓存则跳过，大幅缩短重复运行时间
# - refresh_occ / refresh_env 开关控制是否强制刷新
# =============================================================================
rm(list = ls()); gc()

# =============================================================================
# §0  包加载
# =============================================================================
# install.packages(c("ENMeval","maxnet","terra","sf","usdm","ecospat",
#                    "ggplot2","tidyverse","readxl","RColorBrewer"))

  library(ENMeval)       # v2.0+；maxnet 调参标准工具 [ENMeval]
  library(maxnet)        # R 原生 MaxEnt，无 Java 依赖 [maxnet]
  library(terra)         # 栅格操作
  library(sf)            # 矢量操作
  library(usdm)          # VIF 共线性筛选 [VIF]
  library(ecospat)       # MaxTSS 阈值、Boyce 指数计算
  library(pROC)          # 完整训练集 AUC（与 Gu et al. 2023 一致）
  library(ggplot2)
  library(tidyverse)
  library(readxl)
  library(RColorBrewer)

set.seed(2026)

# =============================================================================
# §1  路径与全局常量
# =============================================================================
work_dir            <- "K:/周建maxent"
setwd(work_dir)

climate_current_dir <- paste0(work_dir, "/30s数据/现在气候数据")
climate_future_dir  <- paste0(work_dir, "/30s数据/未来气候数据")
population_dir      <- paste0(work_dir, "/全球1km下尺度人口密度数据集")
boundary_dir        <- paste0(work_dir, "/1矢量边界")
species_file        <- paste0(work_dir,
  "/2野生经济作物坐标点/0原始数据_精准植物调查数据（用于矢量化）.xls")

output_dir <- paste0(work_dir, "/MaxEnt_Fixed_Results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 物种列表：中文名 → 文件夹代码
species_list <- c(
  "偃松"       = "yangsong",
  "杜香"       = "duxiang",
  "越桔"       = "yuejv",
  "笃斯越桔"   = "dusiyuejv",
  "小黄花菜"   = "xiaohuanghuacai",
  "短瓣金莲花" = "duanbaijinlianhua"
)

scenarios <- c("ssp126", "ssp245", "ssp585")
gcm_list  <- c("BCC-CSM2-MR", "MIROC6", "IPSL-CM6A-LR")

no.Cores  <- 16     # ENMeval 并行核数
null.iter <- 100   # ENMnulls 迭代次数；设为 0 可跳过 null 检验
run_nulls <- TRUE  # FALSE = 跳过 null model（节省约 10–30 分钟/物种）

# 气候情景间保持不变的静态变量（elevation/slope/aspect），
# 或按 SSP 替换的变量（pop_density）
STATIC_VARS <- c("elevation", "pop_density", "slope", "aspect")

# =============================================================================
# §1.5  缓存控制
# =============================================================================
# 设为 TRUE 可强制重新运行预处理步骤（边界裁剪 + 空间稀疏化 + VIF 筛选）
# 设为 FALSE 则使用已有的缓存文件（适合调试建模参数时反复测试）
refresh_occ <- FALSE   # 重新运行 §2（分布点裁剪 + 空间稀疏化）
refresh_env <- FALSE   # 强制重新运行 §3-§4（环境变量加载 + VIF 筛选）

# 缓存文件路径
cache_occ_file  <- paste0(output_dir, "/occurrence_thinned.csv")
cache_vif_file  <- paste0(output_dir, "/vif_selection.csv")
cache_env_file  <- paste0(output_dir, "/env_final.rds")

# =============================================================================
# §2  物种分布点读取、裁剪与空间稀疏化（可缓存）
# =============================================================================

if (file.exists(cache_occ_file) && !refresh_occ) {

  # ---- 从缓存加载 ----
  occ_all <- read_csv(cache_occ_file, show_col_types = FALSE) %>% as.data.frame()
  cat(sprintf("【缓存】加载 %s (%d 条, %d 物种)\n",
              cache_occ_file, nrow(occ_all), length(unique(occ_all$species))))
  cat("各物种记录数：\n"); print(table(occ_all$species))

} else {

  # ---- 完整运行 §2 ----
  cat("===== §2 分布点读取与稀疏化 =====\n")
  species_sheets <- excel_sheets(species_file)
  occ_list       <- list()

# sp_name <- species_sheets[1]

  for (sp_name in species_sheets) {
    occ_df <- read_excel(species_file, sheet = sp_name)

    # 先验证 X/Y 列存在，再转数值——防止列名不一致引发静默错误
    if (!all(c("X", "Y") %in% names(occ_df))) {
      warning(sp_name, ": 缺少 X/Y 列，已跳过"); next
    }

    occ_df$lon     <- suppressWarnings(as.numeric(occ_df$X))
    occ_df$lat     <- suppressWarnings(as.numeric(occ_df$Y))
    occ_df$species <- sp_name

    n_raw  <- nrow(occ_df)
    occ_df <- occ_df[!is.na(occ_df$lon) & !is.na(occ_df$lat), ]
    if (n_raw > nrow(occ_df))
      message(sp_name, ": 去除 NA 坐标 ", n_raw - nrow(occ_df), " 行")

    occ_list[[sp_name]] <- occ_df[, c("lon", "lat", "species")]
  }

  occ_all_raw <- bind_rows(occ_list) %>% distinct(species, lon, lat)

  # 研究区边界（大兴安岭），统一 WGS84
  study_boundary <- st_read(paste0(boundary_dir, "/大兴安岭范围准确.shp"), quiet = TRUE) %>%
    st_transform(crs = 4326)

  # 裁剪至研究区
  occ_sf       <- st_as_sf(occ_all_raw, coords = c("lon", "lat"), crs = 4326)
  occ_in_study <- st_intersection(occ_sf, study_boundary)

  # st_intersection 后须从 geometry 重提坐标：
  # 直接用原始 lon/lat 列时，若记录被边界裁切，坐标与空间位置可能错位
  occ_all <- occ_in_study %>%
    mutate(lon = st_coordinates(.)[, 1],
           lat = st_coordinates(.)[, 2]) %>%
    as.data.frame() %>%
    dplyr::select(lon, lat, species)

  # 空间稀疏化（raster-cell thinning）：每个 30 arcsec 格网仅保留 1 条记录，
  # 消除同格网内重复采样引入的空间自相关
  ref_ras <- rast(
    list.files(climate_current_dir, pattern = "wc2.1_30s_bio_\\d+\\.tif$",
               full.names = TRUE)[1]
  )

  occ_thin_list <- list()
  for (sp in unique(occ_all$species)) {
    sp_df   <- occ_all[occ_all$species == sp, ]
    sp_vect <- terra::vect(sp_df, geom = c("lon", "lat"), crs = "EPSG:4326")
    sp_df$cell <- terra::cellFromXY(ref_ras, terra::crds(sp_vect))
    occ_thin_list[[sp]] <- sp_df %>%
      filter(!is.na(cell)) %>%
      group_by(cell) %>%
      slice_sample(n = 1) %>%
      ungroup() %>%
      dplyr::select(lon, lat, species)
  }

  occ_all <- bind_rows(occ_thin_list)

  # 稀疏化后各物种记录数（供方法描述参考）
  cat("稀疏化后各物种记录数：\n"); print(table(occ_all$species))

  write_csv(occ_all, cache_occ_file)
  cat(sprintf("【缓存已保存】%s\n", cache_occ_file))
}

# ============================================================================
# §3–§4  环境变量加载 / VIF 筛选（带 .tif 缓存）
# ============================================================================

# ---- 缓存路径定义 ----
env_final_tif <- paste0(output_dir, "/env_final.tif")
elev_tif      <- paste0(output_dir, "/elev_masked.tif")
slope_tif     <- paste0(output_dir, "/slope_masked.tif")
aspect_tif    <- paste0(output_dir, "/aspect_masked.tif")

# cache_env_file: .rds（仅存 study_boundary 等 sf 对象）
# cache_vif_file: .csv（VIF 筛选结果表）

cache_ready <- all(file.exists(c(
  cache_env_file, cache_vif_file,
  env_final_tif, elev_tif, slope_tif, aspect_tif
)))

if (cache_ready && !refresh_env) {

  # ===========================================================================
  # 从缓存加载（SpatRaster 一律走 .tif，sf 走 .rds）
  # ===========================================================================
  env_final     <- rast(env_final_tif)
  elev_masked   <- rast(elev_tif)
  slope_masked  <- rast(slope_tif)
  aspect_masked <- rast(aspect_tif)

  vif_table   <- read_csv(cache_vif_file, show_col_types = FALSE)
  kept_vars   <- vif_table %>% filter(status == "保留") %>% pull(variable)
  kept_bio    <- setdiff(kept_vars, STATIC_VARS)
  kept_static <- intersect(kept_vars, STATIC_VARS)

  # study_boundary（sf 对象，无指针问题）
  tmp_env        <- readRDS(cache_env_file)
  study_boundary <- tmp_env$study_boundary
  rm(tmp_env); gc()

  cat(sprintf("【缓存】加载 env_final.tif (%d 层)\n", nlyr(env_final)))
  cat("保留变量:", paste(kept_vars, collapse = ", "), "\n")

  # 背景点采样池
  env_final_complete <- app(env_final, function(x) as.integer(sum(is.na(x)) == 0))
  valid_bg_cells     <- which(values(env_final_complete) == 1)

} else {

  # ===========================================================================
  # 完整运行 §3–§4
  # ===========================================================================
  cat("===== §3 环境变量加载与掩膜 =====\n")

  # --------------------------------------------------------------------------
  # 3.1  BIO1–19（WorldClim 2.1, 30 arcsec）
  # --------------------------------------------------------------------------
  bio_files <- list.files(climate_current_dir,
                          pattern = "wc2.1_30s_bio_\\d+\\.tif$", full.names = TRUE)
  bio_nums  <- as.numeric(gsub(".*bio_(\\d+)\\.tif$", "\\1", basename(bio_files)))
  bio_files <- bio_files[order(bio_nums)]
  bio_stack <- rast(bio_files)
  names(bio_stack) <- paste0("bio", sort(bio_nums))

  # 边界

  if (!exists("study_boundary")) {
    study_boundary <- st_read(paste0(boundary_dir, "/大兴安岭范围准确.shp"), quiet = TRUE) %>%
      st_transform(crs = 4326)
  }
  bio_masked <- mask(crop(bio_stack, study_boundary), study_boundary)

  # --------------------------------------------------------------------------
  # 3.2  高程与地形
  # --------------------------------------------------------------------------
  elev_file <- list.files(climate_current_dir,
                          pattern = "wc2.1_30s_elev\\.tif$", full.names = TRUE)
  elev   <- rast(elev_file); names(elev) <- "elevation"
  slope  <- terrain(elev, v = "slope",  unit = "degrees"); names(slope)  <- "slope"
  aspect <- terrain(elev, v = "aspect", unit = "degrees"); names(aspect) <- "aspect"

  elev_masked   <- mask(crop(elev,   study_boundary), study_boundary)
  slope_masked  <- mask(crop(slope,  study_boundary), study_boundary)
  aspect_masked <- mask(crop(aspect, study_boundary), study_boundary)

  # --------------------------------------------------------------------------
  # 3.3  人口密度（bilinear 重采样，连续变量）
  # --------------------------------------------------------------------------
  pop_2000_file <- paste0(population_dir,
                          "/BaseYear_2000_1km-geotiff/BaseYear_1km/baseYr_total_2000.tif")
  pop_current   <- rast(pop_2000_file)
  pop_current   <- resample(pop_current, bio_stack, method = "bilinear")
  pop_masked    <- mask(crop(pop_current, study_boundary), study_boundary)
  names(pop_masked) <- "pop_density"

  # --------------------------------------------------------------------------
  # 3.4  合并全部候选变量
  # --------------------------------------------------------------------------
  env_vars <- c(bio_masked, elev_masked, pop_masked, slope_masked, aspect_masked)

  # ==========================================================================
  # §4  VIF 共线性筛选
  # ==========================================================================
  cat("===== §4 VIF 共线性筛选 =====\n")

  # 构建完整有效格网掩膜（所有层均非 NA）

  env_vars_complete  <- app(env_vars, function(x) as.integer(sum(is.na(x)) == 0))
  complete_cells_vif <- which(values(env_vars_complete) == 1)

  # 采样（固定种子，可复现）
  set.seed(2026)
  sample_size <- min(10000, length(complete_cells_vif))
  sample_cells <- sample(complete_cells_vif, size = sample_size, replace = FALSE)

  vif_pts <- terra::xyFromCell(env_vars[[1]], sample_cells)
  vif_sv  <- terra::vect(vif_pts, crs = terra::crs(env_vars))

  # terra::extract —— 安全写法（兼容不同 terra 版本）
  vif_env <- terra::extract(env_vars, vif_sv)
  if ("ID" %in% names(vif_env)) vif_env$ID <- NULL
  vif_env <- na.omit(vif_env)

  # Step 1: vifcor (|r| < 0.7)
  vc <- vifcor(vif_env, th = 0.7)
  # Step 2: vifstep (VIF < 10)
  vif_result <- vifstep(vif_env[, vc@results$Variables], th = 10)

  kept_vars     <- as.character(vif_result@results$Variables)
  excluded_vars <- as.character(vif_result@excluded)

  cat("\nVIF 筛选结果：\n")
  cat("  保留 (", length(kept_vars), "):", paste(kept_vars, collapse = ", "), "\n")
  cat("  剔除 (", length(excluded_vars), "):", paste(excluded_vars, collapse = ", "), "\n\n")

  # 导出 VIF 表
  vif_out <- data.frame(
    variable = c(kept_vars, excluded_vars),
    VIF      = c(round(vif_result@results$VIF, 2),
                 rep(NA, length(excluded_vars))),
    status   = c(rep("保留", length(kept_vars)),
                 rep("剔除", length(excluded_vars)))
  )
  write_csv(vif_out, cache_vif_file)
  cat(sprintf("【缓存已保存】%s\n", cache_vif_file))

  # --------------------------------------------------------------------------
  # 构建最终环境栅格 env_final（统一完整掩膜）
  # --------------------------------------------------------------------------
  env_final <- env_vars[[kept_vars]]

  mask_all <- app(env_final, function(x) as.integer(sum(is.na(x)) == 0))
  mask_all[mask_all == 0] <- NA
  env_final <- mask(env_final, mask_all)

  # --------------------------------------------------------------------------
  # 写出 .tif（SpatRaster 持久化唯一可靠方式）
  # --------------------------------------------------------------------------
  terra::writeRaster(env_final,     env_final_tif, overwrite = TRUE)
  terra::writeRaster(elev_masked,   elev_tif,      overwrite = TRUE)
  terra::writeRaster(slope_masked,  slope_tif,     overwrite = TRUE)
  terra::writeRaster(aspect_masked, aspect_tif,    overwrite = TRUE)
  cat(sprintf("【缓存已保存】env_final.tif (%d 层) + 地形 .tif\n", nlyr(env_final)))

  # 写出 .rds（仅 sf / 标量，不含 SpatRaster）
  env_cache <- list(study_boundary = study_boundary)
  saveRDS(env_cache, cache_env_file)
  cat(sprintf("【缓存已保存】%s (study_boundary)\n", cache_env_file))

  # --------------------------------------------------------------------------
  # 后续变量分组
  # --------------------------------------------------------------------------
  kept_bio    <- setdiff(kept_vars, STATIC_VARS)
  kept_static <- intersect(kept_vars, STATIC_VARS)

  # 背景点采样池
  env_final_complete <- app(env_final, function(x) as.integer(sum(is.na(x)) == 0))
  valid_bg_cells     <- which(values(env_final_complete) == 1)
}

# =============================================================================
# §5  ENMeval 调参 + maxnet 建模（逐物种循环）
# =============================================================================

part_info <- list(method = "block", settings = list(orientation = "lat_lon"))

# ENMeval vignette 推荐：以出现点为中心做等积投影缓冲区来定义背景范围，
# 避免将整个研究区都纳入背景（可能引入非目标物种实际可达的环境空间）
bg_buffer_km <- 50  # 背景点缓冲半径（km），可按物种扩散能力调整

model_results <- data.frame()
all_models    <- list()

for (sp_name in names(species_list)) {

  sp_code <- species_list[[sp_name]]
  sp_dir  <- paste0(output_dir, "/", sp_code)
  dir.create(sp_dir, showWarnings = FALSE, recursive = TRUE)
  cat("===== 物种:", sp_name, "=====\n")

  # 提取分布点，再做一次栅格去重（occ_all 已经稀疏化，此处是安全保障）
  sp_data  <- occ_all %>% filter(species == sp_name) %>% dplyr::select(lon, lat)
  cells_sp <- terra::cellFromXY(env_final[[1]], as.matrix(sp_data))
  sp_data  <- sp_data[!duplicated(cells_sp) & !is.na(cells_sp), ] %>% as.data.frame()

  if (nrow(sp_data) < 10) {
    warning(sp_name, ": 有效分布点 ", nrow(sp_data), " < 10，跳过"); next
  }

  # 背景点：ENMeval vignette 推荐 —— 以出现点为中心做缓冲区，仅在缓冲区内采样背景点
  # 步骤：WGS84 → 等积投影（Eckert IV）缓冲 → 合并 → 转回 WGS84 → 与研究边界取交集
  # 等积投影确保 km 单位的缓冲距离在高纬度地区不失真
  eckertIV    <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m"
  occs_sf_sp  <- sf::st_as_sf(sp_data, coords = c("lon", "lat"), crs = 4326)
  occs_buf_sp <- sf::st_transform(occs_sf_sp, crs = eckertIV) %>%
    sf::st_buffer(dist = bg_buffer_km * 1000) %>%
    sf::st_union() %>% sf::st_sf() %>%
    sf::st_transform(crs = 4326) %>%
    sf::st_intersection(study_boundary)

  # 将缓冲区裁剪到 env_final（确保只在有环境数据的区域采样背景点）
  envs_bg_sp <- terra::mask(terra::crop(env_final, occs_buf_sp), occs_buf_sp)

  # terra::spatSample() 自动排除 NA 格网，比手动筛选 valid_cells 更稳健
  set.seed(2026)
  bg_df <- terra::spatSample(envs_bg_sp, size = 10000, na.rm = TRUE,
                              values = FALSE, xy = TRUE) %>% as.data.frame()
  colnames(bg_df) <- c("lon", "lat")
  cat("  背景点缓冲区内有效采样:", nrow(bg_df), "个\n")

# # --- 图1：出现点 + 缓冲区 ---
# plot(envs_bg_sp[[1]], main = paste0(sp_name, " — ", names(envs_bg_sp)[1]))
# points(sp_data$lon, sp_data$lat, pch = 16, col = "red", cex = 1)
# plot(st_geometry(occs_buf_sp), border = "blue", lwd = 2, add = TRUE)

# # --- 图2：加上背景点 ---
# plot(envs_bg_sp[[1]], main = paste0(sp_name, " — Occ + BG"))
# points(bg_df$lon, bg_df$lat, pch = 1, col = "grey50", cex = 0.3)
# points(sp_data$lon, sp_data$lat, pch = 16, col = "red", cex = 1.2)
# plot(st_geometry(occs_buf_sp), border = "blue", lwd = 2, add = TRUE)
# legend("topright", legend = c("Occurrence", "Background", "Buffer"),
#        col = c("red", "grey50", "blue"), pch = c(16, 1, NA),
#        lty = c(NA, NA, 1), lwd = c(NA, NA, 2), cex = 0.8)

# # --- 图3：空间分块着色（运行 get.block 后）---
#   occs = sp_data[, c("lon", "lat")],
#   bg   = bg_df[, c("lon", "lat")],
#   orientation = "lat_lon"
# )

# # 检查一下结果
# str(block_part)

# # --- 然后再画图 ---
# fold_colors <- c("red", "blue", "darkgreen", "purple")
# plot(envs_bg_sp[[1]], main = paste0(sp_name, " — Block Partition"))
# points(bg_df$lon, bg_df$lat,
#        pch = 1, cex = 0.3,
#        col = adjustcolor(fold_colors[block_part$bg.grp], alpha.f = 0.3))
# points(sp_data$lon, sp_data$lat,
#        pch = 17, cex = 1.5,
#        col = fold_colors[block_part$occs.grp])
# plot(sf::st_geometry(occs_buf_sp), border = "blue", lwd = 2, add = TRUE)
# legend("topright", legend = paste("Fold", 1:4),
#        col = fold_colors, pch = 17, cex = 0.8)

  # p10_n：p10 阈值对应的排名位置
  # 出现点 < 10 时用 floor（向下取整），≥ 10 时用 ceiling（向上取整）
  # 详见 Radosavljevic & Anderson (2014) [OR10]
  p10_n <- ifelse(nrow(sp_data) < 10,
                  floor(nrow(sp_data) * 0.9),
                  ceiling(nrow(sp_data) * 0.9))

  cat("  分区方法: block (lat_lon, 4-fold spatial partition)\n")

  tryCatch({

    # ------------------------------------------------------------------
    # §5.1  ENMeval 调参（自适应降级）
    #
    # 特征类型（fc）：L / LQ / H / LQH，与 Gu et al. (2023) 一致
    #   L   = Linear（仅线性特征，最简单）
    #   LQ  = Linear + Quadratic（允许单峰响应）
    #   H   = Hinge only（对少点物种易致矩阵奇异→自动降级移除）
    #   LQH = Linear + Quadratic + Hinge（最常用的综合组合）
    #   共 4 类 × 7 个 rm = 28 个候选参数组合
    # 正则化乘数（rm）：1–4，步长 0.5，共 7 值
    #   rm 越大模型越平滑，过小易过拟合 [Radosavljevic 2014]
    # 空间分区：block，4 折，lat_lon 方向
    #   block CV 保证训练集与测试集空间独立，避免乐观偏差 [Roberts 2017]
    #
    # 降级策略：若 full fc（含 H）因 non-conformable 崩溃，
    # 自动去掉纯 H 重试（LQH 仍保留），保证窄分布物种也能出结果。
    # ------------------------------------------------------------------
    fc_full  <- c("L", "LQ", "H", "LQH")
    fc_safe  <- c("L", "LQ", "LQH")    # 去掉纯 H（铰链过多致矩阵奇异）

    # ---- 第 1 轮：full fc 并行 ----
    enm_eval  <- NULL
    fc_actual <- fc_full

    # envs 传入缓冲区裁剪后的环境栅格（envs_bg_sp），与 bg_df 背景点范围一致
    # vignette 推荐：ENMevaluate 的训练环境范围应与背景点来源范围保持一致
    enm_eval <- tryCatch({
      ENMeval::ENMevaluate(
        occs = sp_data, envs = envs_bg_sp, bg = bg_df,
        algorithm = "maxnet", partitions = part_info$method,
        tune.args = list(rm = seq(0.5, 4, 0.5), fc = fc_actual),
        partition.settings = part_info$settings,
        parallel = TRUE, numCores = no.Cores,
        quiet = FALSE, doClamp = TRUE
      )
    }, error = function(e1) {
      msg1 <- conditionMessage(e1)
      # 若含 H 导致的矩阵错误 → 去掉 H 重试
      if (grepl("non-conformable", msg1)) {
        fc_actual <<- fc_safe
        cat("  fc=H 致矩阵奇异，降级为 L/LQ/LQH\n")
        gc()
        tryCatch({
          ENMeval::ENMevaluate(
            occs = sp_data, envs = envs_bg_sp, bg = bg_df,
            algorithm = "maxnet", partitions = part_info$method,
            tune.args = list(rm = seq(0.5, 4, 0.5), fc = fc_actual),
            partition.settings = part_info$settings,
            parallel = TRUE, numCores = no.Cores,
            quiet = FALSE, doClamp = TRUE
          )
        }, error = function(e2) {
          # 降级后仍失败 → 退串行
          cat(sprintf("  降级并行仍失败（%s），退串行\n", conditionMessage(e2)))
          gc()
          ENMeval::ENMevaluate(
            occs = sp_data, envs = envs_bg_sp, bg = bg_df,
            algorithm = "maxnet", partitions = part_info$method,
            tune.args = list(rm = seq(0.5, 4, 0.5), fc = fc_actual),
            partition.settings = part_info$settings,
            parallel = FALSE, quiet = FALSE, doClamp = TRUE
          )
        })
      } else {
        # 非矩阵错误（内存/端口竞争）→ 直接用 full fc 退串行
        cat(sprintf("  并行调参失败（%s），退回串行\n", msg1))
        gc()
        ENMeval::ENMevaluate(
          occs = sp_data, envs = envs_bg_sp, bg = bg_df,
          algorithm = "maxnet", partitions = part_info$method,
          tune.args = list(rm = seq(0.5, 4, 0.5), fc = fc_full),
          partition.settings = part_info$settings,
          parallel = FALSE, quiet = FALSE, doClamp = TRUE
        )
      }
    })

    res <- enm_eval@results
    write_csv(res, paste0(sp_dir, "/enmeval_results_all.csv"))

    # ------------------------------------------------------------------
    # §5.3  模型选优 —— ENMeval sequential 准则
    #       (Kass et al. 2021, ENMeval 2.0 vignette)
    #
    # Step 1: 选 min(or.10p.avg) —— 遗漏率最低的模型
    # Step 2: 平局 → max(cbi.val.avg) —— 最高验证 CBI
    # ------------------------------------------------------------------

    res_valid <- res %>% filter(!is.na(or.10p.avg), !is.na(cbi.val.avg))
    if (nrow(res_valid) == 0) res_valid <- res

    # Step 1 + Step 2：一体化 sequential 筛选
    best_row <- res_valid %>%
      filter(or.10p.avg == min(or.10p.avg)) %>%
      filter(cbi.val.avg == max(cbi.val.avg)) %>%
      slice(1)   # 极端平局兜底：取第一行

    best_fc      <- as.character(best_row$fc)
    best_rm      <- as.numeric(best_row$rm)
    # tune.args 列（如 "fc.LQH_rm.2"）是 ENMeval 内部索引键，
    # 用于从 eval.models() / eval.predictions() 中精确定位最优模型
    best_tune_key <- as.character(best_row$tune.args)

    cat("  最优参数: FC =", best_fc, "| RM =", best_rm,
        "| OR10_cv =",   round(best_row$or.10p.avg,  4),
        "| AUC_cv =",    round(best_row$auc.val.avg,  4),
        "| CBI_cv =",    round(best_row$cbi.val.avg,  4), "\n")

    auc_diff <- if ("auc.diff.avg" %in% names(best_row))
                  round(as.numeric(best_row$auc.diff.avg), 4) else NA_real_

    # ------------------------------------------------------------------
    # §5.4  Null model 显著性检验（可选）
    #
    # ENMnulls 将出现点随机置换 null.iter 次，每次用相同参数建模并评估，
    # 得到 null 分布的 OR10 和 CBI；若真实值显著优于 null 分布，
    # 说明模型捕捉到了真实的生态信号，而非随机地理格局
    # p(OR10) = null 模型中 OR10 ≤ 实际值的比例（越小越好，期望 < 0.05）
    # p(CBI)  = null 模型中 CBI ≥ 实际值的比例（越小越好，期望 < 0.05）
    # ------------------------------------------------------------------
    null_pval_or <- null_pval_cbi <- NA_real_

    if (run_nulls && null.iter > 0) {
      cat("  运行 Null model 检验（", null.iter, "次迭代）...\n")
      tryCatch({
        null_mod <- ENMeval::ENMnulls(
          e            = enm_eval,
          mod.settings = list(fc = best_fc, rm = best_rm),
          no.iter      = null.iter,
          parallel     = TRUE,
          numCores     = no.Cores,
          quiet        = TRUE
        )

        null_res      <- null_mod@null.results
        null_pval_or  <- mean(null_res$or.10p.avg  <= best_row$or.10p.avg,  na.rm = TRUE)
        null_pval_cbi <- mean(null_res$cbi.val.avg >= best_row$cbi.val.avg, na.rm = TRUE)
        cat("  Null model: p(OR10) =", round(null_pval_or, 3),
            "| p(CBI) =",             round(null_pval_cbi, 3), "\n")

        saveRDS(null_mod, paste0(sp_dir, "/null_model.rds"))
      }, error = function(e) {
        warning("Null model 失败（", conditionMessage(e), "）")
      })
    }

    # ------------------------------------------------------------------
    # §5.5  提取训练数据表（供 §5.9 全训练集评估使用）
    #
    # 从 env_final（全研究区）提取出现点 + 背景点的环境值，
    # bg_df 来自缓冲区（envs_bg_sp ? env_final），两者在相同位置的值一致，
    # 使用 SpatVector 路径：vect() → extract(ID=FALSE)
    # ------------------------------------------------------------------
    pres_sv   <- terra::vect(as.data.frame(sp_data),
                              geom = c("lon", "lat"), crs = terra::crs(env_final))
    pres_vals <- terra::extract(env_final, pres_sv, ID = FALSE)
    pres_swd  <- na.omit(cbind(as.data.frame(sp_data), pres_vals))

    bg_sv   <- terra::vect(bg_df, geom = c("lon", "lat"), crs = terra::crs(env_final))
    bg_vals <- terra::extract(env_final, bg_sv, ID = FALSE)
    bg_swd  <- na.omit(cbind(bg_df, bg_vals))

    # all_p：1 = 出现点，0 = 背景点（maxnet 格式）
    all_p    <- c(rep(1L, nrow(pres_swd)), rep(0L, nrow(bg_swd)))
    all_data <- rbind(pres_swd[, kept_vars], bg_swd[, kept_vars])

    # ------------------------------------------------------------------
    # §5.6  提取最终 maxnet 模型（ENMeval vignette 推荐）
    #
    # ENMevaluate() 内部已用全部出现点 + 背景点为每个参数组合训练了完整模型，
    # 直接用 eval.models()[[best_tune_key]] 提取，避免用相同参数重复训练，
    # 保证预测与评估指标来自同一模型对象（与 vignette §"Final model" 一致）
    # ------------------------------------------------------------------
    final_model <- ENMeval::eval.models(enm_eval)[[best_tune_key]]
    saveRDS(final_model, paste0(sp_dir, "/final_model_maxnet.rds"))

    # ------------------------------------------------------------------
    # §5.7  变量重要性（排列重要性,作废）
    #
    # eval.variable.importance() 基于 ENMeval 内部的最优折模型，
    # 对每个变量随机置换后重新预测，记录 AUC 下降量作为重要性得分
    # 结果是"置换重要性"而非线性系数，适用于 MaxEnt 等非参数模型
    # ------------------------------------------------------------------
    # tryCatch({
    #   vi_df   <- ENMeval::eval.variable.importance(enm_eval)
    #   vi_best <- vi_df %>%
    #     filter(fc == best_fc, rm == best_rm) %>%
    #     arrange(desc(permutation.importance))
    #   write_csv(vi_best, paste0(sp_dir, "/variable_importance.csv"))

    # }, error = function(e) NULL)

    # ------------------------------------------------------------------
    # §5.8  响应曲线（边际响应）→ PNG
    #
    # ENMeval vignette 推荐：直接调用 maxnet::response.plot()，
    # 该函数利用模型内部存储的训练数据自动计算边际响应，
    # 比手动固定中位数的方式更准确（基于训练样本真实分布）
    # ------------------------------------------------------------------
    tryCatch({
      rc_df <- bind_rows(lapply(kept_vars, function(v) {
        mrc <- maxnet::response.plot(
          x    = final_model,
          v    = v,
          type = "cloglog",
          plot = FALSE
        )
        data.frame(variable = v, value = mrc[[v]], suitability = mrc$pred)
      }))

    }, error = function(e) NULL)

    # ------------------------------------------------------------------
    # §5.9  完整训练集性能评估 + 阈值计算
    #
    # 对应 Gu et al. (2023) / R1_bird.R 的 SDMs.metrics.full：
    # 在全部出现点 + 背景点上预测，计算：
    #   AUC_full  = pROC::auc（训练集，存在乐观偏差，仅作参考）
    #   TSS_full  = MaxTSS（ecospat.max.tss）
    #   Boyce_full= 连续 Boyce 指数（ecospat.boyce，仅用出现点预测值）
    #   p10       = 训练出现点预测值的第 10 百分位阈值（用于二值化）
    #   MaxTSS_thr= 使 TSS 最大化的概率阈值（用于二值化）
    #
    # 注：论文主要评价指标应引用 CV 版本（AUC_cv / OR10_cv / CBI_cv），
    #     full 版本受训练集自评影响，但与 Gu et al. 2023 报告口径一致
    # ------------------------------------------------------------------
    pred_all  <- as.numeric(predict(final_model, all_data,              type = "cloglog"))
    pred_pres <- as.numeric(predict(final_model, pres_swd[, kept_vars], type = "cloglog"))

    # 完整训练集 AUC（与 R1_bird.R 第 246 行一致）
    auc_full <- tryCatch(
      as.numeric(pROC::auc(pROC::roc(response  = all_p,
                                     predictor = pred_all,
                                     quiet     = TRUE))),
      error = function(e) NA_real_
    )

    # MaxTSS 及其对应阈值（与 R1_bird.R 第 247/250 行一致）
    tss_res    <- ecospat::ecospat.max.tss(Pred = pred_all, Sp.occ = all_p)
    tss_full   <- as.numeric(tss_res$max.TSS)
    tss_thresh <- as.numeric(tss_res$max.threshold)

    # 连续 Boyce 指数（与 R1_bird.R 第 248 行一致）
    # fit = 全体数据预测值，obs = 仅出现点预测值；PEplot=FALSE 不弹图
    boyce_full <- tryCatch(
      as.numeric(ecospat::ecospat.boyce(fit    = pred_all,
                                        obs    = pred_pres,
                                        PEplot = FALSE)$cor),
      error = function(e) NA_real_
    )

    # p10 阈值（与 R1_bird.R 第 249 行一致）
    p10_thresh <- as.numeric(rev(sort(pred_pres))[p10_n])

    # ------------------------------------------------------------------
    # §5.10  性能汇总表
    #
    # 分两组指标，来源不同，用途不同：
    #
    # [CV 指标] 来自 ENMeval 空间块交叉验证，空间独立，无训练偏差，
    #   适合作为论文主要评价指标（Methods 和 Table 中报告）：
    #   AUC_cv   = 验证集平均 AUC；AUC_diff = 过拟合程度（越小越好）
    #   OR10_cv  = CV 遗漏率（期望 ≤ 0.1）；CBI_cv = CV Boyce 指数
    #
    # [Full 指标] 在全部训练数据上预测，存在乐观偏差，
    #   与 Gu et al. (2023) / R1_bird.R 报告口径一致（可用于对比验证）：
    #   AUC_full = 训练集 AUC；TSS_full = 最大 TSS；Boyce_full = Boyce 指数
    #
    # [阈值] 仅用于栅格二值化，不纳入模型评价：
    #   p10_threshold = 出现点预测值第 10 百分位
    #   tss_threshold = MaxTSS 对应概率切点
    # ------------------------------------------------------------------
    perf_df <- data.frame(
      species        = sp_name,
      n_occ          = nrow(sp_data),
      best_FC        = best_fc,
      best_RM        = best_rm,
      # CV 指标（ENMeval block CV
      AUC_cv         = round(as.numeric(best_row$auc.val.avg),   4),
      AUC_diff       = auc_diff,
      OR10_cv        = round(as.numeric(best_row$or.10p.avg),    4),
      CBI_cv         = round(as.numeric(best_row$cbi.val.avg),   4),
      AICc           = round(as.numeric(best_row$AICc),           2),
      delta_AICc     = round(as.numeric(best_row$delta.AICc),     2),
      null_p_OR10    = round(null_pval_or,  3),
      null_p_CBI     = round(null_pval_cbi, 3),
      # Full 指标（全训练集，与 Gu et al. 2023 口径一致
      AUC_full       = round(auc_full,   4),
      TSS_full       = round(tss_full,   4),
      Boyce_full     = round(boyce_full, 4),
      # 二值化阈值（训练集计算
      p10_threshold  = round(p10_thresh,  4),
      tss_threshold  = round(tss_thresh,  4)
    )

    model_results <- bind_rows(model_results, perf_df)
    all_models[[sp_name]] <- list(
      model       = final_model,
      performance = perf_df,
      enm_eval    = enm_eval
    )

    cat("  [CV]   AUC_cv:", round(perf_df$AUC_cv, 3),
        "| OR10_cv:", round(perf_df$OR10_cv, 3),
        "| CBI_cv:",  round(perf_df$CBI_cv,  3),
        "| AUC_diff:", round(perf_df$AUC_diff, 3), "\n")
    cat("  [Full] AUC_full:", round(perf_df$AUC_full, 3),
        "| TSS_full:", round(perf_df$TSS_full, 3),
        "| Boyce_full:", round(perf_df$Boyce_full, 3), "\n\n")

  }, error = function(e) {
    # 去掉 tryCatch 后可获得完整报错堆栈，便于单物种调试
    cat(sprintf("  ? %s 建模失败: %s\n", sp_name, conditionMessage(e)))
  })

  gc()
}

# 汇总表写盘，并打印供核查
if (nrow(model_results) > 0) {
  write_csv(model_results, paste0(output_dir, "/model_performance_summary.csv"))
  # CV 指标（主表）：AUC_cv / OR10_cv / CBI_cv / AUC_diff + null model p 值
  # Full 指标（补充）：AUC_full / TSS_full / Boyce_full（与 Gu et al. 2023 一致）
  show_cols <- intersect(
    c("species", "n_occ", "best_FC", "best_RM",
      "AUC_cv", "OR10_cv", "CBI_cv", "AUC_diff",
      "null_p_OR10", "null_p_CBI",
      "AUC_full", "TSS_full", "Boyce_full"),
    names(model_results)
  )
  cat("=== 模型性能汇总 ===\n")
  print(model_results[, show_cols], row.names = FALSE)
} else {
  stop("所有物种均建模失败。\n",
       "调试方法：以 sp_name <- '偃松' 单独运行 §5 循环体（去掉外层 tryCatch）\n",
       "或运行 001b_test_species_debug.R 获取详细错误信息")
}

# =============================================================================
# §6  当前 + 未来分布预测（3 SSP × 3 GCM）
# =============================================================================
if (length(all_models) == 0)
  stop("all_models 为空，请确认 §5 已成功运行。")

# 二值写盘工具函数
# INT1U（1 字节无符号整型）存储 0/1，NAflag=255 明确区分 NA 与 0，
# 便于后续用栅格运算（+、|）合并多情景二值图
write_binary <- function(cont, thr, file) {
  terra::writeRaster(cont >= thr, file,
                     overwrite = TRUE, datatype = "INT1U", NAflag = 255)
}

# §6.0  预加载未来人口密度（SSP 对应 2050 年情景）
# 文件命名规则：ssp1 → SSP1-2.6，ssp2 → SSP2-4.5，ssp5 → SSP5-8.5
ssp_to_pop   <- c(ssp126 = "ssp1", ssp245 = "ssp2", ssp585 = "ssp5")
pop_fut_list <- list()
for (ssp in scenarios) {
  pop_ssp  <- ssp_to_pop[[ssp]]
  pop_file <- paste0(population_dir, "/2050-", pop_ssp, "-geotiff/",
                     pop_ssp, "_total_2050.tif")
  if (file.exists(pop_file)) {
    pop_r <- rast(pop_file)
    pop_r <- resample(pop_r, env_final[[1]], method = "bilinear")
    pop_r <- mask(crop(pop_r, study_boundary), study_boundary)
    names(pop_r) <- "pop_density"
    pop_fut_list[[ssp]] <- pop_r
  } else {
    warning("未来人口文件不存在，将跳过含 pop_density 的物种预测：", pop_file)
    pop_fut_list[[ssp]] <- NULL
  }
}

# §6.1  预测主循环
for (sp_name in names(all_models)) {

  sp_code     <- species_list[[sp_name]]
  sp_dir      <- paste0(output_dir, "/", sp_code)
  final_model <- all_models[[sp_name]]$model
  perf        <- all_models[[sp_name]]$performance

  cat("预测物种:", sp_name, "\n")

  # -- 当前气候预测 --
  # maxnet.predictRaster() 为 ENMeval 官方预测函数，
  # 确保特征构建与训练时完全一致；doClamp=TRUE 将预测值钳制在训练范围内，
  # 防止外推区域产生异常高/低适宜性（在气候情景预测中尤为重要）
  pred_current <- ENMeval::maxnet.predictRaster(
    mod       = final_model,
    envs      = env_final[[kept_vars]],
    pred.type = "cloglog",
    doClamp   = TRUE
  )
  names(pred_current) <- "suitability"

  terra::writeRaster(pred_current,
                     paste0(sp_dir, "/continuous_current.tif"), overwrite = TRUE)
  write_binary(pred_current, perf$p10_threshold,
               paste0(sp_dir, "/binary_current_p10.tif"))
  write_binary(pred_current, perf$tss_threshold,
               paste0(sp_dir, "/binary_current_tss.tif"))

  # -- 未来气候预测：3 SSP × 3 GCM --
  for (ssp in scenarios) {
    for (gcm in gcm_list) {

      fut_file <- paste0(climate_future_dir,
                         "/wc2.1_30s_bioc_", gcm, "_", ssp, "_2041-2060.tif")
      if (!file.exists(fut_file)) {
        warning("气候文件不存在，跳过：", basename(fut_file)); next
      }
      if ("pop_density" %in% kept_vars && is.null(pop_fut_list[[ssp]])) {
        warning("缺少 ", ssp, " 人口数据，跳过 ", gcm, "_", ssp); next
      }

      # 图层命名：WorldClim 2.1 未来数据固定 19 波段，按位置重命名最稳健；
      # 若波段数非 19（罕见），尝试从图层名末尾数字解析
      env_fut <- rast(fut_file)
      if (nlyr(env_fut) == 19) {
        names(env_fut) <- paste0("bio", 1:19)
      } else {
        names(env_fut) <- sub(".*_(\\d+)$", "bio\\1", names(env_fut))
      }

      missing_bio <- setdiff(kept_bio, names(env_fut))
      if (length(missing_bio) > 0) {
        warning("缺少气候层 ", paste(missing_bio, collapse = ", "),
                "，跳过 ", gcm, "_", ssp); next
      }

      # 裁剪气候层 + 完整掩膜（与训练集保持一致的 NA 处理）
      env_fut_masked <- mask(crop(env_fut[[kept_bio]], study_boundary), study_boundary)
      m_fut          <- app(env_fut_masked, function(x) as.integer(sum(is.na(x)) == 0))
      m_fut[m_fut == 0] <- NA
      env_fut_masked <- mask(env_fut_masked, m_fut)

      # 拼接静态变量（elevation/slope/aspect 跨情景不变；pop_density 用 SSP 版本）
      parts <- list(env_fut_masked)
      for (sv in kept_static) {
        layer <- switch(sv,
          "elevation"   = elev_masked,
          "slope"       = slope_masked,
          "aspect"      = aspect_masked,
          "pop_density" = pop_fut_list[[ssp]],
          NULL
        )
        if (!is.null(layer)) parts <- c(parts, list(layer))
      }
      # do.call(c, parts)[[kept_vars]] 确保层顺序与训练时严格一致
      env_fut_final <- do.call(c, parts)[[kept_vars]]

      pred_future <- ENMeval::maxnet.predictRaster(
        mod       = final_model,
        envs      = env_fut_final,
        pred.type = "cloglog",
        doClamp   = TRUE
      )
      names(pred_future) <- "suitability"

      proj_id <- paste0(gcm, "_", ssp)
      terra::writeRaster(pred_future,
                         paste0(sp_dir, "/continuous_", proj_id, ".tif"),
                         overwrite = TRUE)
      write_binary(pred_future, perf$p10_threshold,
                   paste0(sp_dir, "/binary_", proj_id, "_p10.tif"))
      write_binary(pred_future, perf$tss_threshold,
                   paste0(sp_dir, "/binary_", proj_id, "_tss.tif"))
    }
  }

  cat("  完成\n")
  gc()
}

cat("\n全部预测完成。输出目录：", output_dir, "\n")
