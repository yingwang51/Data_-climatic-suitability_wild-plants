# ============================================================================
# 004.GCM_Replicate_Analysis.R
# 大兴安岭6种野生经济植物 — 投影不确定性方差分解
#
# 分析设计 (两层互补):
#   - 物种 (6) × SSP情景 (3) × GCM (3) 全交叉, 单元内 n=1, 共 54 个集合投影输出
#   Tier 1 — 主推断 (物种差异 / 情景差异)  [Step 5]:
#     重复测量ANOVA (afex::aov_ez), GCM 作随机重复单元; 物种↔物种×GCM、
#     情景↔情景×GCM、交互↔三阶残差 分层检验, 使结论对 GCM 选择稳健;
#     事后 emmeans(model="univariate") 走正确误差层; SSP 有序 → 线性趋势检验。
# ============================================================================

rm(list = ls())
gc()

# ============================================================================
# TAP Step 0: 环境准备
# ============================================================================

library(terra)        # 栅格处理与面积计算
library(sf)           # 矢量边界
library(tidyverse)    # 数据处理
library(car)          # Levene检验, Type III SS
library(emmeans)      # 边际均值与事后检验
library(multcomp)     # 多重比较 (emmeans cld 依赖)
library(ggplot2)      # 可视化
library(RColorBrewer) # 配色
library(knitr)        # 表格输出
library(DescTools)    # Kendall W 一致性系数
library(afex)         # 重复测量ANOVA (aov_ez) — Tier1 主推断, GCM作随机重复
library(patchwork)    # 多panel拼图 (Fig5)

set.seed(2026)

# 完全平衡设计 (每单元 n=1), Type I = Type II = Type III SS, 方差分解唯一确定;
# 仍设 sum-to-zero 对比, 使 emmeans 边际均值解释为真正的边际效应。
options(contrasts = c("contr.sum", "contr.poly"))

# ============================================================================
# TAP Step 1: 路径与参数设置
# ============================================================================

work_dir    <- "K:/周建maxent"
output_dir  <- file.path(work_dir, "MaxEnt_Fixed_Results")
res_dir     <- file.path(output_dir, "gcm_replicate_analysis")
fig_dir     <- file.path(res_dir, "figures")
tbl_dir     <- file.path(res_dir, "tables")

dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tbl_dir, showWarnings = FALSE, recursive = TRUE)

# 物种列表
species_map <- c(
  "偃松"       = "yangsong",
  "杜香"       = "duxiang",
  "越桔"       = "yuejv",
  "笃斯越桔"   = "dusiyuejv",
  "小黄花菜"   = "xiaohuanghuacai",
  "短瓣金莲花" = "duanbaijinlianhua"
)

# 物种拉丁名 (英文投稿用; 按 species_code 索引)
# 请务必核对确认 (尤其 杜香 / 短瓣金莲花 存在同物异名), 直接改这里即可
species_latin <- c(
  yangsong          = "Pinus pumila",
  duxiang           = "Rhododendron tomentosum",   # syn. Ledum palustre
  yuejv             = "Vaccinium vitis-idaea",
  dusiyuejv         = "Vaccinium uliginosum",
  xiaohuanghuacai   = "Hemerocallis minor",
  duanbaijinlianhua = "Trollius ledebourii"
)

# 图轴用属名缩写 (P. pumila / V. vitis-idaea ...); 全称 species_latin 留作图注/表/首次出现
species_latin_abbr <- sub("^([A-Z])[a-z]+ ", "\\1. ", species_latin)

# 未来情景
scenarios     <- c("ssp126", "ssp245", "ssp585")
scenario_label <- c(ssp126 = "SSP1-2.6", ssp245 = "SSP2-4.5", ssp585 = "SSP5-8.5")

# GCM (主推断中作随机重复单元; 方差分解中作不确定性来源)
gcm_list <- c("BCC-CSM2-MR", "MIROC6", "IPSL-CM6A-LR")

# 研究区边界 (用于面积计算的参考)
boundary_file <- file.path(work_dir, "1矢量边界/大兴安岭范围准确.shp")
boundary <- st_read(boundary_file, quiet = TRUE)
boundary_wgs <- st_transform(boundary, 4326)

# 模型表现/质量汇总 (来自 MaxEnt_Fixed_Results 建模结果)
# species 列为中文名, 与 species_map 的键一致, 末尾用于补充模型质量信息
perf_file <- file.path(output_dir, "model_performance_summary.csv")
perf <- read.csv(perf_file, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
str(perf)

# ============================================================================
# TAP Step 2: 面积提取 (从60个TSS二值栅格)
# ============================================================================
# 适生区面积用 terra::expanse() 计算 (按像元值分类汇总), 不自建面积函数:
#   - 栅格为 WGS84 经纬度坐标, 像元并非等面积 (本研究区像元面积随纬度相差约 14%)
#   - expanse() 对经纬度栅格自动按椭球计算每个像元的真实地表面积, 无需等面积假设
#   - byValue = TRUE 分别返回值 0/1 的面积; 取 value == 1 即适生区面积 (km²)
#   - 研究区边界外像元先 mask 为 NA, expanse 自动忽略
area_df <- data.frame()

for (sp_cn in names(species_map)) {
  sp_code <- species_map[[sp_cn]]
  sp_dir  <- file.path(output_dir, sp_code)

  cat(sprintf("  [%s] %s\n", sp_code, sp_cn))

  # --- 当前适生面积 ---
  cur_file <- file.path(sp_dir, "binary_current_tss.tif")
  if (file.exists(cur_file)) {
    r_cur    <- mask(crop(rast(cur_file), boundary_wgs), boundary_wgs)
    ex_cur   <- terra::expanse(r_cur, unit = "km", byValue = TRUE)
    cur_area <- sum(ex_cur$area[ex_cur$value == 1])
    cat(sprintf("    当前面积: %.2f km²\n", cur_area))
  } else {
    warning(sprintf("当前文件缺失: %s", cur_file))
    next
  }

  # --- 未来适生面积 (3 SSP × 3 GCM) ---
  for (ssp in scenarios) {
    for (gcm in gcm_list) {
      fut_file <- file.path(sp_dir, sprintf("binary_%s_%s_tss.tif", gcm, ssp))
      if (!file.exists(fut_file)) {
        warning(sprintf("未来文件缺失: %s", fut_file))
        next
      }

      r_fut    <- mask(crop(rast(fut_file), boundary_wgs), boundary_wgs)
      ex_fut   <- terra::expanse(r_fut, unit = "km", byValue = TRUE)
      fut_area <- sum(ex_fut$area[ex_fut$value == 1])

      # --- 适生区地理动态分解 (当前 vs 未来 二值图逐像元叠加) ---
      #   收缩 unfilling = 当前有、未来无;  稳定 filling = 两期都有;
      #   新增 overfilling = 当前无、未来有。
      #   单栅格编码 cur*10 + fut, 一次 expanse 取出: 10=收缩 11=稳定 1=新增 (0=两期都无)
      #   解读: 稳定=无扩散假设下的避难所; 稳定+新增=完全扩散下的潜在分布。
      if (!terra::compareGeom(r_cur, r_fut, stopOnError = FALSE)) {
        r_fut <- terra::resample(r_fut, r_cur, method = "near")  # 网格不一致时对齐(类别用近邻)
      }
      chg     <- r_cur * 10 + r_fut
      ex_chg  <- terra::expanse(chg, unit = "km", byValue = TRUE)
      pick    <- function(v) { s <- ex_chg$area[ex_chg$value == v]; if (length(s)) sum(s) else 0 }
      lost_km2   <- pick(10)   # 收缩 unfilling
      stable_km2 <- pick(11)   # 稳定 filling
      gained_km2 <- pick(1)    # 新增 overfilling

      # 计算变化 (净变化 = 新增 - 收缩)
      delta_km2    <- fut_area - cur_area
      pct_change   <- (fut_area - cur_area) / cur_area * 100

      area_df <- rbind(area_df, data.frame(
        species_cn  = sp_cn,
        species_code = sp_code,
        scenario    = ssp,
        scenario_label = scenario_label[[ssp]],
        gcm         = gcm,
        current_area_km2 = cur_area,
        future_area_km2  = fut_area,
        delta_km2   = delta_km2,
        pct_change  = pct_change,
        lost_km2    = lost_km2,                       # 收缩 unfilling (km²)
        stable_km2  = stable_km2,                     # 稳定 filling   (km²)
        gained_km2  = gained_km2,                     # 新增 overfilling (km²)
        unfilling_pct = lost_km2   / cur_area * 100,  # 收缩占当前面积 %
        stability_pct = stable_km2 / cur_area * 100,  # 稳定占当前面积 % (与 unfilling 之和≈100)
        expansion_pct = gained_km2 / cur_area * 100,  # 新增占当前面积 % (pct_change≈expansion-unfilling)
        stringsAsFactors = FALSE
      ))
    }
  }
}

cat(sprintf("\n  共提取 %d 条记录\n\n", nrow(area_df)))

# 保存原始数据
write_csv(area_df, file.path(tbl_dir, "area_raw_data.csv"))

# 物种顺序 = 当前适生面积从小到大 (数据驱动, 全图统一)
area_rank     <- area_df %>%
  distinct(species_code, species_cn, current_area_km2) %>%
  arrange(current_area_km2)

species_order <- area_rank$species_code   # 物种代码, 面积升序
cn_order      <- area_rank$species_cn      # 对应中文名, 面积升序

# 将三因子转为 factor (aov_ez 需要; scenario 按强迫强度排序以便趋势检验)
area_df$species_cn <- factor(area_df$species_cn, levels = cn_order)  # 按当前面积从小到大
area_df$scenario   <- factor(area_df$scenario, levels = c("ssp126", "ssp245", "ssp585"))
area_df$gcm        <- factor(area_df$gcm)

# ============================================================================
# 适生区地理动态: 收缩(unfilling) / 稳定(filling) / 新增(overfilling)
# ============================================================================
# 净变化(pct_change)会掩盖周转; 此处给出增减分量, 区分"稳定避难所"与"需扩散的新增区"。

range_dyn <- area_df %>%
  group_by(species_cn, scenario_label) %>%
  summarise(
    current_km2  = mean(current_area_km2),
    lost_km2     = mean(lost_km2),      # 收缩 unfilling
    stable_km2   = mean(stable_km2),    # 稳定 filling
    gained_km2   = mean(gained_km2),    # 新增 overfilling
    unfilling_pct = mean(unfilling_pct),
    stability_pct = mean(stability_pct),
    expansion_pct = mean(expansion_pct),
    .groups = "drop"
  ) %>%
  arrange(species_cn, scenario_label)

print(knitr::kable(range_dyn, digits = 1,
  caption = "适生区地理动态 (GCM均值): 收缩/稳定/新增 (km² 及占当前面积%)"))

write.csv(range_dyn, file.path(tbl_dir, "range_dynamics_loss_stable_gain.csv"), row.names = FALSE)

# ============================================================================
# TAP Step 3: 投影不确定性方差分解 (Tier 2 — 全因子固定效应ANOVA)
# ============================================================================
# 注: 此处只做"GCM不确定性贡献多少"的方差分解(η²); 物种/情景的假设检验在 Step 5
#     用 GCM 作随机重复的重复测量ANOVA完成 (误差项更诚实)。

cat("\n>>> [Step 3] 投影不确定性方差分解 (Tier 2)\n\n")

# 模型: 物种×SSP×GCM 三因子, 主效应 + 全部二阶交互;
#   单元内 n=1, 三阶交互(species:scenario:gcm)自动成为残差/误差项;
#   平衡设计下 SS 唯一可加, 各项 SS 占总 SS 比例即该来源的方差贡献(η²)。
mod1 <- aov(pct_change ~ (species_cn + scenario + gcm)^2, data = area_df)
# mod1 <- aov(pct_change ~ species_cn + scenario + gcm + species_cn:scenario + species_cn:gcm + scenario:gcm, data = area_df)

cat("───────────────────────────────────────────────────────\n")
cat("  Model: pct_change ~ (Species + SSP + GCM)^2\n")
cat("  (三阶交互作为残差; SS 占比 = 不确定性方差分解)\n")
cat("───────────────────────────────────────────────────────\n\n")

aov_tab <- summary(mod1)[[1]]
print(aov_tab)

# 不确定性方差分解表: 各来源 SS 占总变异比例 (即 η²)
var_part <- data.frame(
  source       = trimws(rownames(aov_tab)),
  df           = aov_tab[["Df"]],
  SS           = aov_tab[["Sum Sq"]],
  F_value      = aov_tab[["F value"]],
  p_value      = aov_tab[["Pr(>F)"]],
  stringsAsFactors = FALSE
)
var_part$pct_variance <- var_part$SS / sum(var_part$SS) * 100

cat("\n不确定性方差分解 (各来源占总变异 %):\n")
print(knitr::kable(var_part, digits = c(0, 0, 1, 2, 4, 2),
  caption = "投影面积变化率的不确定性方差分解 (pct_variance = η²)"))

# GCM 总贡献 = 主效应 + 含 GCM 的全部交互
gcm_contrib  <- sum(var_part$pct_variance[grepl("gcm", var_part$source)])
spp_contrib  <- var_part$pct_variance[var_part$source == "species_cn"]
resid_contrib <- var_part$pct_variance[var_part$source == "Residuals"]
cat(sprintf("\n  物种贡献: %.1f%% | GCM相关贡献(主+交互): %.1f%% | 残差: %.1f%%\n",
            spp_contrib, gcm_contrib, resid_contrib))

write.csv(var_part, file.path(tbl_dir, "uncertainty_partition.csv"), row.names = FALSE)

# 交互项 p 值 (供 Step 5 简单效应判断)
interaction_p <- var_part$p_value[var_part$source == "species_cn:scenario"]

# 保存方差分解结果
sink(file.path(tbl_dir, "variance_partition_results.txt"))
cat("==============================================\n")
cat("  Uncertainty Variance Partitioning\n")
cat("  Model: pct_change ~ (Species + SSP + GCM)^2\n")
cat("  (3-way interaction = residual; balanced design, n=1/cell)\n")
cat("==============================================\n\n")
cat("--- ANOVA table ---\n")
print(aov_tab)
cat("\n--- Variance partition (% of total SS = eta^2) ---\n")
print(var_part, row.names = FALSE)
cat(sprintf("\nSpecies: %.1f%% | GCM-related: %.1f%% | Residual: %.1f%%\n",
            spp_contrib, gcm_contrib, resid_contrib))
sink()

# ============================================================================
# TAP Step 4: 诊断检验
# ============================================================================

cat("\n>>> [Step 4] 模型诊断\n\n")

# ANOVA 残差诊断图 (2×2): 残差-拟合、QQ、尺度-位置、残差-杠杆
png(file.path(fig_dir, "FigS1_ANOVA_Diagnostics.png"),
    width = 2000, height = 2000, res = 300)
par(mfrow = c(2, 2))
plot(mod1)
dev.off()

# 4.1 残差正态性检验
shapiro_res <- shapiro.test(residuals(mod1))
cat(sprintf("  Shapiro-Wilk 正态性检验: W = %.4f, p = %.4f\n",
            shapiro_res$statistic, shapiro_res$p.value))
if (shapiro_res$p.value > 0.05) {
  cat("  → 残差服从正态分布 (p > 0.05)\n\n")
} else {
  cat("  → 残差不服从正态分布 (p < 0.05)，将使用非参数方法验证\n\n")
}

# 4.2 方差齐性检验 (Levene)
levene_res <- leveneTest(pct_change ~ species_cn * scenario, data = area_df)
cat("  Levene 方差齐性检验:\n")
print(levene_res)
if (levene_res$`Pr(>F)`[1] > 0.05) {
  cat("  → 方差齐性满足 (p > 0.05)\n\n")
} else {
  cat("  → 方差不齐 (p < 0.05)，注意解释时需谨慎\n\n")
}

# 4.3 检查异常值 (Cook 距离, 阈值 4/n)
cooks_d <- cooks.distance(mod1)
influential <- which(cooks_d > 4 / length(cooks_d))
if (length(influential) > 0) {
  cat(sprintf("  潜在异常观测 (Cook's D > 4/n): %d 条\n", length(influential)))
  print(area_df[influential, c("species_cn", "scenario", "gcm", "pct_change")])
} else {
  cat("  未检测到异常观测点\n")
}

# ============================================================================
# TAP Step 5: 主推断 — 物种/情景差异 (GCM 作随机重复, 重复测量ANOVA)
# ============================================================================
#
# 与 Step 3 的分工:
#   - Step 3 (mod1, 全固定) = Tier2 不确定性"方差分解": 量化各来源贡献多少变异(η²);
#   - Step 5 (m_main, aov_ez) = Tier1 主推断: 检验物种差异/情景差异, 把 GCM 当随机
#     重复单元, 使结论对"换一个GCM"是否稳健得到诚实评估。
#   误差自动分层: 物种↔物种×GCM, 情景↔情景×GCM, 物种×情景↔三阶残差。
#   事后用 emmeans(model = "univariate") 走每个效应的正确误差层 (而非最小三阶残差);
#   默认的多元(multivariate)解法在仅3个GCM下自由度=2、过度保守, 故显式指定 univariate。

cat("\n>>> [Step 5] 主推断: 物种/情景差异 (GCM作随机重复)\n\n")

# --- 5.0 重复测量ANOVA总表 (主结果) ---
m_main <- aov_ez(id = "gcm", dv = "pct_change",
                 within = c("species_cn", "scenario"), data = area_df,
                 include_aov = TRUE)
main_tab <- m_main$anova_table
cat("───────────────────────────────────────────────────────\n")
cat("  主推断: pct_change ~ species_cn * scenario + Error(gcm/...)\n")
cat("  (物种↔物种×GCM | 情景↔情景×GCM | 交互↔三阶; GG球形校正)\n")
cat("───────────────────────────────────────────────────────\n")
print(main_tab)

sink(file.path(tbl_dir, "main_anova_gcm_random.txt"))
cat("Main inference: repeated-measures ANOVA, GCM as random replicate\n")
cat("Model: pct_change ~ species_cn * scenario + Error(gcm/(species_cn*scenario))\n")
cat("Each effect tested against its own (effect x GCM) error stratum.\n")
cat("==============================================================\n\n")
print(main_tab)
sink()

# --- 5.1 物种主效应: 边际均值 + Tukey + CLD (误差 = 物种×GCM) ---
emm_species <- emmeans(m_main, ~ species_cn, model = "univariate")
cat("\n物种边际均值 (EMM, 误差=物种×GCM):\n")
print(emm_species)

tukey_species <- contrast(emm_species, method = "tukey", adjust = "tukey")
cat("\n物种间Tukey HSD多重比较:\n")
print(tukey_species)

cld_species <- cld(emm_species, Letters = letters, adjust = "tukey")
cat("\n物种紧凑字母显示 (CLD):\n")
print(cld_species)

# --- 5.2 SSP情景: 有序因子, 以单调趋势为主、两两为辅 (误差 = 情景×GCM) ---
emm_ssp <- emmeans(m_main, ~ scenario, model = "univariate")
cat("\nSSP情景边际均值:\n")
print(emm_ssp)

# 有序因子: 随强迫增强的单调趋势 (linear) 比两两Tukey更有功效、更切题
# (poly 假定等间距序; 此处作"是否随SSP序单调变化"的趋势近似)
trend_ssp <- contrast(emm_ssp, "poly")
cat("\nSSP趋势检验 (正交多项式; linear = 随强迫的单调变化):\n")
print(trend_ssp)

tukey_ssp <- contrast(emm_ssp, method = "tukey", adjust = "tukey")
cat("\nSSP间Tukey HSD多重比较 (辅助):\n")
print(tukey_ssp)

cld_ssp <- cld(emm_ssp, Letters = letters, adjust = "tukey")
cat("\nSSP紧凑字母显示:\n")
print(cld_ssp)

# --- 5.3 交互显著时的简单主效应分析 (误差层同走 univariate) ---
if (interaction_p < 0.05) {
  cat("\n--- 物种×情景交互显著，进行简单主效应分析 ---\n\n")

  # 每个SSP下物种的简单效应
  simple_species_ssp <- emmeans(m_main, ~ species_cn | scenario, model = "univariate")
  tukey_simple <- contrast(simple_species_ssp, method = "tukey", adjust = "tukey")

  sink(file.path(tbl_dir, "simple_effects_species_by_ssp.txt"))
  cat("简单主效应: 每个SSP情景下物种间差异 (GCM作随机重复)\n")
  cat("==============================================\n\n")
  print(tukey_simple)
  sink()

  # 每个物种下SSP的简单效应
  simple_ssp_species <- emmeans(m_main, ~ scenario | species_cn, model = "univariate")
  tukey_simple2 <- contrast(simple_ssp_species, method = "tukey", adjust = "tukey")

  sink(file.path(tbl_dir, "simple_effects_ssp_by_species.txt"))
  cat("简单主效应: 每个物种下SSP情景间差异 (GCM作随机重复)\n")
  cat("==============================================\n\n")
  print(tukey_simple2)
  sink()
}

# 保存主推断与事后检验结果
sink(file.path(tbl_dir, "posthoc_results.txt"))
cat("==================================================================\n")
cat("  Main inference & post-hoc (GCM as random replicate; afex::aov_ez)\n")
cat("==================================================================\n\n")
cat("--- Repeated-measures ANOVA (omnibus) ---\n")
print(main_tab)
cat("\n--- Species Main Effect: Tukey HSD (error = species x GCM) ---\n")
print(tukey_species)
cat("\n--- Species CLD ---\n")
print(cld_species)
cat("\n--- SSP Scenario: linear/quadratic trend (ordered factor) ---\n")
print(trend_ssp)
cat("\n--- SSP Scenario: Tukey HSD (auxiliary) ---\n")
print(tukey_ssp)
cat("\n--- SSP CLD ---\n")
print(cld_ssp)
sink()

# ============================================================================
# TAP Step 6: GCM 集合排序一致性 (Kendall W)
# ============================================================================
#
# 方差分解已量化 GCM 的"方差贡献"; 此处用 Kendall W 补充另一面:
# 三个 GCM 对各物种变化率的"排序"是否一致 (集合预测的稳健性)。
# 已删除: (a) Friedman 物种检验 — 与方差分解中物种效应重复, 且残差正态(无需非参);
#         (b) 当前vs未来 Wilcoxon — n=3 理论上无法显著(最小p≈0.25), 结果误导。

cat("\n>>> [Step 6] GCM 集合排序一致性 (Kendall W)\n\n")

kendall_results <- data.frame()
for (ssp in scenarios) {
  ssp_wide <- area_df %>%
    filter(scenario == ssp) %>%
    dplyr::select(gcm, species_cn, pct_change) %>%
    pivot_wider(names_from = species_cn, values_from = pct_change)
  mat <- as.matrix(ssp_wide[, -1])   # 行=GCM, 列=物种

  if (nrow(mat) >= 2) {
    # DescTools::KendallW 要求 行=被评对象(物种)、列=评分者(GCM), 故转置 t(mat)
    # W ∈ [0,1], 越接近 1 表示各 GCM 对物种变化率的排序越一致
    kres <- DescTools::KendallW(t(mat), correct = FALSE, test = TRUE)
    cat(sprintf("  Kendall W [%s]: W = %.3f (GCM间排序一致性, Friedman p = %.4f)\n",
                scenario_label[[ssp]], unname(kres$estimate), kres$p.value))
    kendall_results <- rbind(kendall_results, data.frame(
      scenario = scenario_label[[ssp]],
      kendall_W = unname(kres$estimate),
      p_value   = kres$p.value,
      stringsAsFactors = FALSE
    ))
  }
}
write.csv(kendall_results, file.path(tbl_dir, "kendall_w_results.csv"), row.names = FALSE)

# ============================================================================
# TAP Step 7: 出版级可视化 (正文 Fig1-5 + 附件 FigS1-S4; 不确定性分析归附件)
# ============================================================================

cat("\n>>> [Step 7] 生成出版级图表...\n\n")

# 配色
species_colors <- brewer.pal(6, "Set2")
names(species_colors) <- names(species_map)

ssp_colors <- c("SSP1-2.6" = "#2c7bb6", "SSP2-4.5" = "#fdae61", "SSP5-8.5" = "#d7191c")

# 公用主题
theme_pub <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3, color = "gray90"),
    strip.background = element_rect(fill = "gray95"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

# -----------------------------------------------------------------------
# Fig S3 (附件): 不确定性方差分解 (Tier 2)
# -----------------------------------------------------------------------
cat("  [附件] FigS3 — 不确定性方差分解\n")

vp_plot <- var_part %>%
  mutate(
    source_label = dplyr::recode(source,
      "species_cn"          = "物种",
      "scenario"            = "SSP情景",
      "gcm"                 = "GCM",
      "species_cn:scenario" = "物种 × SSP",
      "species_cn:gcm"      = "物种 × GCM",
      "scenario:gcm"        = "SSP × GCM",
      "Residuals"           = "残差(三阶交互)"),
    grp = case_when(
      grepl("gcm", source)  ~ "GCM相关",
      source == "Residuals" ~ "残差",
      TRUE                  ~ "物种/SSP")
  ) %>%
  arrange(pct_variance)
vp_plot$source_label <- factor(vp_plot$source_label, levels = vp_plot$source_label)

p_vp <- ggplot(vp_plot, aes(x = source_label, y = pct_variance, fill = grp)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", pct_variance)), hjust = -0.15, size = 3.5) +
  coord_flip() +
  scale_fill_manual(values = c("GCM相关" = "#d7191c", "物种/SSP" = "#2c7bb6", "残差" = "gray60"),
                    name = "来源类别") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Fig.S3 投影面积变化率的不确定性方差分解",
    subtitle = "各来源 SS 占总变异比例 (η²) | 红=GCM相关, 蓝=物种/SSP",
    x = "", y = "方差贡献 (% of total)"
  ) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave(file.path(fig_dir, "FigS3_Uncertainty_Partition.png"), p_vp,
       width = 9, height = 5, dpi = 300)

# -----------------------------------------------------------------------
# Fig 1: 适生区面积变化率箱线图 (核心图)
# -----------------------------------------------------------------------
cat("  [1/5] Fig1 — 面积变化率箱线图\n")

area_df$scenario_label <- factor(area_df$scenario_label,
                                  levels = c("SSP1-2.6", "SSP2-4.5", "SSP5-8.5"))

p1 <- ggplot(area_df, aes(x = species_cn, y = pct_change, fill = species_cn)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(shape = gcm), width = 0.15, size = 2.5, alpha = 0.8) +
  facet_wrap(~ scenario_label, ncol = 3) +
  scale_fill_manual(values = species_colors, guide = "none") +
  scale_shape_manual(values = c(16, 17, 15), name = "GCM") +
  labs(
    title = "Fig.1 适生区面积变化率 (TSS阈值, n=3 GCM)",
    subtitle = expression("面积变化率 % = (Future - Current) / Current × 100%"),
    x = "物种", y = "面积变化率 (%)"
  ) +
  theme_pub

ggsave(file.path(fig_dir, "Fig1_AreaChange_Boxplot.png"), p1,
       width = 13, height = 6, dpi = 300)

# -----------------------------------------------------------------------
# Fig 5 (核心): Range dynamics — unfilling / filling / overfilling, mean ± SE
#   3 panels (one per component), 误差棒 = ±1 SE 跨 3 个 GCM; 物种用拉丁名; 英文(投稿用)
# -----------------------------------------------------------------------
cat("  [核心] Fig5 — range dynamics 3-panel ±SE (English)\n")

# 跨 GCM 求每个分量(占当前面积%)的均值与标准误 (n = 3 GCMs)
dyn_se <- area_df %>%
  mutate(sp_latin = factor(species_latin_abbr[species_code], levels = unname(species_latin_abbr[species_order])),
         scen     = factor(scenario_label, levels = c("SSP1-2.6", "SSP2-4.5", "SSP5-8.5"))) %>%
  group_by(sp_latin, scen) %>%
  summarise(across(c(unfilling_pct, stability_pct, expansion_pct),
                   list(m = ~mean(.x), se = ~sd(.x) / sqrt(n())), .names = "{.col}_{.fn}"),
            .groups = "drop")
write.csv(dyn_se, file.path(tbl_dir, "range_dynamics_mean_se.csv"), row.names = FALSE)

scen_cols <- c("SSP1-2.6" = "#2c7bb6", "SSP2-4.5" = "#fdae61", "SSP5-8.5" = "#d7191c")

# 一个分量一个 panel: 物种(x) × 情景(分组柱) + ±SE 误差棒
# show_x: 仅最后一行(c)保留物种轴标签, 上方 a/b 隐去以省空间
make_se_panel <- function(mcol, secol, ttl, show_x = TRUE) {
  p <- ggplot(dyn_se, aes(x = sp_latin, y = .data[[mcol]], fill = scen)) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    geom_errorbar(aes(ymin = .data[[mcol]] - .data[[secol]],
                      ymax = .data[[mcol]] + .data[[secol]]),
                  position = position_dodge(0.8), width = 0.25, linewidth = 0.4) +
    scale_fill_manual(values = scen_cols, name = "Scenario") +
    labs(title = ttl, x = NULL, y = "% of current area") +
    theme_pub +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "italic"))
  if (!show_x) p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p
}

p5a <- make_se_panel("unfilling_pct_m", "unfilling_pct_se", "(a) Unfilling (range loss)", show_x = FALSE)
p5b <- make_se_panel("stability_pct_m", "stability_pct_se", "(b) Stability (retained)",   show_x = FALSE)
p5c <- make_se_panel("expansion_pct_m", "expansion_pct_se", "(c) Overfilling (range gain)", show_x = TRUE)

p_dyn <- (p5a / p5b / p5c) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Fig.5 Suitable-area range dynamics under climate change",
    subtitle = "Mean ± SE across 3 GCMs | % of current suitable area"
  ) & theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "Fig5_RangeDynamics_LossStableGain.png"), p_dyn,
       width = 8, height = 11, dpi = 300)

# -----------------------------------------------------------------------
# Fig 6 (Part1: 当前适生分布图 + 面积柱形图) 已拆为独立脚本:
#   004b.Fig6_CurrentSuitability.R  (读 area_raw_data.csv + binary_current_tss.tif,
#   无需重跑本脚本即可快速调整出图)
# -----------------------------------------------------------------------

# -----------------------------------------------------------------------
# Fig 2: 绝对适生面积对比 (Current vs Future × 3 SSP)
# -----------------------------------------------------------------------
cat("  [2/5] Fig2 — 绝对适生面积对比图\n")

# 长格式: current + 3 future SSP均值
area_long <- area_df %>%
  group_by(species_cn, scenario_label) %>%
  summarise(
    current = mean(current_area_km2),
    future  = mean(future_area_km2),
    sd_future = sd(future_area_km2),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(current, future), names_to = "time", values_to = "area_km2") %>%
  mutate(
    time = factor(time, levels = c("current", "future"), labels = c("当前", "未来")),
    sd = ifelse(time == "未来", sd_future, 0)
  )

p2 <- ggplot(area_long, aes(x = species_cn, y = area_km2, fill = time)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = area_km2 - sd, ymax = area_km2 + sd),
    position = position_dodge(width = 0.8), width = 0.2,
    data = area_long %>% filter(time == "未来")
  ) +
  geom_text(
    aes(label = sprintf("%.0f", area_km2)),
    position = position_dodge(width = 0.8),
    vjust = -0.5, size = 2.8
  ) +
  facet_wrap(~ scenario_label, ncol = 3) +
  scale_fill_manual(values = c("当前" = "#66c2a5", "未来" = "#fc8d62")) +
  labs(
    title = "Fig.2 当前 vs 未来适生区面积 (km²)",
    subtitle = "误差棒: ±1 SD (n=3 GCM)",
    x = "物种", y = expression("适生区面积 (km"^2*")"), fill = "时期"
  ) +
  theme_pub

ggsave(file.path(fig_dir, "Fig2_AbsoluteArea_Barplot.png"), p2,
       width = 14, height = 6, dpi = 300)

# -----------------------------------------------------------------------
# Fig 3: 物种 × SSP 面积变化率热图
# -----------------------------------------------------------------------
cat("  [3/5] Fig3 — 面积变化率热图\n")

heat_data <- area_df %>%
  group_by(species_cn, scenario_label) %>%
  summarise(
    pct_mean = mean(pct_change),
    pct_sd   = sd(pct_change),
    .groups  = "drop"
  ) %>%
  mutate(
    label = sprintf("%.1f%%\n(±%.1f)", pct_mean, pct_sd)
  )

p3 <- ggplot(heat_data, aes(x = scenario_label, y = species_cn, fill = pct_mean)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = label), size = 3.2, fontface = "bold") +
  scale_fill_gradient2(
    low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
    midpoint = 0, name = "变化率 (%)"
  ) +
  labs(
    title = "Fig.3 适生区面积变化率热图",
    subtitle = "值 = 均值(GCM n=3) ± SD | 蓝色=收缩, 红色=扩张",
    x = "SSP情景", y = "物种"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.text = element_text(size = 11)
  )

ggsave(file.path(fig_dir, "Fig3_AreaChange_Heatmap.png"), p3,
       width = 9, height = 6, dpi = 300)

# -----------------------------------------------------------------------
# Fig 4: 物种间多重比较 (Tukey HSD CLD)
# -----------------------------------------------------------------------
cat("  [4/5] Fig4 — 物种间Tukey HSD多重比较\n")

cld_plot_data <- as.data.frame(cld_species)

p4 <- ggplot(cld_plot_data, aes(x = factor(species_cn, levels = cn_order), y = emmean)) +
  geom_point(size = 4, color = "#2c7bb6") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.15,
                linewidth = 0.8, color = "#2c7bb6") +
  geom_text(aes(label = .group), vjust = -1.2, size = 5, fontface = "bold") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "Fig.4 物种边际均值 (Tukey HSD, α = 0.05)",
    subtitle = "GCM作随机重复 (误差=物种×GCM) | 不同字母示差异显著 | 误差线: 95% CI",
    x = "物种", y = "面积变化率边际均值 (%)"
  ) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 11))

ggsave(file.path(fig_dir, "Fig4_Species_TukeyHSD.png"), p4,
       width = 8, height = 5, dpi = 300)

# -----------------------------------------------------------------------
# Fig S4 (附件): GCM变异性展示 (各情景下GCM间一致性)
# -----------------------------------------------------------------------
cat("  [附件] FigS4 — GCM间变异性\n")

p5 <- ggplot(area_df, aes(x = species_cn, y = future_area_km2,
                            fill = gcm, group = gcm)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(aes(yintercept = current_area_km2), linetype = "dashed",
             color = "gray30", linewidth = 0.7) +
  facet_wrap(~ scenario_label, ncol = 3) +
  scale_fill_brewer(palette = "Set1", name = "GCM") +
  labs(
    title = "Fig.S4 各GCM预测的适生区面积对比",
    subtitle = "虚线: 当前适生面积 | 柱: 各GCM未来预测",
    x = "物种", y = expression("适生区面积 (km"^2*")")
  ) +
  theme_pub

ggsave(file.path(fig_dir, "FigS4_GCM_Variability.png"), p5,
       width = 14, height = 7, dpi = 300)

# ============================================================================
# 补充图: 当前各物种面积对比
# ============================================================================

current_areas <- area_df %>%
  group_by(species_cn) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(species_cn = reorder(species_cn, current_area_km2))

p_current <- ggplot(current_areas, aes(x = species_cn, y = current_area_km2, fill = species_cn)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.0f km²", current_area_km2)),
            hjust = -0.1, size = 3.5) +
  scale_fill_manual(values = species_colors, guide = "none") +
  coord_flip() +
  labs(
    title = "当前气候条件下各物种适生区面积 (TSS阈值)",
    x = "", y = expression("适生区面积 (km"^2*")")
  ) +
  theme_pub +
  theme(axis.text.y = element_text(size = 11))

ggsave(file.path(fig_dir, "FigS2_CurrentArea_AllSpecies.png"), p_current,
       width = 8, height = 5, dpi = 300)

# ============================================================================
# 结果汇总表格导出
# ============================================================================

# 表1: 综合结果表 (跨情景)
summary_table <- area_df %>%
  group_by(species_cn) %>%
  summarise(
    current_area_km2 = mean(current_area_km2),
    ssp126_mean_pct  = mean(pct_change[scenario == "ssp126"]),
    ssp245_mean_pct  = mean(pct_change[scenario == "ssp245"]),
    ssp585_mean_pct  = mean(pct_change[scenario == "ssp585"]),
    overall_mean_pct = mean(pct_change),
    overall_sd_pct   = sd(pct_change),
    .groups = "drop"
  ) %>%
  arrange(overall_mean_pct)

# 补充模型质量信息 (来自 MaxEnt_Fixed_Results/model_performance_summary.csv)
# 便于在解释面积变化时同时审视各物种模型的可靠性
perf_quality <- perf %>%
  transmute(
    species_cn    = species,
    n_occ         = n_occ,
    best_FC       = best_FC,
    best_RM       = best_RM,
    AUC_cv        = AUC_cv,
    TSS_full      = TSS_full,
    Boyce_full    = Boyce_full,
    CBI_cv        = CBI_cv,
    tss_threshold = tss_threshold
  )
write.csv(perf_quality, file.path(tbl_dir, "model_performance_used.csv"), row.names = FALSE)

summary_table <- summary_table %>%
  left_join(
    perf_quality %>%
      dplyr::select(species_cn, n_occ, AUC_cv, TSS_full, Boyce_full, tss_threshold),
    by = "species_cn"
  )

cat("\n综合结果表 (面积变化率 + 模型质量):\n")
print(knitr::kable(summary_table, digits = 2,
  caption = "各物种适生区面积变化率与模型表现汇总"))

write.csv(summary_table, file.path(tbl_dir, "summary_results.csv"), row.names = FALSE)

cat("\n===========================================================\n")
cat("  分析完成!\n")
cat("  结果保存路径: ", res_dir, "\n")
cat("  图表保存路径: ", fig_dir, "\n")
cat("  表格保存路径: ", tbl_dir, "\n")
cat("===========================================================\n\n")

# 输出Session信息
cat("R Session Info:\n")
sessionInfo()
