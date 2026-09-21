# ============================================================================
# 005b.ANOVA_Figures.R
# 重复测量ANOVA 事后检验 + 两大图
#
# Fig 3: 总体变化 (pct_change) — 物种边际均值 + SSP趋势
# Fig 4: 面积分量 (lost / stable / gained) — 物种×SSP 交互
# ============================================================================

rm(list = ls()); gc()
suppressMessages({
  library(tidyverse)
  library(afex); library(emmeans); library(multcomp)
  library(patchwork); library(scales); library(knitr)
})

options(contrasts = c("contr.sum", "contr.poly"))

# --- 路径 ---
work_dir <- "K:/周建maxent"
res_dir  <- file.path(work_dir, "MaxEnt_Fixed_Results", "gcm_replicate_analysis")
tbl_dir  <- file.path(res_dir, "tables")
fig_dir  <- file.path(res_dir, "figures")
out_dir  <- file.path(res_dir, "anova_results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- 颜色 ---
COL_LOST   <- "#D32F2F"
COL_STABLE <- "#2E7D32"
COL_GAINED <- "#1976D2"
COL_SSP    <- c("SSP1-2.6" = "#1B7837", "SSP2-4.5" = "#E6AB02",
                "SSP5-8.5" = "#D73027")

# --- 数据 ---
area_df <- read.csv(file.path(tbl_dir, "area_raw_data.csv"), fileEncoding = "UTF-8")

area_rank <- area_df %>%
  distinct(species_code, species_cn, current_area_km2) %>%
  arrange(current_area_km2)

area_df$species_cn <- factor(area_df$species_cn, levels = area_rank$species_cn)
area_df$scenario   <- factor(area_df$scenario, levels = c("ssp126", "ssp245", "ssp585"))
area_df$gcm        <- factor(area_df$gcm)
area_df$scenario_label <- factor(area_df$scenario_label,
  levels = c("SSP1-2.6", "SSP2-4.5", "SSP5-8.5"))

# ============================================================================
# 1. ANOVA: pct_change (总体变化)
# ============================================================================
cat(">>> ANOVA: pct_change\n")
m_pct <- aov_ez(id = "gcm", dv = "pct_change",
                within = c("species_cn", "scenario"),
                data = area_df, include_aov = TRUE)
print(m_pct$anova_table)

emm_spp_pct <- emmeans(m_pct, ~ species_cn, model = "univariate")
cld_spp_pct <- cld(emm_spp_pct, Letters = letters, adjust = "tukey")

emm_ssp_pct <- emmeans(m_pct, ~ scenario, model = "univariate")
trend_pct    <- contrast(emm_ssp_pct, "poly")

# 交互简单效应
emm_interact <- emmeans(m_pct, ~ species_cn | scenario, model = "univariate")
cld_interact <- cld(emm_interact, Letters = letters, adjust = "tukey")

# ============================================================================
# 2. ANOVA: lost / stable / gained
# ============================================================================
m_lost   <- aov_ez(id = "gcm", dv = "lost_km2",   within = c("species_cn", "scenario"), data = area_df, include_aov = TRUE)
m_stable <- aov_ez(id = "gcm", dv = "stable_km2", within = c("species_cn", "scenario"), data = area_df, include_aov = TRUE)
m_gained <- aov_ez(id = "gcm", dv = "gained_km2", within = c("species_cn", "scenario"), data = area_df, include_aov = TRUE)

emm_lost_spp   <- cld(emmeans(m_lost,   ~ species_cn, model = "univariate"), Letters = letters)
emm_stable_spp <- cld(emmeans(m_stable, ~ species_cn, model = "univariate"), Letters = letters)
emm_gained_spp <- cld(emmeans(m_gained, ~ species_cn, model = "univariate"), Letters = letters)

# 交互 (lost/stable 有相同结构)
emm_lost_int   <- emmeans(m_lost,   ~ species_cn | scenario, model = "univariate")
emm_stable_int <- emmeans(m_stable, ~ species_cn | scenario, model = "univariate")
emm_gained_int <- emmeans(m_gained, ~ species_cn | scenario, model = "univariate")

# ============================================================================
# 3. 保存事后检验文本
# ============================================================================
sink(file.path(out_dir, "posthoc_complete.txt"))
cat("==========================================================\n")
cat("  Complete Post-hoc: aov_ez (GCM random replicate)\n")
cat("==========================================================\n\n")

cat("--- 1. pct_change: Species main effect (CLD) ---\n")
print(as.data.frame(cld_spp_pct))

cat("\n--- 2. pct_change: SSP linear trend ---\n")
print(trend_pct)

cat("\n--- 3. pct_change: Species x SSP simple effects (CLD) ---\n")
print(as.data.frame(cld_interact))

cat("\n--- 4. Area components: Species CLD ---\n")
cat("Lost:\n"); print(as.data.frame(emm_lost_spp))
cat("Stable:\n"); print(as.data.frame(emm_stable_spp))
cat("Gained:\n"); print(as.data.frame(emm_gained_spp))
sink()

# ============================================================================
# 准备: 拉丁名缩写 (与 Fig6 一致) + CLD 字母
# ============================================================================
species_latin_map <- c(
  "小黄花菜" = "H. minor",
  "短瓣金莲花" = "T. ledebourii",
  "笃斯越桔" = "V. uliginosum",
  "偃松" = "P. pumila",
  "越桔" = "V. vitis-idaea",
  "杜香" = "R. tomentosum"
)
area_df$latin <- factor(species_latin_map[as.character(area_df$species_cn)],
                        levels = species_latin_map[levels(area_df$species_cn)])

# --- CLD: 物种间 (大写) ---
# --- CLD: 物种间 (大写) ---
cld_spp_pct    <- multcomp::cld(emmeans(m_pct,   ~ species_cn, model = "univariate"), Letters = LETTERS)
cld_spp_lost   <- multcomp::cld(emmeans(m_lost,   ~ species_cn, model = "univariate"), Letters = LETTERS)
cld_spp_stable <- multcomp::cld(emmeans(m_stable, ~ species_cn, model = "univariate"), Letters = LETTERS)
cld_spp_gained <- multcomp::cld(emmeans(m_gained, ~ species_cn, model = "univariate"), Letters = LETTERS)

cld_spp_upper <- data.frame(
  latin = cld_spp_pct$species_cn,  # same order for all
  pct_change   = as.character(cld_spp_pct$.group),
  lost         = as.character(cld_spp_lost$.group),
  stable       = as.character(cld_spp_stable$.group),
  gained       = as.character(cld_spp_gained$.group)
)
cld_spp_upper$latin <- factor(species_latin_map[as.character(cld_spp_upper$latin)],
                              levels = levels(area_df$latin))

# --- CLD: SSP 内物种 (小写) ---
cld_ssp_within <- function(dv_name) {
  ssp_map <- c(ssp126 = "SSP1-2.6", ssp245 = "SSP2-4.5", ssp585 = "SSP5-8.5")
  result <- list()
  for (sp in levels(area_df$species_cn)) {
    sp_data <- subset(area_df, species_cn == sp)
    sp_model <- aov_ez(id = "gcm", dv = dv_name, within = "scenario",
                       data = sp_data, include_aov = TRUE)
    e <- emmeans(sp_model, ~ scenario, model = "univariate")
    d <- as.data.frame(multcomp::cld(e, Letters = letters, adjust = "tukey"))
    d$species_cn <- sp
    result[[sp]] <- d
  }
  d_all <- do.call(rbind, result)
  d_all$latin <- species_latin_map[as.character(d_all$species_cn)]
  d_all$scenario_label <- ssp_map[as.character(d_all$scenario)]
  d_all
}
cld_lower_pct    <- cld_ssp_within("pct_change")
cld_lower_lost   <- cld_ssp_within("lost_km2")
cld_lower_stable <- cld_ssp_within("stable_km2")
cld_lower_gained <- cld_ssp_within("gained_km2")

# ============================================================================
# Fig 3: 总体变化 (pct_change) — 物种×SSP + CLD 字母
# ============================================================================
pct_summary <- area_df %>%
  group_by(latin, scenario_label) %>%
  summarise(
    mean_pct = mean(pct_change),
    se_pct   = sd(pct_change) / sqrt(n()),
    .groups = "drop"
  ) %>%
  left_join(cld_spp_upper %>% dplyr::select(latin, upper = pct_change), by = "latin") %>%
  left_join(cld_lower_pct %>% dplyr::select(latin, scenario_label, lower = .group),
            by = c("latin", "scenario_label"))

fig3 <- ggplot(pct_summary, aes(latin, mean_pct, fill = scenario_label)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean_pct - se_pct, ymax = mean_pct + se_pct),
                position = position_dodge(0.8), width = 0.2, linewidth = 0.3) +
  geom_text(data = pct_summary %>% distinct(latin, upper) %>% mutate(scenario_label = NA_character_),
            aes(x = latin, label = upper,
                y = max(pct_summary$mean_pct + pct_summary$se_pct) + 8),
            size = 4.5, vjust = 0) +
  geom_text(aes(label = lower,
                y = ifelse(mean_pct >= 0, mean_pct + se_pct + 5, mean_pct - se_pct - 5),
                group = scenario_label),
            position = position_dodge(0.8), size = 3.5, vjust = 0) +
  scale_fill_manual(values = COL_SSP, name = NULL, na.translate = FALSE) +
  scale_x_discrete(limits = levels(area_df$latin)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
    scale_y_continuous(breaks = seq(-100, 200, by = 50),
                       limits = c(-100, 200)) +
  labs(x = NULL, y = expression(Delta*" Suitable area (%)")) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(face = "italic", hjust = 0.5, size = 11),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        legend.position = "top",
        legend.text = element_text(size = 11),
        legend.key.size = unit(0.6, "cm"),
        plot.title = element_blank())

ggsave(file.path(fig_dir, "Fig3_OverallChange.png"), fig3,
       width = 9, height = 6, dpi = 300)

# ============================================================================
# Fig 4: 面积分量百分比 — 三个独立子图纵排 + patchwork + 共享图例
# ============================================================================
components <- c("Lost", "Stable", "Gained")
comp_labels <- c("A  Lost", "B  Stable", "C  Gained")
pct_comp_names <- c("unfilling_pct", "stability_pct", "expansion_pct")

fig4_panels <- list()
for (ci in seq_along(components)) {
  comp <- components[ci]
  cname <- pct_comp_names[ci]

  comp_data <- area_df %>%
    dplyr::select(latin, scenario_label, gcm, pct = all_of(cname)) %>%
    group_by(latin, scenario_label) %>%
    summarise(mean_pct = mean(pct), se_pct = sd(pct) / sqrt(n()), .groups = "drop") %>%
    arrange(latin)

  # CLD letters: uppercase (between-species) + lowercase (within-SSP)
  upper_col <- tolower(comp)  # "lost", "stable", "gained"
  comp_upper <- cld_spp_upper %>% dplyr::select(latin, upper = all_of(upper_col))
  if (comp == "Lost") {
    comp_lower <- cld_lower_lost %>% dplyr::select(latin, scenario_label, lower = .group)
  } else if (comp == "Stable") {
    comp_lower <- cld_lower_stable %>% dplyr::select(latin, scenario_label, lower = .group)
  } else {
    comp_lower <- cld_lower_gained %>% dplyr::select(latin, scenario_label, lower = .group)
  }

  comp_data <- comp_data %>%
    mutate(ymin = mean_pct - se_pct, ymax = mean_pct + se_pct) %>%
    left_join(comp_lower, by = c("latin", "scenario_label"))

  upper_dat <- comp_data %>%
    left_join(comp_upper, by = "latin") %>%
    distinct(latin, upper) %>%
    mutate(scenario_label = NA_character_,
           y_pos = max(comp_data$ymax) * 1.15)

  brks <- if (comp == "Gained") seq(0, 200, by = 50) else seq(0, 120, by = 20)

  p <- local({
    brks <- brks
    ggplot(comp_data, aes(latin, mean_pct, fill = scenario_label)) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax),
                  position = position_dodge(0.8), width = 0.2, linewidth = 0.35) +
    geom_text(data = upper_dat,
              aes(x = latin, label = upper, y = y_pos),
              size = 7, vjust = 0.5) +
    geom_text(aes(label = lower, y = ymax + 5, group = scenario_label),
              position = position_dodge(0.8), size = 5, vjust = 0) +
    scale_fill_manual(values = COL_SSP, name = NULL, na.translate = FALSE,
                      guide = if (ci == 1) "legend" else "none") +
    scale_x_discrete(limits = levels(area_df$latin)) +
    scale_y_continuous(breaks = brks, limits = range(brks)) +
    labs(x = NULL, y = "% of current area", title = comp_labels[ci]) +
    theme_classic(base_size = 22) +
    theme(axis.text.x = if (ci == 3)
            element_text(face = "italic", hjust = 0.5, size = 22)
          else element_blank(),
          axis.ticks.x = if (ci == 3) element_line() else element_blank(),
          axis.text.y = element_text(size = 18),
          axis.title.y = element_text(size = 21),
          plot.title = element_text(size = 24, hjust = 0.02),
          legend.position = if (ci == 1) "top" else "none",
          legend.text = element_text(size = 19),
          legend.key.size = unit(1, "cm"))
  })

  fig4_panels[[ci]] <- p
}

fig4 <- wrap_plots(fig4_panels, ncol = 1)

ggsave(file.path(fig_dir, "Fig4_AreaComponents.png"), fig4,
       width = 14, height = 20, dpi = 300)

# ============================================================================
# 汇总 ANOVA 表格 (用于论文)
# ============================================================================
extract_anova <- function(m, dv_name) {
  tab <- m$anova_table
  data.frame(
    DV = dv_name,
    Effect = rownames(tab),
    df_num = tab[["num Df"]],
    df_den = tab[["den Df"]],
    F = round(tab[["F"]], 2),
    p = round(tab[["Pr(>F)"]], 4),
    sig = ifelse(tab[["Pr(>F)"]] < 0.001, "***",
          ifelse(tab[["Pr(>F)"]] < 0.01, "**",
          ifelse(tab[["Pr(>F)"]] < 0.05, "*", "ns"))),
    stringsAsFactors = FALSE
  )
}

full_tab <- rbind(
  extract_anova(m_pct,   "pct_change"),
  extract_anova(m_lost,   "lost_km2"),
  extract_anova(m_stable, "stable_km2"),
  extract_anova(m_gained, "gained_km2")
)
write.csv(full_tab, file.path(out_dir, "anova_table_paper_ready.csv"), row.names = FALSE)

cat("\n\n=== 完成 ===\n")
cat("Fig 3:", file.path(fig_dir, "Fig3_OverallChange.png"), "\n")
cat("Fig 4:", file.path(fig_dir, "Fig4_AreaComponents.png"), "\n")
cat("posthoc:", file.path(out_dir, "posthoc_complete.txt"), "\n")
