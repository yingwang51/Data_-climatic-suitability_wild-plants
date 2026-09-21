# ============================================================================
# 008.DriverResponse_Analysis.R
# Fig. 24 — ComplexHeatmap: permutation importance + embedded PDP curves
# All 10 variables, 6 species, single integrated figure
# ============================================================================

rm(list = ls())
gc()

library(maxnet)
library(terra)
library(dplyr)
library(tidyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)
set.seed(2026)

# ============================================================================
# Part 0: Paths & config
# ============================================================================
work_dir <- Sys.getenv("MAXENT_PROJECT_DIR", unset = "K:/周建maxent")
res_dir  <- file.path(work_dir, "MaxEnt_Fixed_Results")
fix_dir  <- file.path(work_dir, "MaxEnt_Fixed_Results")
ana_dir  <- file.path(fix_dir, "gcm_replicate_analysis")
fig_dir  <- file.path(ana_dir, "figures")
out_dir  <- file.path(ana_dir, "driver_analysis")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Species (abbreviated genus)
species_map <- c(
  yangsong          = "P. pumila",
  duxiang           = "R. tomentosum",
  yuejv             = "V. vitis-idaea",
  dusiyuejv         = "V. uliginosum",
  xiaohuanghuacai   = "H. minor",
  duanbaijinlianhua = "T. ledebourii"
)

# All 10 variable labels
var_labels <- c(
  bio4        = "BIO-4",
  bio3        = "BIO-3",
  bio9        = "BIO-9",
  bio12       = "BIO-12",
  bio14       = "BIO-14",
  bio15       = "BIO-15",
  elevation   = "Elevation",
  pop_density = "Pop. density",
  slope       = "Slope",
  aspect      = "Aspect"
)

# ============================================================================
# Part 1: Load models & env
# ============================================================================
cat("Loading models and environment data...\n")

model_yangsong  <- readRDS(file.path(res_dir, "yangsong",           "final_model_maxnet.rds"))
model_duanbai   <- readRDS(file.path(res_dir, "duanbaijinlianhua", "final_model_maxnet.rds"))
model_duxiang   <- readRDS(file.path(res_dir, "duxiang",           "final_model_maxnet.rds"))
model_yuejv     <- readRDS(file.path(res_dir, "yuejv",             "final_model_maxnet.rds"))
model_dusiyuejv <- readRDS(file.path(res_dir, "dusiyuejv",         "final_model_maxnet.rds"))
model_xiaohuang <- readRDS(file.path(res_dir, "xiaohuanghuacai",   "final_model_maxnet.rds"))

model_list <- list(
  "P. pumila"         = model_yangsong,
  "R. tomentosum"     = model_duxiang,
  "V. vitis-idaea"    = model_yuejv,
  "V. uliginosum"     = model_dusiyuejv,
  "H. minor"          = model_xiaohuang,
  "T. ledebourii"     = model_duanbai
)

env_final <- terra::rast(file.path(fix_dir, "env_final.tif"))
kept_vars <- names(env_final)
cat(sprintf("  %d environmental variables loaded\n", length(kept_vars)))

# ============================================================================
# Part 2: Permutation importance matrix (ALL 10 variables)
# ============================================================================
n_perm <- 100L
perm_summary_file <- file.path(out_dir, "permutation_importance_summary.csv")
perm_runs_file   <- file.path(out_dir, "permutation_importance_100runs.csv")

# Reuse a complete 100-permutation cache; otherwise recompute from fixed models.
cache_ok <- FALSE
if (file.exists(perm_summary_file) && file.exists(perm_runs_file)) {
  cached_summary <- tryCatch(read.csv(perm_summary_file), error = function(e) NULL)
  cached_runs <- tryCatch(read.csv(perm_runs_file), error = function(e) NULL)
  required_summary <- c("species", "variable", "perm_imp", "perm_imp_raw",
                        "perm_imp_sd", "n_perm")
  required_runs <- c("species", "variable", "permutation", "score")
  expected_keys <- as.vector(outer(unname(species_map), kept_vars, paste, sep = "::"))

  if (!is.null(cached_summary) && !is.null(cached_runs) &&
      all(required_summary %in% names(cached_summary)) &&
      all(required_runs %in% names(cached_runs))) {
    summary_keys <- paste(cached_summary$species, cached_summary$variable, sep = "::")
    cache_ok <- nrow(cached_summary) == length(expected_keys) &&
      !anyDuplicated(summary_keys) && setequal(summary_keys, expected_keys) &&
      all(cached_summary$n_perm == n_perm) &&
      all(is.finite(cached_summary$perm_imp)) &&
      nrow(cached_runs) == length(expected_keys) * n_perm
  }
}

if (cache_ok) {
  cat("Loading cached 100-permutation importance results...\n")
  perm_df <- cached_summary
} else {
  cat("Calculating 100-repeat model-agnostic permutation importance...\n")

  set.seed(2026)
  # maxnet has no native ENMeval importance table. Use prediction sensitivity:
  # importance = 1 - Pearson correlation(original prediction, permuted prediction).
  perm_data <- terra::spatSample(
    env_final, size = 10000, method = "random",
    na.rm = TRUE, values = TRUE, xy = FALSE
  ) %>% as.data.frame()

  perm_list <- list()
  perm_runs_list <- list()
  for (sp_code in names(species_map)) {
    sp_name <- species_map[[sp_code]]
    model   <- model_list[[sp_name]]
    pred_0  <- as.numeric(predict(model, perm_data, type = "cloglog"))
    sp_index <- match(sp_code, names(species_map))
    cat(sprintf("  [%d/%d] %s: 10 variables x %d permutations\n",
                sp_index, length(species_map), sp_name, n_perm))

    score_mat <- sapply(seq_along(kept_vars), function(j) {
      var_name <- kept_vars[j]
      vapply(seq_len(n_perm), function(r) {
        set.seed(2026 + sp_index * 100000 + j * 1000 + r)
        permuted <- perm_data
        permuted[[var_name]] <- sample(permuted[[var_name]])
        pred_p <- as.numeric(predict(model, permuted, type = "cloglog"))
        score <- 1 - cor(pred_0, pred_p, use = "complete.obs")
        if (is.finite(score)) max(0, score) else 0
      }, numeric(1))
    })
    colnames(score_mat) <- kept_vars
    scores_raw <- colMeans(score_mat)
    scores_sd  <- apply(score_mat, 2, sd)
    scores_pct <- if (sum(scores_raw) > 0) {
      100 * scores_raw / sum(scores_raw)
    } else {
      scores_raw
    }

    perm_list[[sp_code]] <- data.frame(
      species = sp_name, variable = kept_vars,
      perm_imp = as.numeric(scores_pct),
      perm_imp_raw = as.numeric(scores_raw),
      perm_imp_sd = as.numeric(scores_sd),
      n_perm = n_perm,
      stringsAsFactors = FALSE
    )
    perm_runs_list[[sp_code]] <- data.frame(
      species = sp_name,
      variable = rep(kept_vars, each = n_perm),
      permutation = rep(seq_len(n_perm), times = length(kept_vars)),
      score = as.numeric(score_mat),
      stringsAsFactors = FALSE
    )
  }

  perm_df <- bind_rows(perm_list)
  perm_runs_df <- bind_rows(perm_runs_list)
  write.csv(perm_runs_df, perm_runs_file, row.names = FALSE)
  write.csv(perm_df, perm_summary_file, row.names = FALSE)
  cat("  Saved raw runs and summary cache.\n")
}

# Species order: by current area (descending)
area_data <- data.frame(
  species = c("R. tomentosum", "V. vitis-idaea", "P. pumila",
              "V. uliginosum", "T. ledebourii", "H. minor"),
  area_km2 = c(32938, 32709, 29085, 28849, 21479, 20721))
sp_order <- area_data %>% arrange(desc(area_km2)) %>% pull(species)

# Variable order: BIO3-BIO19 numerically, followed by non-BIO predictors
bio_order <- paste0("bio", 3:19)
var_order <- c(
  intersect(bio_order, kept_vars),
  setdiff(kept_vars, bio_order)
)

cat(sprintf("  Species (%d): %s\n", length(sp_order), paste(sp_order, collapse=", ")))
cat(sprintf("  Variables (%d): %s\n", length(var_order), paste(var_order, collapse=", ")))

# Build matrix
perm_mat <- perm_df %>%
  dplyr::select(species, variable, perm_imp) %>%
  pivot_wider(names_from = variable, values_from = perm_imp) %>%
  tibble::column_to_rownames("species") %>%
  as.matrix()
perm_mat <- perm_mat[sp_order, var_order]

# ============================================================================
# Part 3: PDP curves for ALL species × ALL variables
# ============================================================================
cat("Computing PDP curves for all species × all variables...\n")

make_response <- function(model, var_name, env_stack, n = 100) {
  r <- env_stack[[var_name]]
  v_range <- tryCatch({
    c(terra::global(r, "min", na.rm = TRUE)[1, 1],
      terra::global(r, "max", na.rm = TRUE)[1, 1])
  }, error = function(e) c(-5, 35))
  v_seq <- seq(v_range[1], v_range[2], length.out = n)
  
  fixed_vals <- sapply(names(env_stack), function(vn) {
    if (vn == var_name) return(NA)
    tryCatch(terra::global(env_stack[[vn]], "mean", na.rm = TRUE)[1, 1],
             error = function(e) 0)
  })
  
  pred_df <- data.frame(matrix(rep(fixed_vals, each = n), nrow = n,
                                dimnames = list(NULL, names(fixed_vals))))
  pred_df[[var_name]] <- v_seq
  pred_df <- pred_df[, names(env_stack), drop = FALSE]
  
  pred <- suppressWarnings(predict(model, pred_df, type = "cloglog"))
  data.frame(v = v_seq, p = pred)
}

sig_colors <- c("Positive" = "#d7191c", "Negative" = "#2c7bb6", "NS" = "gray70")

pdp_lookup  <- list()
pdp_raw     <- list()  # store raw p before scaling
sig_summary <- list()
counter     <- 0
total       <- length(sp_order) * length(var_order)

for (sp_name in sp_order) {
  pdp_lookup[[sp_name]] <- list()
  model <- model_list[[sp_name]]
  
  for (var_name in var_order) {
    counter <- counter + 1
    rdf <- tryCatch(make_response(model, var_name, env_final),
                    error = function(e) NULL)
    if (is.null(rdf)) next
    
    lm_fit <- lm(p ~ v, data = rdf)
    slope  <- coef(lm_fit)[2]
    
    pi_val <- perm_mat[sp_name, var_name]
    if (is.na(pi_val)) pi_val <- 0
    
    if (pi_val > 5 && slope > 0) {
      cls <- "Positive"
    } else if (pi_val > 5 && slope < 0) {
      cls <- "Negative"
    } else {
      cls <- "NS"
    }
    
    # Store raw p values for global scaling
    pdp_raw[[paste0(sp_name, "_", var_name)]] <- data.frame(
      species = sp_name, variable = var_name,
      v = rdf$v, p = rdf$p, color = sig_colors[cls],
      stringsAsFactors = FALSE)
    
    sig_summary[[paste0(sp_name, "_", var_name)]] <- data.frame(
      species = sp_name, variable = var_name,
      sig_class = cls, perm_imp = pi_val, slope = slope,
      stringsAsFactors = FALSE)
  }
  if (counter %% 10 == 0) cat(sprintf("  ...%d/%d\n", counter, total))
}
sig_df <- bind_rows(sig_summary)

# Global normalization: all curves share [0,1] based on global min/max p
x_ranges  <- list()          # var_name → c(min, max) of v
raw_all_p <- c()

for (key in names(pdp_raw)) {
  raw_all_p <- c(raw_all_p, pdp_raw[[key]]$p)
  var_name <- pdp_raw[[key]]$variable[1]
  if (is.null(x_ranges[[var_name]])) {
    x_ranges[[var_name]] <- range(pdp_raw[[key]]$v, na.rm = TRUE)
  }
}

global_min <- min(raw_all_p, na.rm = TRUE)
global_max <- max(raw_all_p, na.rm = TRUE)
cat(sprintf("  Global suitability range: [%.3f, %.3f]\n", global_min, global_max))

# Populate pdp_lookup with globally-scaled values
for (key in names(pdp_raw)) {
  rdf      <- pdp_raw[[key]]
  sp_name  <- rdf$species[1]
  var_name <- rdf$variable[1]
  p_scaled <- (rdf$p - global_min) / (global_max - global_min + 1e-10)
  pdp_lookup[[sp_name]][[var_name]] <- list(
    v = rdf$v, p_scaled = p_scaled, color = rdf$color[1])
}

cat(sprintf("  Computed %d PDP curves (global Y scaling 0–1)\n", nrow(sig_df)))

# ============================================================================
# Part 4: ComplexHeatmap — single integrated figure
# ============================================================================
cat("Building ComplexHeatmap with embedded PDP curves...\n")

col_fun <- colorRamp2(
  c(0, 10, 20, 35, 60),
  c("#f7fcf5", "#c7e9c0", "#74c476", "#238b45", "#00441b")
)

cell_fun <- function(j, i, x, y, width, height, fill) {
  sp_name  <- rownames(perm_mat)[i]
  var_name <- colnames(perm_mat)[j]
  
  pdp <- pdp_lookup[[sp_name]][[var_name]]
  if (!is.null(pdp) && length(pdp$v) > 1) {
    n <- length(pdp$v)
    
    # Convert unit positions to numeric (npc coordinates within cell)
    x_left  <- convertX(x - width/2,  "npc", valueOnly = TRUE)
    x_right <- convertX(x + width/2,  "npc", valueOnly = TRUE)
    y_bot   <- convertY(y - height/2, "npc", valueOnly = TRUE)
    y_top   <- convertY(y + height/2, "npc", valueOnly = TRUE)
    
    margin_x <- 0.12
    margin_y_bottom <- 0.26
    margin_y_top    <- 0.16
    
    cell_x <- seq(x_left + margin_x * (x_right - x_left),
                  x_right - margin_x * (x_right - x_left),
                  length.out = n)
    cell_y <- (y_bot + margin_y_bottom * (y_top - y_bot)) +
              pdp$p_scaled * ((y_top - y_bot) * (1 - margin_y_bottom - margin_y_top))
    
    grid.lines(unit(cell_x, "npc"), unit(cell_y, "npc"),
               gp = gpar(col = pdp$color, lwd = 1.0))
  }
  
  pi_val <- perm_mat[i, j]
  if (!is.na(pi_val)) {
    grid.text(sprintf("%.1f", pi_val),
              x = x, y = y - height/2 + unit(0.09, "npc"),
              gp = gpar(fontsize = 10, fontface = "bold", col = "#808080"),
              just = "bottom")
  }
}

# X-axis scale labels: keep min and max on one line.
x_scale_lab <- sapply(var_order, function(vn) {
  r <- x_ranges[[vn]]
  if (is.null(r)) return("")
  sprintf("%.0f \u2013 %.0f", r[1], r[2])
})

# Y-axis: global 0–1, single label for all rows
y_scale_lab <- sprintf("\u251C\u2500 0 \u2500\u2500\u2500 1 \u2524")

# Bottom annotation
ha_bottom <- HeatmapAnnotation(
  x_scale = anno_text(x_scale_lab[var_order],
                       gp = gpar(fontsize = 6.5, col = "#555555"),
                       location = 0.5, just = "centre",
                       rot = 0),
  which = "column",
  height = unit(0.55, "cm"),
  show_annotation_name = FALSE
)

# Right annotation: Y-axis (vertical scale, rotated 90°)
ha_right <- rowAnnotation(
  y_scale = anno_text(rep("0 \u2014 1", length(sp_order)),
                       gp = gpar(fontsize = 9, col = "#555555"),
                       location = 0.5, just = "centre",
                       rot = 90),
  width = unit(0.7, "cm"),
  show_annotation_name = FALSE
)

ht <- Heatmap(
  perm_mat,
  name = "Permutation importance (%)",
  col = col_fun,
  
  cluster_rows    = FALSE,
  cluster_columns = FALSE,
  
  row_labels    = sp_order,
  column_labels = var_labels[var_order],
  row_names_side   = "left",
  column_names_side = "top",
  row_names_gp    = gpar(fontsize = 9, fontface = "italic"),
  column_names_gp = gpar(fontsize = 6.5),
  column_names_rot = 45,
  
  width  = ncol(perm_mat) * unit(1.7, "cm"),
  height = nrow(perm_mat) * unit(1.6, "cm"),
  
  top_annotation    = NULL,
  bottom_annotation = ha_bottom,
  right_annotation  = ha_right,
  
  cell_fun = cell_fun,
  rect_gp  = gpar(col = "white", lwd = 1.5),
  
  show_heatmap_legend = FALSE
)

png(file.path(fig_dir, "Fig24_DriverResponse_Combined.png"),
    width = 9.5, height = 5.4, units = "in", res = 300)
draw(ht, padding = unit(c(3, 6, 3, 6), "mm"))
dev.off()

cat("  Fig. 24 saved.\n")

# ============================================================================
# Part 5: Significance summary
# ============================================================================
cat("Saving significance summary...\n")

sig_wide <- sig_df %>%
  mutate(label = case_when(
    sig_class == "Positive" ~ sprintf("+ PI=%.1f", perm_imp),
    sig_class == "Negative" ~ sprintf("- PI=%.1f", perm_imp),
    TRUE                   ~ sprintf("~ PI=%.1f", perm_imp)
  )) %>%
  dplyr::select(species, variable, label) %>%
  pivot_wider(names_from = variable, values_from = label)

write.csv(sig_wide, file.path(out_dir, "pdp_significance_summary.csv"), row.names = FALSE)
cat("\nSignificance summary:\n")
print.data.frame(sig_wide)

# Permutation importance matrix. If the existing CSV is open in Excel, keep the
# completed analysis and save the current matrix under a fallback name.
perm_matrix_file <- file.path(out_dir, "permutation_importance_matrix.csv")
tryCatch(
  write.csv(as.data.frame(perm_mat), perm_matrix_file),
  error = function(e) {
    fallback_file <- file.path(out_dir, "permutation_importance_matrix_new.csv")
    write.csv(as.data.frame(perm_mat), fallback_file)
    warning(
      "Could not overwrite ", perm_matrix_file,
      " (it may be open in Excel). Saved the current matrix to ", fallback_file,
      ". Close the old CSV before the next run."
    )
  }
)

cat("\n============================================================\n")
cat("  Done. Fig. 24 → ", file.path(fig_dir, "Fig24_DriverResponse_Combined.png"), "\n")
cat("============================================================\n")
