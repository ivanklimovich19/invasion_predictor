# запуск: Rscript R/01_model.R ambrosia
args <- commandArgs(trailingOnly = TRUE)
sp_id <- if (length(args)) args[1] else "ambrosia"

source("R/config.R")
source("R/functions.R")

cfg <- species[[sp_id]]
out <- file.path("output", sp_id)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

occ <- read_occ(cfg)
r <- load_bioclim(cfg$res)
bg <- sample_bg(r, cfg$n_bg, occ, seed)
train <- make_training(occ, bg, r)
cat(cfg$latin, ": присутствий", sum(train$class == 1),
    ", фоновых", sum(train$class == 0), "\n")

rf <- fit_rf(train, ntree, seed)

metrics <- data.frame(
  species = sp_id,
  n_presence = sum(train$class == 1),
  n_background = sum(train$class == 0),
  oob_error = rf$err.rate[ntree, "OOB"],
  auc_oob = auc_of(train$class, rf$votes[, "1"]),
  auc_block_cv = block_cv_auc(train, ntree, seed, block_km, cv_folds))
print(metrics)
write.csv(metrics, file.path(out, "metrics.csv"), row.names = FALSE)

imp <- importance_table(rf)
write.csv(imp, file.path(out, "importance.csv"), row.names = FALSE,
          fileEncoding = "UTF-8")
save_plot(importance_plot(imp), file.path(out, "importance.png"), 8, 6.5)

top <- imp$var[1:3]
ranges <- lapply(top, function(v) {
  pd <- response_curve(rf, train, v)
  save_plot(curve_plot(pd, v, train), file.path(out, paste0("response_", v, ".png")), 6, 4.5)
  cbind(var = v, curve_range(pd))
})
ranges <- do.call(rbind, ranges)
print(ranges)
write.csv(ranges, file.path(out, "response_ranges.csv"), row.names = FALSE)

belarus <- get_belarus()
pred <- predict_map(r, rf, belarus)
writeRaster(pred, file.path(out, "prediction.tif"), overwrite = TRUE)
save_plot(map_plot(pred, occ, belarus, cfg$latin), file.path(out, "prediction_map.png"))

zones <- zone_areas(pred, zone_breaks)
print(zones)
write.csv(zones, file.path(out, "zone_areas.csv"), row.names = FALSE,
          fileEncoding = "UTF-8")
