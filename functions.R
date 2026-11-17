library(terra)
library(sf)
library(dplyr)
library(readr)
library(lubridate)
library(randomForest)
library(pROC)
library(pdp)
library(ggplot2)
library(geodata)
library(rnaturalearth)

num <- function(x) as.numeric(gsub(",", ".", x))

read_occ <- function(cfg) {
  d <- read_delim(cfg$file, delim = ";", col_types = cols(.default = "c"),
                  locale = locale(encoding = cfg$encoding),
                  show_col_types = FALSE)
  names(d) <- tolower(trimws(names(d)))
  names(d)[names(d) == "longtude"] <- "longitude"

  d$longitude <- num(d$longitude)
  d$latitude <- num(d$latitude)
  if ("abundance" %in% names(d)) d$abundance <- num(d$abundance)
  if ("date" %in% names(d)) {
    y <- year(dmy(d$date, quiet = TRUE))
    y_only <- suppressWarnings(as.numeric(ifelse(grepl("^[0-9]{4}$", d$date), d$date, NA)))
    d$year <- ifelse(is.na(y), y_only, y)
  }

  d %>% filter(!is.na(longitude), !is.na(latitude),
               between(longitude, bbox[1], bbox[2]),
               between(latitude, bbox[3], bbox[4]))
}

load_bioclim <- function(res) {
  r <- worldclim_global(var = "bio", res = res, path = "data/worldclim")
  r <- crop(r, ext(bbox))
  names(r) <- paste0("bio", 1:19)
  r
}

get_belarus <- function() {
  ne_countries(scale = "medium", country = "Belarus", returnclass = "sf")
}

sample_bg <- function(r, n, occ, seed) {
  set.seed(seed)
  bg <- spatSample(r[[1]], size = n, method = "random", na.rm = TRUE,
                   xy = TRUE, values = FALSE)
  bg <- as.data.frame(bg)
  names(bg) <- c("longitude", "latitude")
  busy <- cellFromXY(r, as.matrix(occ[, c("longitude", "latitude")]))
  bg[!(cellFromXY(r, as.matrix(bg)) %in% busy), ]
}

make_training <- function(occ, bg, r) {
  pts <- bind_rows(
    occ %>% distinct(longitude, latitude) %>% mutate(class = 1),
    bg %>% mutate(class = 0))
  env <- terra::extract(r, as.matrix(pts[, c("longitude", "latitude")]), ID = FALSE)
  d <- cbind(pts, env)
  d <- d[complete.cases(d), ]
  d$class <- factor(d$class, levels = c(0, 1))
  d
}

predictors <- function(d) d[, grep("^bio", names(d))]

fit_rf <- function(d, ntree, seed) {
  set.seed(seed)
  n1 <- sum(d$class == 1)
  randomForest(x = predictors(d), y = d$class, ntree = ntree,
               importance = TRUE, strata = d$class, sampsize = c(n1, n1))
}

auc_of <- function(y, p) {
  as.numeric(auc(roc(y, p, levels = c("0", "1"), direction = "<", quiet = TRUE)))
}

block_cv_auc <- function(d, ntree, seed, size_km, k) {
  v <- vect(d[, c("longitude", "latitude")], geom = c("longitude", "latitude"),
            crs = "EPSG:4326")
  xy <- crds(project(v, "EPSG:32635"))
  block <- paste(floor(xy[, 1] / (size_km * 1000)), floor(xy[, 2] / (size_km * 1000)))

  set.seed(seed)
  ub <- unique(block)
  fold_of <- setNames(sample(rep(seq_len(k), length.out = length(ub))), ub)
  fold <- unname(fold_of[block])

  res <- sapply(seq_len(k), function(i) {
    tr <- d[fold != i, ]
    te <- d[fold == i, ]
    if (sum(tr$class == 1) < 5 || length(unique(te$class)) < 2) return(NA_real_)
    auc_of(te$class, predict(fit_rf(tr, ntree, seed), predictors(te), type = "prob")[, "1"])
  })
  mean(res, na.rm = TRUE)
}

importance_table <- function(rf) {
  imp <- importance(rf, type = 1)
  data.frame(var = rownames(imp), mda = imp[, 1]) %>%
    mutate(pct = 100 * mda / max(mda), name = unname(bio_ru[var])) %>%
    arrange(desc(mda))
}

response_curve <- function(rf, d, var) {
  pd <- partial(rf, pred.var = var, train = predictors(d), prob = TRUE,
                which.class = 2, grid.resolution = 50)
  data.frame(x = pd[[var]], p = pd$yhat)
}

curve_range <- function(pd, thr = 0.5) {
  a <- pd[pd$p > thr, ]
  data.frame(from = if (nrow(a)) min(a$x) else NA,
             to = if (nrow(a)) max(a$x) else NA,
             peak_x = pd$x[which.max(pd$p)], peak_p = max(pd$p))
}

predict_map <- function(r, rf, belarus) {
  p <- predict(r, rf, type = "prob", index = 2, na.rm = TRUE)
  names(p) <- "p"
  mask(p, vect(belarus))
}

zone_areas <- function(p, breaks) {
  m <- cbind(breaks[-length(breaks)], breaks[-1], 1:(length(breaks) - 1))
  cls <- classify(p, m, include.lowest = TRUE)
  a <- expanse(cls, unit = "km", byValue = TRUE)
  data.frame(zone = c("низкая", "средняя", "высокая")[a$value], area_km2 = round(a$area))
}

map_plot <- function(p, occ, belarus, latin) {
  df <- as.data.frame(p, xy = TRUE)
  ggplot() +
    geom_raster(data = df, aes(x, y, fill = p)) +
    scale_fill_viridis_c(option = "plasma", limits = c(0, 1),
                         name = "Вероятность\nприсутствия") +
    geom_sf(data = belarus, fill = NA, colour = "black", linewidth = 0.5) +
    geom_point(data = occ, aes(longitude, latitude), colour = "red",
               size = 1, alpha = 0.6) +
    coord_sf(xlim = c(23, 33), ylim = c(51, 56.5)) +
    labs(title = bquote("Прогноз распространения" ~ italic(.(latin))),
         x = "Долгота", y = "Широта") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

importance_plot <- function(imp) {
  imp$name <- factor(imp$name, levels = rev(imp$name))
  ggplot(imp, aes(pct, name)) +
    geom_segment(aes(x = 0, xend = pct, yend = name), colour = "grey50") +
    geom_point(size = 2.5) +
    scale_x_continuous(limits = c(0, 105), expand = c(0, 0)) +
    labs(x = "Важность переменной, % от максимума", y = NULL) +
    theme_classic(base_size = 12)
}

curve_plot <- function(pd, var, d) {
  ggplot(pd, aes(x, p)) +
    geom_line(linewidth = 1) +
    geom_rug(data = data.frame(x = d[[var]]), aes(x), inherit.aes = FALSE,
             sides = "b", alpha = 0.2) +
    labs(title = paste("Кривая отклика:", var), x = var,
         y = "Вероятность присутствия") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

save_plot <- function(p, file, w = 9, h = 7) {
  ggsave(file, p, width = w, height = h, dpi = 300)
}
