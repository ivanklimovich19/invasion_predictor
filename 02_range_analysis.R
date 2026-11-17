# запуск: Rscript R/02_range_analysis.R ambrosia
args <- commandArgs(trailingOnly = TRUE)
sp_id <- if (length(args)) args[1] else "ambrosia"

source("R/config.R")
source("R/functions.R")
library(sp)
library(adehabitatHR)

cfg <- species[[sp_id]]
out <- file.path("output", sp_id)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

occ <- read_occ(cfg) %>% filter(!is.na(year))
pts <- st_as_sf(occ, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
utm <- st_transform(pts, 32635)
belarus <- get_belarus()
to_wgs <- function(x) st_transform(x, 4326)

base_map <- function() {
  ggplot() +
    geom_sf(data = belarus, fill = "grey85", colour = "black") +
    coord_sf(xlim = c(23, 33), ylim = c(51, 56.5), expand = FALSE) +
    labs(x = "Долгота", y = "Широта") +
    theme_minimal()
}

# общий ареал и ядра концентрации
xy <- as(st_geometry(utm), "Spatial")
kde <- kernelUD(xy, h = "href")
k95 <- to_wgs(st_as_sf(getverticeshr(kde, percent = 95)))
k50 <- to_wgs(st_as_sf(getverticeshr(kde, percent = 50)))
mcp_all <- st_convex_hull(st_union(utm))

p <- base_map() +
  geom_sf(data = k95, fill = "green", alpha = 0.2, colour = NA) +
  geom_sf(data = k50, fill = "darkgreen", alpha = 0.4, colour = NA) +
  geom_sf(data = to_wgs(mcp_all), fill = NA, colour = "blue", linetype = "dashed") +
  geom_sf(data = pts, colour = "red", size = 1.5, alpha = 0.7) +
  labs(title = bquote("Ареал и ядра распространения" ~ italic(.(cfg$latin))))
save_plot(p, file.path(out, "range_map.png"), 8, 6)

# динамика по периодам
last_year <- max(occ$year)
labs_p <- c("до 2000", "2000-2009", "2010-2019", paste0("2020-", last_year))
utm$period <- cut(utm$year, c(-Inf, 1999, 2009, 2019, Inf), labels = labs_p)
pts$period <- utm$period

polys <- list()
stats <- lapply(labs_p, function(l) {
  s <- utm[utm$period == l, ]
  n <- nrow(s)
  if (n < 5) return(data.frame(period = l, n = n, area_km2 = NA, lon = NA, lat = NA))
  h <- st_convex_hull(st_union(s))
  polys[[l]] <<- st_sf(period = l, geometry = h)
  cen <- st_coordinates(to_wgs(st_centroid(h)))
  data.frame(period = l, n = n, area_km2 = round(as.numeric(st_area(h)) / 1e6),
             lon = round(cen[1], 2), lat = round(cen[2], 2))
})
stats <- do.call(rbind, stats)
print(stats)
write.csv(stats, file.path(out, "range_periods.csv"), row.names = FALSE,
          fileEncoding = "UTF-8")

chi <- chisq.test(stats$n)
cat("chi2 =", round(chi$statistic, 1), " df =", chi$parameter,
    " p =", format.pval(chi$p.value, digits = 3), "\n")
ok <- stats[!is.na(stats$area_km2), ]
if (nrow(ok) >= 3) {
  sp_cor <- suppressWarnings(cor.test(seq_len(nrow(ok)), ok$area_km2, method = "spearman"))
  cat("spearman rho =", round(sp_cor$estimate, 2),
      " p =", format.pval(sp_cor$p.value, digits = 3), "\n")
}

if (length(polys)) {
  pl <- to_wgs(do.call(rbind, polys))
  p <- base_map() +
    geom_sf(data = pl, aes(fill = period), alpha = 0.25, colour = "black") +
    geom_sf(data = pts, aes(colour = period), size = 1.8, alpha = 0.8) +
    scale_fill_brewer(palette = "Set2", name = "Период") +
    scale_colour_brewer(palette = "Set2", name = "Период") +
    labs(title = "Динамика ареала по периодам")
  save_plot(p, file.path(out, "range_dynamics.png"), 10, 7)
}

# ядра последнего периода
recent <- utm[utm$period == labs_p[4], ]
if (nrow(recent) >= 5) {
  kr <- kernelUD(as(st_geometry(recent), "Spatial"), h = "href")
  r95 <- to_wgs(st_as_sf(getverticeshr(kr, percent = 95)))
  r50 <- to_wgs(st_as_sf(getverticeshr(kr, percent = 50)))
  p <- base_map() +
    geom_sf(data = r95, fill = "green", alpha = 0.2, colour = NA) +
    geom_sf(data = r50, fill = "darkgreen", alpha = 0.4, colour = NA) +
    geom_sf(data = pts[pts$period == labs_p[4], ], colour = "red", size = 1.8) +
    labs(title = paste("Ядра концентрации,", labs_p[4]))
  save_plot(p, file.path(out, "kde_latest.png"), 8, 6)
}
