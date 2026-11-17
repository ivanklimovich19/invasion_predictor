# запуск: Rscript R/03_occurrence_summary.R ambrosia
args <- commandArgs(trailingOnly = TRUE)
sp_id <- if (length(args)) args[1] else "ambrosia"

source("R/config.R")
source("R/functions.R")

cfg <- species[[sp_id]]
out <- file.path("output", sp_id)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

occ <- read_occ(cfg)
belarus <- get_belarus()

summ <- data.frame(
  records = nrow(occ),
  unique_coords = nrow(distinct(occ, longitude, latitude)),
  year_min = suppressWarnings(min(occ$year, na.rm = TRUE)),
  year_max = suppressWarnings(max(occ$year, na.rm = TRUE)))
if ("abundance" %in% names(occ)) {
  summ$abundance_mean <- round(mean(occ$abundance, na.rm = TRUE), 2)
  summ$abundance_median <- median(occ$abundance, na.rm = TRUE)
}
print(summ)
write.csv(summ, file.path(out, "summary.csv"), row.names = FALSE)

# накопление находок
by_year <- occ %>% filter(!is.na(year)) %>% count(year) %>% mutate(cum = cumsum(n))
p <- ggplot(by_year, aes(year, cum)) +
  geom_line(linewidth = 1.1, colour = "darkblue") +
  geom_point(colour = "red") +
  labs(title = bquote("Накопление находок" ~ italic(.(cfg$latin)) ~ "по годам"),
       x = "Год", y = "Число находок (накопленное)") +
  theme_minimal()
save_plot(p, file.path(out, "finds_cumulative.png"), 8, 5)

# обилие
if ("abundance" %in% names(occ)) {
  p <- ggplot() +
    geom_sf(data = belarus, fill = "grey85", colour = "black") +
    geom_point(data = occ, aes(longitude, latitude, size = abundance, colour = abundance),
               alpha = 0.7) +
    scale_size_continuous(range = c(1, 5), name = "Обилие") +
    scale_colour_gradient(low = "blue", high = "red", name = "Обилие") +
    coord_sf(xlim = c(23, 33), ylim = c(51, 56.5)) +
    labs(x = "Долгота", y = "Широта",
         title = bquote("Обилие" ~ italic(.(cfg$latin)))) +
    theme_minimal()
  save_plot(p, file.path(out, "abundance_map.png"), 8, 6)
}

# местообитания и антропогенность
share <- function(col) {
  occ %>%
    filter(!is.na(.data[[col]]), .data[[col]] != "NA") %>%
    group_by(across(all_of(col))) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>%
    arrange(desc(n))
}
for (col in c("habitat", "anthropogenization")) {
  if (col %in% names(occ)) {
    tab <- share(col)
    print(head(tab, 10))
    write.csv(tab, file.path(out, paste0(col, ".csv")), row.names = FALSE,
              fileEncoding = "UTF-8")
  }
}
