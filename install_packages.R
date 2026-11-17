pkgs <- c("geodata", "terra", "sf", "randomForest", "pROC", "pdp",
          "dplyr", "readr", "lubridate", "ggplot2",
          "rnaturalearth", "rnaturalearthdata", "adehabitatHR")

new <- setdiff(pkgs, rownames(installed.packages()))
if (length(new)) install.packages(new)
