# =============================================================================
# 00_configuracion.R
# Configuración común. Todos los scripts lo cargan al inicio con
#   source("scripts/00_configuracion.R")
#
# - Revisa que estén instalados los paquetes y, si falta alguno, lo instala
#   (una sola vez; diveRsity desde GitHub necesita Rtools en Windows).
# - Crea las carpetas outputs/ y figuras/.
# - No carga datos ni deja objetos de análisis en memoria.
#
# Estructura de la carpeta Tesis_Ceroxylon:
#   Matriz_datos_8poblaciones.xlsx   datos originales (las dos fuentes)
#   scripts/                         estos scripts
#   outputs/                         tablas y archivos intermedios
#   figuras/                         figuras de la tesis
#
# Orden de ejecución: 01, 01b, 02, 03, 04, 05, 06 y luego los anexos A1 a A4.
# 99_control.R compara los números clave con los de la tesis.
#
# RStudio: en Tools > Global Options > General, desmarcar
# "Restore .RData into workspace at startup" y poner "Never" en
# "Save workspace to .RData on exit", para que cada sesión arranque limpia.
# =============================================================================

paquetes_cran <- c("adegenet", "hierfstat", "poppr", "mmod", "vegan", "pegas",
                   "PopGenReport", "dunn.test", "ade4", "MASS", "igraph",
                   "expm", "geosphere", "readxl", "writexl", "openxlsx",
                   "dplyr", "tidyr", "ggplot2", "ggrepel", "patchwork",
                   "ggdendro", "ragg", "systemfonts", "sf", "geodata",
                   "rnaturalearth", "rnaturalearthdata", "ggspatial", "remotes")
faltan <- paquetes_cran[!paquetes_cran %in% rownames(installed.packages())]
if (length(faltan) > 0) install.packages(faltan, dependencies = TRUE)

if (!requireNamespace("LEA", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("LEA")
}
if (!requireNamespace("diveRsity", quietly = TRUE)) {
  remotes::install_github("kkeenan02/diveRsity")   # no está en CRAN para R 4.3
}

dir.create("outputs", showWarnings = FALSE)
dir.create("figuras", showWarnings = FALSE)
rm(paquetes_cran, faltan)

