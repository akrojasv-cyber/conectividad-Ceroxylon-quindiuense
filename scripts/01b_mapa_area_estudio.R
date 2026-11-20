# =============================================================================
# 01b_mapa_area_estudio.R
# Mapa del área de estudio con las ocho poblaciones.
#
# Sección de la tesis: Metodología (área de estudio)
# Lee:      coordenadas escritas en el script; límites de GADM y Natural Earth (internet la primera vez)
# Produce:  figuras/mapa_area_estudio.png / .pdf
#
# Correr desde la carpeta Tesis_Ceroxylon, en una sesión limpia de R
# (Session > Restart R), después de los scripts con número menor.
# =============================================================================

source("scripts/00_configuracion.R")


# =============================================================================
# BLOQUE: Mapa del área de estudio
# (antes: 11_mapa_area_estudio.R)
# =============================================================================
# =============================================================================
# 11_mapa_area_estudio.R
# Mapa del área de estudio en tres paneles (A. América del Sur, B. Colombia,
# C. departamentos con las poblaciones), con límites oficiales y coordenadas
# reales. Necesita internet la primera vez (descarga GADM y Natural Earth).
# =============================================================================

library(sf); library(geodata); library(rnaturalearth)
library(ggplot2); library(ggspatial); library(ggrepel)
library(patchwork); library(dplyr)

# --- 1. Coordenadas de las poblaciones ---------------------------------------
# Usar EXACTAMENTE las del análisis de aislamiento por distancia (geosphere).
# Si ya tienes un data.frame con ellas, reemplaza este bloque.
dms_a_decimal <- function(g, m, s) g + m/60 + s/3600   # grados, minutos, segundos
pobs <- data.frame(
  codigo = c("SF","CX","CO","TO","LI","PI","AN","TE"),
  lat = c(dms_a_decimal(5,15,18), dms_a_decimal(4,38,15), dms_a_decimal(4,38,15),
          dms_a_decimal(4,31,0),  dms_a_decimal(4,30,47), dms_a_decimal(4,20,2),
          dms_a_decimal(4,19,6),  dms_a_decimal(3,43,50)),
  lon = -c(dms_a_decimal(75,24,45), dms_a_decimal(75,29,18), dms_a_decimal(75,31,40),
           dms_a_decimal(75,28,30), dms_a_decimal(75,35,37), dms_a_decimal(75,42,35),
           dms_a_decimal(75,29,48), dms_a_decimal(76,5,20)),
  fuente = c("GR","CV","GR","CV","CV","GR","GR","GR")
)
stopifnot(!anyNA(pobs$lon), !anyNA(pobs$lat))
pobs_sf <- st_as_sf(pobs, coords = c("lon","lat"), crs = 4326)

# --- 2. Límites ---------------------------------------------------------------
dir.create("capas", showWarnings = FALSE)
col_dep <- st_as_sf(gadm("COL", level = 1, path = "capas"))
sudam   <- ne_countries(continent = "South America", scale = "medium", returnclass = "sf")

deps_estudio <- c("Quindío","Tolima","Caldas","Valle del Cauca")
col_dep <- col_dep %>%
  mutate(depto = ifelse(NAME_1 %in% deps_estudio, NAME_1, "Otros departamentos"))

# Revisión: ¿en qué departamento cae cada punto? (debe coincidir con la Tabla)
print(st_join(pobs_sf, col_dep["NAME_1"]) %>% st_drop_geometry())

colores <- c("Quindío" = "#a8d5c2", "Tolima" = "#d6e2f5",
             "Caldas" = "#f0d5ec", "Valle del Cauca" = "#fbf1cf",
             "Otros departamentos" = "#e6e6e6")
col_fuente <- c(GR = "#2c3e9e", CV = "#c0392b")

# --- 3. Panel A: América del Sur ---------------------------------------------
pA <- ggplot() +
  geom_sf(data = sudam, fill = "#d9d9d9", colour = "white", linewidth = 0.2) +
  geom_sf(data = filter(sudam, admin == "Colombia"), fill = "#4f9a7f", colour = "grey30") +
  annotation_scale(location = "bl", width_hint = 0.4) +
  annotation_north_arrow(location = "tr", height = unit(0.8, "cm"), width = unit(0.6, "cm"),
                         style = north_arrow_orienteering(text_size = 7)) +
  labs(title = "A. América del Sur") + theme_void() +
  theme(plot.title = element_text(face = "bold", size = 10),
        panel.border = element_rect(colour = "grey50", fill = NA))

# --- 4. Panel B: Colombia con el recuadro del área ---------------------------
area <- st_bbox(filter(col_dep, NAME_1 %in% deps_estudio))
pB <- ggplot() +
  geom_sf(data = col_dep, aes(fill = depto), colour = "white", linewidth = 0.15,
          show.legend = FALSE) +
  scale_fill_manual(values = colores) +
  annotate("rect", xmin = area["xmin"], xmax = area["xmax"],
           ymin = area["ymin"], ymax = area["ymax"], fill = NA, colour = "grey20") +
  annotation_scale(location = "bl", width_hint = 0.4) +
  annotation_north_arrow(location = "tr", height = unit(0.8, "cm"), width = unit(0.6, "cm"),
                         style = north_arrow_orienteering(text_size = 7)) +
  labs(title = "B. Colombia") + theme_void() +
  theme(plot.title = element_text(face = "bold", size = 10),
        panel.border = element_rect(colour = "grey50", fill = NA))

# --- 5. Panel C: departamentos y poblaciones ---------------------------------
pC <- ggplot() +
  geom_sf(data = col_dep, aes(fill = depto), colour = "grey40", linewidth = 0.25) +
  scale_fill_manual(values = colores, name = "Departamento",
                    breaks = c(deps_estudio, "Otros departamentos")) +
  geom_sf(data = pobs_sf, aes(colour = fuente), size = 2.6) +
  scale_colour_manual(values = col_fuente, name = "Fuente",
                      labels = c(GR = "González-Rivillas et al. (2018)",
                                 CV = "Chacón-Vargas et al. (2019)")) +
  geom_text_repel(data = cbind(pobs, st_coordinates(pobs_sf)),
                  aes(X, Y, label = codigo), size = 3.4, fontface = "bold",
                  min.segment.length = 0, box.padding = 0.4, seed = 1) +
  coord_sf(xlim = c(area["xmin"], area["xmax"]), ylim = c(area["ymin"], area["ymax"]),
           expand = TRUE) +
  annotation_scale(location = "br", width_hint = 0.25) +
  annotation_north_arrow(location = "tr", height = unit(1, "cm"), width = unit(0.8, "cm"),
                         style = north_arrow_orienteering(text_size = 8)) +
  labs(title = "C. Departamentos del área de estudio", x = NULL, y = NULL,
       caption = "Sistema de coordenadas WGS84. Límites, GADM 4.1; continente, Natural Earth.") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 10),
        legend.position = "inside",
        legend.position.inside = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = "white", colour = "grey60", linewidth = 0.3),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        panel.grid = element_line(colour = "grey90", linewidth = 0.2))
# --- 6. Composición y guardado -----------------------------------------------
mapa <- (pA / pB | pC) + plot_layout(widths = c(1, 2.6))
ggsave("figuras/mapa_area_estudio.png", mapa, width = 16, height = 10.5, units = "cm", dpi = 400)
ggsave("figuras/mapa_area_estudio.pdf", mapa, width = 16, height = 10.5, units = "cm")
