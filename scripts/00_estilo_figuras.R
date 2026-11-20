# =============================================================================
# 00_estilo_figuras.R
# Estilo único para todas las figuras de la tesis. Cada script de figuras
# empieza con  source("00_estilo_figuras.R")  y guarda con  guardar().
#
# Reglas:
#   - Fuente Latin Modern Roman (la misma del documento, \usepackage{lmodern}).
#   - Títulos de ejes a 10 pt; números, leyendas y rótulos a 9 pt.
#   - Nada en negrilla. Sin título dentro de la figura (va en el \caption).
#   - Letras de panel "(a)", "(b)" en texto normal.
#   - Dos anchos: completo 15.4 cm (= \textwidth) y medio 11 cm.
#   - En LaTeX la figura se incluye a su tamaño real:
#       \includegraphics{figuras/nombre.pdf}        (sin width=)
#
# Fuente: instalar una sola vez "Latin Modern Roman" en el computador
# (archivos .otf de https://www.gust.org.pl/projects/e-foundry/latin-modern,
# clic derecho > Instalar). Si no está instalada, se usa "serif".
# =============================================================================

library(ggplot2)

fuentes_sistema <- unique(systemfonts::system_fonts()$family)
candidatas <- c("LMRoman10", "LM Roman 10", "Latin Modern Roman", "Latin Modern Roman 10")
FUENTE <- candidatas[candidatas %in% fuentes_sistema][1]
if (is.na(FUENTE)) {
  warning("Latin Modern no está instalada; se usa 'serif'.")
  FUENTE <- "serif"
}
cat("Fuente de las figuras:", FUENTE, "\n")

ANCHO_COMPLETO <- 15.4   # cm, igual a \textwidth
ANCHO_MEDIO    <- 11     # cm
TAM_BASE       <- 10     # pt, títulos de ejes
TAM_PEQ        <- 9      # pt, números, leyendas, rótulos

TINTA      <- "#3a3a38"
GRIS_GUIA  <- "grey88"
COL_FUENTE <- c(GR = "#2a78d6", CV = "#eb6834")   # validada para daltonismo
ETIQ_FUENTE <- c(GR = "González-Rivillas et al. (2018)",
                 CV = "Chacón-Vargas et al. (2019)")
AZUL       <- "#2a78d6"   # para figuras de una sola serie

tema_tesis <- function(base = TAM_BASE, peq = TAM_PEQ) {
  theme_minimal(base_size = base, base_family = FUENTE) +
    theme(
      text             = element_text(colour = TINTA, face = "plain"),
      axis.title       = element_text(size = base, face = "plain"),
      axis.text        = element_text(size = peq, colour = TINTA, face = "plain"),
      legend.title     = element_text(size = peq, face = "plain"),
      legend.text      = element_text(size = peq, face = "plain"),
      strip.text       = element_text(size = peq, face = "plain", colour = TINTA),
      plot.title       = element_text(size = base, face = "plain", hjust = 0),
      plot.tag         = element_text(size = base, face = "plain"),
      plot.caption     = element_text(size = peq - 1, face = "plain", colour = "grey40"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = GRIS_GUIA, linewidth = 0.3),
      axis.line        = element_blank(),
      legend.position  = "bottom",
      legend.key.size  = unit(0.4, "cm")
    )
}
theme_set(tema_tesis())
# Textos dibujados con geom_text / annotate en el mismo tamaño y fuente
update_geom_defaults("text",  list(family = FUENTE, size = TAM_PEQ / .pt, colour = TINTA))
update_geom_defaults("label", list(family = FUENTE, size = TAM_PEQ / .pt, colour = TINTA))

# Guardar siempre igual: PDF con fuentes incrustadas + PNG de respaldo
guardar <- function(fig, nombre, ancho = ANCHO_COMPLETO, alto) {
  dir.create("figuras", showWarnings = FALSE)
  ggsave(file.path("figuras", paste0(nombre, ".pdf")), fig,
         width = ancho, height = alto, units = "cm", device = cairo_pdf)
  ggsave(file.path("figuras", paste0(nombre, ".png")), fig,
         width = ancho, height = alto, units = "cm", dpi = 400)
  invisible(fig)
}

# Para figuras hechas con gráficos base de R (por ejemplo scatter de adegenet):
#   abrir_pdf("nombre", ancho, alto); ...dibujo...; dev.off()
abrir_pdf <- function(nombre, ancho = ANCHO_COMPLETO, alto) {
  dir.create("figuras", showWarnings = FALSE)
  cairo_pdf(file.path("figuras", paste0(nombre, ".pdf")),
            width = ancho / 2.54, height = alto / 2.54, family = FUENTE, pointsize = TAM_PEQ)
  par(family = FUENTE, cex.lab = TAM_BASE / TAM_PEQ, font.lab = 1, font.main = 1,
      col.axis = TINTA, col.lab = TINTA, mar = c(4, 4, 1, 1))
}

