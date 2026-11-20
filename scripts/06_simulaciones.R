# =============================================================================
# 06_simulaciones.R
# Pérdida de poblaciones (aleatoria y dirigida) y restauración de corredores.
#
# Sección de la tesis: 3.4 Escenarios de pérdida y restauración
# Lee:      edges_df_10loci.rds, outputs/metricas_red_10loci.xlsx
# Produce:  outputs/simulaciones_conservacion_10loci.xlsx, robustez_estatica, corredores_parejas
#
# Correr desde la carpeta Tesis_Ceroxylon, en una sesión limpia de R
# (Session > Restart R), después de los scripts con número menor.
# =============================================================================

source("scripts/00_configuracion.R")


# =============================================================================
# BLOQUE: Simulaciones
# (antes: simulaciones.R)
# =============================================================================
# ============================================================
#  SIMULACIÓN DE ESCENARIOS DE CONSERVACIÓN — Ceroxylon quindiuense
#  Script 06 — Objetivo 4: remoción de nodos, pérdida progresiva
#  de hábitat, y corredores de restauración
#
#  REQUIERE haber corrido antes: 03_red_conectividad_10loci.R
#  (debe existir: "edges_df_10loci.rds")
#  Y tener en memoria (o recargar) "gmat" (matriz completa de
#  migración relativa, de divMigrate — incluye pares NO significativos,
#  necesaria para simular corredores desde cero)
# ============================================================


# ── PASO 0 — PAQUETES ─────────────────────────────────────────

library(igraph)
library(openxlsx)
library(ggplot2)


# ── PASO 1 — RECONSTRUIR LA MATRIZ DE PESOS W ─────────────────

edges_df <- readRDS("edges_df_10loci.rds")
pop_order <- c("AN","CO","CX","LI","PI","SF","TE","TO")
n <- length(pop_order)

W0 <- matrix(0, nrow = n, ncol = n, dimnames = list(pop_order, pop_order))
for (i in seq_len(nrow(edges_df))) {
  W0[edges_df$from[i], edges_df$to[i]] <- edges_df$weight[i]
}

# gmat = matriz COMPLETA de migración relativa (28 pares, de divMigrate),
# necesaria para saber el valor base de un corredor que hoy no existe en W0.
# Si no la tienes en memoria, vuelve a cargarla como se hizo en el script 03.
if (!exists("gmat")) stop("Falta 'gmat' en memoria — vuelve a correr el bloque de divMigrate del script 03.")


# ── PASO 2 — FUNCIÓN DE EFICIENCIA GLOBAL (Latora & Marchiori, 2001) ──

eficiencia_global <- function(W_mat, pop_order) {
  n <- length(pop_order)
  g <- graph_from_adjacency_matrix(W_mat, mode = "directed", weighted = TRUE, diag = FALSE)
  if (ecount(g) == 0) return(0)
  E(g)$dist <- 1 / E(g)$weight
  d <- distances(g, mode = "out", weights = E(g)$dist)
  diag(d) <- NA
  d[!is.finite(d)] <- NA
  total <- sum(1 / d, na.rm = TRUE)
  total / (n * (n - 1))
}

E0 <- eficiencia_global(W0, pop_order)
cat("Eficiencia global de la red actual (E0):", round(E0, 4), "\n")


# Número de clases de comunicación de la cadena = componentes
# fuertemente conexos del grafo con arista donde W > 0
n_clases <- function(W) {
  g <- graph_from_adjacency_matrix((W > 0) * 1, mode = "directed")
  components(g, mode = "strong")$no
}

# Tiempo de relajación de la cadena en tiempo continuo:
# tau = 1 / |Re(lambda_2)|, con lambda_2 el valor propio no nulo
# del generador Q = W - D_out con mayor parte real
t_relajacion <- function(W) {
  Q  <- W - diag(rowSums(W))
  ev <- sort(Re(eigen(Q, only.values = TRUE)$values), decreasing = TRUE)
  if (abs(ev[2]) < 1e-10) return(Inf)   # más de una clase cerrada
  1 / abs(ev[2])
}


# ============================================================
# PASO 3 — ROBUSTEZ ESTÁTICA (reemplaza el bucle anterior)
# ============================================================

robustez <- data.frame()
for (nodo in pop_order) {
  W_removido <- W0
  W_removido[nodo, ] <- 0
  W_removido[, nodo] <- 0
  E_k     <- eficiencia_global(W_removido, pop_order)
  delta_E <- (E0 - E_k) / E0 * 100
  # el nodo quitado queda como una clase sola; se descuenta
  clases_k <- n_clases(W_removido) - 1
  robustez <- rbind(robustez, data.frame(
    Nodo                = nodo,
    Eficiencia_sin_nodo = round(E_k, 4),
    Delta_E_pct         = round(delta_E, 2),
    Clases_restantes    = clases_k))
}
robustez <- robustez[order(-robustez$Delta_E_pct), ]
cat("Clases en la red completa:", n_clases(W0), "\n")
print(robustez)

# ============================================================
#  PASO 4 — RESILIENCIA: ATAQUES DIRIGIDOS FRENTE A ALEATORIOS
#  (reemplaza el análisis que se hizo en Python)
#
#  Pegar en el script 06 después del PASO 3 (robustez estática).
#  Usa W0 (matriz de pesos de la red final) y pop_order, que ya
#  están creados en los PASOS 1 y 2.
# ============================================================

N0      <- length(pop_order)
W_total <- sum(W0)                 # suma de las 14 aristas dirigidas
cat("Peso total W_0:", round(W_total, 3), "\n")

# ── Indicadores para un conjunto de poblaciones que sobreviven ──

# S: fracción de las N0 poblaciones en el componente fuertemente
#    conexo más grande de la red que queda
indicador_S <- function(vivas) {
  if (length(vivas) == 0) return(0)
  W_sub <- W0[vivas, vivas, drop = FALSE]
  g <- graph_from_adjacency_matrix((W_sub > 0) * 1, mode = "directed")
  max(components(g, mode = "strong")$csize) / N0
}

# F: fracción del peso total que se conserva
indicador_F <- function(vivas) {
  if (length(vivas) == 0) return(0)
  sum(W0[vivas, vivas]) / W_total
}

# Área bajo la curva respecto a q = k/N0 (regla del trapecio)
auc <- function(y) sum((head(y, -1) + tail(y, -1)) / (2 * N0))


# ── Escenario dirigido: quitar la población de mayor fuerza total,
#    recalculada después de cada eliminación ─────────────────────

vivas   <- pop_order
orden   <- character(0)
S_dir   <- indicador_S(vivas)
F_dir   <- indicador_F(vivas)

for (k in 1:N0) {
  W_sub  <- W0[vivas, vivas, drop = FALSE]
  fuerza <- rowSums(W_sub) + colSums(W_sub)      # s_p = s_out + s_in
  quitar <- vivas[which.max(fuerza)]             # en empate, la primera en orden
  orden  <- c(orden, quitar)
  vivas  <- setdiff(vivas, quitar)
  S_dir  <- c(S_dir, indicador_S(vivas))
  F_dir  <- c(F_dir, indicador_F(vivas))
}
cat("Orden del ataque dirigido:", paste(orden, collapse = ", "), "\n")


# ── Escenario aleatorio: enumeración exhaustiva de todos los
#    subconjuntos de poblaciones eliminadas para cada k ───────────

aleatorio <- data.frame()
for (k in 0:N0) {
  quitadas <- if (k == 0) list(character(0)) else combn(pop_order, k, simplify = FALSE)
  S_k <- sapply(quitadas, function(q) indicador_S(setdiff(pop_order, q)))
  F_k <- sapply(quitadas, function(q) indicador_F(setdiff(pop_order, q)))
  aleatorio <- rbind(aleatorio, data.frame(
    k = k,
    S_media = mean(S_k), S_inf = quantile(S_k, 0.025), S_sup = quantile(S_k, 0.975),
    F_media = mean(F_k), F_inf = quantile(F_k, 0.025), F_sup = quantile(F_k, 0.975)))
}
rownames(aleatorio) <- NULL

# Comprobación: E[F(k)] = (N0-k)(N0-k-1) / (N0(N0-1))
F_teorico <- (N0 - 0:N0) * (N0 - 0:N0 - 1) / (N0 * (N0 - 1))
cat("Diferencia máxima con la fórmula cerrada:",
    max(abs(aleatorio$F_media - F_teorico)), "\n")


# ── Tablas de resultados ──────────────────────────────────────

resiliencia <- data.frame(
  k = 0:N0, q = round((0:N0) / N0, 3),
  S_dirigido  = round(S_dir, 3), S_aleatorio = round(aleatorio$S_media, 3),
  F_dirigido  = round(F_dir, 4), F_aleatorio = round(aleatorio$F_media, 4))
print(resiliencia)

auc_tabla <- data.frame(
  Indicador = c("Conectividad (S)", "Peso conservado (F)"),
  Dirigido  = round(c(auc(S_dir), auc(F_dir)), 3),
  Aleatorio = round(c(auc(aleatorio$S_media), auc(aleatorio$F_media)), 3))
print(auc_tabla)


# ── Figura: (a) S(k) y (b) F(k) ──────────────────────────────


library(patchwork)

# ── Tipo de cada corredor según las clases que une ───────────
nombre_clase <- c("Núcleo", "LI/TO", "CX")[match(
  clases0, c(clases0["CO"], clases0["LI"], clases0["CX"]))]
names(nombre_clase) <- pop_order

tipo_corredor <- function(par) {
  cl <- sort(nombre_clase[strsplit(par, "-")[[1]]])
  if (all(cl == c("CX", "Núcleo")))    return("CX con el núcleo")
  if (all(cl == c("LI/TO", "Núcleo"))) return("LI o TO con el núcleo")
  "CX con LI o TO"
}

uno <- subset(res_cor, n_corredores == 1)
uno$Tipo <- factor(sapply(uno$Corredores, tipo_corredor),
                   levels = c("LI o TO con el núcleo", "CX con el núcleo", "CX con LI o TO"))

colores_tipo <- c("LI o TO con el núcleo" = "#2171b5",
                  "CX con el núcleo"      = "#238b45",
                  "CX con LI o TO"        = "#cb181d")

tema_fig <- theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype = "dashed", linewidth = 0.3),
        legend.position = "inside",                 # ggplot2 >= 3.5
        legend.background = element_rect(color = "grey60", linewidth = 0.3),
        legend.text = element_text(size = 8),
        legend.key.width = unit(1.2, "lines"))

eje_t <- scale_x_log10(breaks = niveles_t, labels = niveles_t)

# ── (a) Un corredor ───────────────────────────────────────────
panel_a <- ggplot(uno, aes(t, Eficiencia, group = Corredores, color = Tipo)) +
  geom_hline(yintercept = E0, linetype = "dotted", color = "grey40") +
  annotate("text", x = 0.05, y = E0, label = "Red actual", hjust = 0, vjust = -0.5,
           size = 2.8, color = "grey40") +
  geom_line(linewidth = 0.6, alpha = 0.8) +
  geom_point(size = 1.4) +
  scale_color_manual(values = colores_tipo) +
  eje_t +
  labs(title = "(a) Un corredor", x = "Peso del corredor (t, escala log)",
       y = "Eficiencia global", color = NULL) +
  tema_fig +
  theme(legend.position.inside = c(0.02, 0.98), legend.justification = c(0, 1))

# ── (b) Dos corredores que unen toda la red ──────────────────
dos <- subset(res_cor, n_corredores == 2 & Clases == 1)
rango <- aggregate(Tau ~ t, data = dos, FUN = function(x) c(min = min(x), max = max(x)))
rango <- data.frame(t = rango$t, tau_min = rango$Tau[, "min"], tau_max = rango$Tau[, "max"])

mejor <- "CO-CX + CO-TO"     # menor tau con t = 0.05--0.2
peor  <- "CX-TO + LI-TE"     # mayor tau con t = 0.1--1
lineas_b <- rbind(
  data.frame(subset(dos, Corredores == mejor)[, c("t", "Tau")], Pareja = paste("Mejor:", mejor)),
  data.frame(subset(dos, Corredores == peor)[,  c("t", "Tau")], Pareja = paste("Peor:", peor)))

panel_b <- ggplot() +
  geom_ribbon(data = rango, aes(t, ymin = tau_min, ymax = tau_max),
              fill = "grey85") +
  geom_line(data = lineas_b, aes(t, Tau, color = Pareja), linewidth = 0.9) +
  geom_point(data = lineas_b, aes(t, Tau, color = Pareja, shape = Pareja), size = 2) +
  scale_color_manual(values = c("#2171b5", "#cb181d")) +
  scale_shape_manual(values = c(16, 15)) +
  eje_t +
  scale_y_log10() +
  labs(title = "(b) Dos corredores que unen toda la red",
       x = "Peso de cada corredor (t, escala log)",
       y = expression("Tiempo de relajación, " * tau * " (escala log)"),
       color = NULL, shape = NULL) +
  tema_fig +
  theme(legend.position.inside = c(0.98, 0.98), legend.justification = c(1, 1))

fig_corredores <- panel_a + panel_b
print(fig_corredores)

ggsave("outputs/corredores.png", fig_corredores, device = ragg::agg_png,
       width = 20, height = 9, units = "cm", dpi = 300)

# Para revisar: mejor y peor pareja con cada peso
print(rango)
dos_ord <- dos[order(dos$t, dos$Tau), ]
print(dos_ord[!duplicated(dos_ord$t), c("t", "Corredores", "Tau")])
print(dos_ord[!duplicated(dos_ord$t, fromLast = TRUE), c("t", "Corredores", "Tau")])

library(patchwork)

q <- (0:N0) / N0

panel_resiliencia <- function(y_dir, y_ale, titulo, eje_y) {
  datos <- rbind(
    data.frame(q = q, y = y_ale,
               escenario = sprintf("Aleatorio: media exacta (AUC = %.3f)", auc(y_ale))),
    data.frame(q = q, y = y_dir,
               escenario = sprintf("Dirigido por fuerza (AUC = %.3f)", auc(y_dir))))
  datos$escenario <- factor(datos$escenario, levels = unique(datos$escenario))
  
  ggplot(datos, aes(q, y, color = escenario, shape = escenario)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.2) +
    scale_color_manual(values = c("#2171b5", "#cb181d")) +
    scale_shape_manual(values = c(16, 15)) +
    scale_x_continuous(
      breaks = q, labels = c("0", "0.125", "0.25", "0.375", "0.5",
                             "0.625", "0.75", "0.875", "1"),
      expand = expansion(mult = 0.02),
      sec.axis = sec_axis(~ . * N0, breaks = 0:N0,
                          name = "Número de poblaciones eliminadas")) +
    scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(title = titulo, x = expression("Fracción de poblaciones eliminadas, " * q == k / N[0]),
         y = eje_y, color = NULL, shape = NULL) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, size = 11),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linetype = "dashed", linewidth = 0.3),
          legend.position = "inside",            # ggplot2 >= 3.5; en versiones
          legend.position.inside = c(0.98, 0.98),# anteriores: legend.position = c(0.98, 0.98)
          legend.justification = c(1, 1),
          legend.background = element_rect(color = "grey60", linewidth = 0.3),
          legend.text = element_text(size = 8),
          legend.key.width = unit(1.2, "lines"))
}

panel_a <- panel_resiliencia(S_dir, aleatorio$S_media,
                             "(a) Conservación de la conectividad",
                             expression(S(k) == "|" * C[max](k) * "|" / N[0]))
panel_b <- panel_resiliencia(F_dir, aleatorio$F_media,
                             "(b) Conservación del peso de las conexiones",
                             expression(F(k) == W(k) / W[0]))

fig_resiliencia <- panel_a + panel_b
print(fig_resiliencia)

# ragg suaviza las líneas (antialiasing); la ventana de RStudio en Windows
# no lo hace, por eso allí se ven pixeladas. Revisar el PNG guardado.
ggsave("outputs/resiliencia.png", fig_resiliencia, device = ragg::agg_png,
       width = 20, height = 9, units = "cm", dpi = 300)

# ── Agregar al Excel del PASO 6 ───────────────────────────────
#   addWorksheet(wb, "Resiliencia");     writeData(wb, "Resiliencia", resiliencia)
#   addWorksheet(wb, "Resiliencia_AUC"); writeData(wb, "Resiliencia_AUC", auc_tabla)
# ============================================================
# PASO 5 — CORREDORES DE RESTAURACIÓN (reemplaza el anterior)
# ============================================================

# Clases de la cadena observada
g0      <- graph_from_adjacency_matrix((W0 > 0) * 1, mode = "directed")
clases0 <- components(g0, mode = "strong")$membership
names(clases0) <- pop_order

# Candidatos: todos los pares de poblaciones en clases distintas
pares      <- t(combn(pop_order, 2))
candidatos <- pares[clases0[pares[, 1]] != clases0[pares[, 2]], , drop = FALSE]
cat("Número de corredores candidatos:", nrow(candidatos), "\n")

# Peso del corredor: fracción t del flujo máximo, igual en los dos sentidos
niveles_t <- c(0.05, 0.1, 0.2, 0.5, 1)

agregar <- function(W, cor, t) {
  W[cor[1], cor[2]] <- t
  W[cor[2], cor[1]] <- t
  W
}

evaluar <- function(W_t, t, n_cor, etiqueta) {
  k <- n_clases(W_t)
  data.frame(t            = t,
             n_corredores = n_cor,
             Corredores   = etiqueta,
             Clases       = k,
             Eficiencia   = round(eficiencia_global(W_t, pop_order), 4),
             Tau          = if (k == 1) round(t_relajacion(W_t), 2) else NA)
}

res_cor <- data.frame()
combs   <- combn(nrow(candidatos), 2)

for (t in niveles_t) {
  # Un corredor
  for (i in seq_len(nrow(candidatos))) {
    a <- candidatos[i, ]
    res_cor <- rbind(res_cor,
                     evaluar(agregar(W0, a, t), t, 1, paste(a, collapse = "-")))
  }
  # Dos corredores
  for (j in seq_len(ncol(combs))) {
    a <- candidatos[combs[1, j], ]
    b <- candidatos[combs[2, j], ]
    W_t <- agregar(agregar(W0, a, t), b, t)
    res_cor <- rbind(res_cor,
                     evaluar(W_t, t, 2, paste(paste(a, collapse = "-"),
                                              paste(b, collapse = "-"), sep = " + ")))
  }
}

# Resumen: escenarios que dejan la cadena con una sola clase,
# ordenados por tiempo de relajación (menor = integra más rápido)
irreducibles <- subset(res_cor, Clases == 1)
irreducibles <- irreducibles[order(irreducibles$t, irreducibles$Tau), ]
print(head(irreducibles[irreducibles$t == 0.1, ], 10))
print(tail(irreducibles[irreducibles$t == 0.1, ], 5))

# Gráfico: eficiencia según t para los corredores individuales
ggplot(subset(res_cor, n_corredores == 1),
       aes(x = t, y = Eficiencia, color = Corredores)) +
  geom_line(linewidth = 0.8) + geom_point() +
  geom_hline(yintercept = E0, linetype = "dashed", color = "gray40") +
  labs(x = "Peso del corredor (t)", y = "Eficiencia global") +
  theme_minimal()


# ── PASO 6 — EXPORTAR (versión nueva) ────────────────────────
wb <- createWorkbook()
addWorksheet(wb, "Robustez_estatica"); writeData(wb, "Robustez_estatica", robustez)
addWorksheet(wb, "Corredores");        writeData(wb, "Corredores", res_cor)
addWorksheet(wb, "Resiliencia");     writeData(wb, "Resiliencia", resiliencia)
addWorksheet(wb, "Resiliencia_AUC"); writeData(wb, "Resiliencia_AUC", auc_tabla)
if (!dir.exists("outputs")) dir.create("outputs")
saveWorkbook(wb, "outputs/simulaciones_conservacion_10loci.xlsx", overwrite = TRUE)


# =============================================================================
# BLOQUE: Figuras de simulaciones
# (antes: 19_figuras_simulaciones.R)
# =============================================================================
rm(list = ls())   # memoria limpia para este bloque
# =============================================================================
# 19_figuras_simulaciones.R
# Figuras nuevas de las simulaciones, con el estilo común.
#
#   1. robustez_estatica  : (a) caída de la eficiencia al quitar cada población
#                           (b) caída frente a la fuerza total
#   2. corredores_parejas : mapa de calor del tiempo de relajación tau para las
#                           136 parejas de corredores (t = 0.1)
#
# Entradas: outputs/simulaciones_conservacion_10loci.xlsx (script 06)
#           outputs/metricas_red_10loci.xlsx (métricas)
# Correr desde Tesis_Ceroxylon.
# =============================================================================

source("scripts/00_estilo_figuras.R")
library(readxl); library(dplyr); library(tidyr); library(ggrepel)
library(patchwork)

rob  <- read_excel("outputs/simulaciones_conservacion_10loci.xlsx", sheet = "Robustez_estatica")
cor  <- read_excel("outputs/simulaciones_conservacion_10loci.xlsx", sheet = "Corredores")
nodo <- read_excel("outputs/metricas_red_10loci.xlsx", sheet = "Metricas_nodo")

# ------------------------------------------------------- 1. Robustez ---------
rob <- rob %>%
  left_join(select(nodo, Nodo = Poblacion, Fuerza_total), by = "Nodo") %>%
  mutate(efecto = ifelse(Clases_restantes > 3, "Parte la red", "No parte la red"))
col_efecto <- c("Parte la red" = "#c0392b", "No parte la red" = AZUL)

pa <- ggplot(rob, aes(x = Delta_E_pct, y = reorder(Nodo, Delta_E_pct), fill = efecto)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.1f", Delta_E_pct)), hjust = -0.2,
            family = FUENTE, size = TAM_PEQ / .pt, colour = TINTA) +
  scale_fill_manual(values = col_efecto, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Caída de la eficiencia global (%)", y = NULL, tag = "(a)") +
  theme(panel.grid.major.y = element_blank())

# Etiquetas: LI y TO caen en el mismo punto, se juntan; posición manual sin líneas
etq <- rob %>%
  group_by(Fuerza_total, Delta_E_pct, efecto) %>%
  summarise(lab = paste(sort(Nodo), collapse = ", "), .groups = "drop") %>%
  mutate(hj = ifelse(lab == "CO", 1.25, -0.25),          # CO a la izquierda (borde)
         vj = case_when(lab == "CX" ~ 1.2, lab == "LI, TO" ~ -0.4, TRUE ~ 0.5))

pb <- ggplot(rob, aes(x = Fuerza_total, y = Delta_E_pct)) +
  geom_point(aes(fill = efecto), shape = 21, size = 2.6, colour = "white", stroke = 0.5) +
  geom_text(data = etq, aes(label = lab, hjust = hj, vjust = vj),
            family = FUENTE, size = TAM_PEQ / .pt, colour = TINTA) +
  scale_fill_manual(values = col_efecto, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.05))) +
  labs(x = "Fuerza total de la población quitada",
       y = "Caída de eficiencia (%)", tag = "(b)")

fig1 <- (pa | pb) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom", plot.tag = element_text(size = TAM_BASE))
guardar(fig1, "robustez_estatica", ANCHO_COMPLETO, 7.5)

# ------------------------------------------ 2. Parejas de corredores ---------
t_sel  <- 0.1
nucleo <- c("AN", "CO", "PI", "SF", "TE")
tipo_corredor <- function(x) {
  p <- strsplit(x, "-")[[1]]
  if ("CX" %in% p && any(p %in% nucleo)) "CX–núcleo"
  else if (any(p %in% c("LI", "TO")) && any(p %in% nucleo)) "LI/TO–núcleo"
  else "CX–LI/TO"
}

uno <- cor %>% filter(n_corredores == 1, t == t_sel) %>% pull(Corredores)
orden <- tibble(c = uno, tipo = sapply(uno, tipo_corredor)) %>%
  mutate(tipo = factor(tipo, levels = c("CX–núcleo", "LI/TO–núcleo", "CX–LI/TO"))) %>%
  arrange(tipo, c) %>% pull(c)
pos <- setNames(seq_along(orden), orden)
n_tipo <- table(sapply(orden, tipo_corredor))[c("CX–núcleo", "LI/TO–núcleo", "CX–LI/TO")]

# Cada pareja en el triángulo inferior: fila = el corredor que va después en el orden
dos <- cor %>%
  filter(n_corredores == 2, t == t_sel) %>%
  separate(Corredores, into = c("c1", "c2"), sep = " \\+ ") %>%
  mutate(fila = ifelse(pos[c1] > pos[c2], c1, c2),
         col  = ifelse(pos[c1] > pos[c2], c2, c1))

etiq <- function(x) gsub("-", "–", x)
lim_x <- orden[-length(orden)]
lim_y <- rev(orden[-1])
cortes_x <- cumsum(n_tipo)[1:2] + 0.5                       # 5.5 y 15.5
cortes_y <- length(lim_y) - cumsum(n_tipo)[1:2] + 1.5       # 12.5 y 2.5
# Bordes de bloque solo dentro del triángulo
bordes <- data.frame(x    = c(cortes_x, 0.5, 0.5),
                     xend = c(cortes_x, cortes_x),
                     y    = c(0.5, 0.5, cortes_y),
                     yend = c(cortes_y, cortes_y))

fig2 <- ggplot(dos, aes(x = col, y = fila, fill = Tau)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_segment(data = bordes, aes(x = x, xend = xend, y = y, yend = yend),
               inherit.aes = FALSE, colour = "grey40", linewidth = 0.3) +
  scale_x_discrete(limits = lim_x, labels = etiq) +
  scale_y_discrete(limits = lim_y, labels = etiq) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(9, "YlOrRd")[3:9],
                       na.value = "grey88", name = expression(tau)) +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = TAM_PEQ - 1),
        axis.text.y = element_text(size = TAM_PEQ - 1),
        panel.grid.major = element_blank(),
        legend.position = "right",
        legend.key.height = unit(1.2, "cm"))
guardar(fig2, "corredores_parejas", ANCHO_COMPLETO, 13)

# Comprobación: deben salir 80 parejas que unen la red y 56 que no
cat("Parejas que unen la red:", sum(!is.na(dos$Tau)), " | que no:", sum(is.na(dos$Tau)), "\n")
