# =============================================================================
# 05_metricas_markov.R
# Métricas de nodo y de red, Mantel (F_ST y red contra distancia) y modelo de Markov.
#
# Sección de la tesis: 3.3 Métricas del grafo, aislamiento por distancia y dinámica de Markov
# Lee:      red_10loci.rds, edges_df_10loci.rds
# Produce:  outputs/metricas_red_10loci.xlsx, outputs/modelo_dinamico_10loci.xlsx, P_markov_10loci.rds, Q_laplaciano_10loci.rds, metricas_perfil_nodos, markov_visitas_permanencia, mantel_ibd
#
# Correr desde la carpeta Tesis_Ceroxylon, en una sesión limpia de R
# (Session > Restart R), después de los scripts con número menor.
# =============================================================================

source("scripts/00_configuracion.R")


# =============================================================================
# BLOQUE: Métricas de la red y aislamiento por distancia
# (antes: METRICAS_RED Y AISLAMIENTO POR DISTANCIA.R)
# =============================================================================
# ============================================================
#  MÉTRICAS DE LA RED Y AISLAMIENTO POR DISTANCIA
#  Script 03b — Objetivo 3 (10 loci corregidos)
#
#  REQUIERE haber corrido antes: 03a_construccion_red_10loci.R
#  (debe existir: "red_10loci.rds")
# ============================================================

library(igraph); library(ade4); library(writexl)

red <- readRDS("red_10loci.rds")
list2env(red, envir = .GlobalEnv)   # pop_order, gmat, edges_df, aristas_finales, ...


# ── 0. VERSIONES DE LA RED ───────────────────────────────────

# Dirigida: peso = migración relativa; longitud = 1/peso
g_dir <- graph_from_data_frame(edges_df, directed = TRUE, vertices = pop_order)
E(g_dir)$dist <- 1 / E(g_dir)$weight

# No dirigida: las 7 aristas de E_f, con dos pesos posibles
g_top <- graph_from_data_frame(aristas_finales[, c("Pop1", "Pop2")],
                               directed = FALSE, vertices = pop_order)
E(g_top)$r       <- aristas_finales$r
E(g_top)$m_media <- mapply(function(a, b) (gmat[a, b] + gmat[b, a]) / 2,
                           aristas_finales$Pop1, aristas_finales$Pop2)

# Componentes conexos (débil = fuerte, porque cada arista va en los dos sentidos)
comp <- components(g_top)
print(split(names(comp$membership), comp$membership))
nucleo <- names(comp$membership)[comp$membership == which.max(comp$csize)]
g_nuc  <- induced_subgraph(g_dir, nucleo)


# ============================================================
# MÉTRICAS DE NODO
# ============================================================

# Fuerza de nodo
fuerza_total <- fuerza_sal + fuerza_ent

# PageRank ponderado (alfa = 0.85 por defecto)
pagerank        <- page_rank(g_dir, weights = E(g_dir)$weight)$vector
pagerank_nucleo <- page_rank(g_nuc, weights = E(g_nuc)$weight)$vector

# Comparación con la distribución estacionaria discreta del núcleo (pi^(d))
W_nuc <- as_adjacency_matrix(g_nuc, attr = "weight", sparse = FALSE)
P_nuc <- W_nuc / rowSums(W_nuc)
ev    <- eigen(t(P_nuc))
pi_markov <- Re(ev$vectors[, which.min(abs(ev$values - 1))])
pi_markov <- pi_markov / sum(pi_markov)
names(pi_markov) <- rownames(P_nuc)
round(rbind(PageRank = pagerank_nucleo[nucleo], pi_d = pi_markov[nucleo]), 4)

# Betweenness ponderado (longitud = 1/peso, dirigido, sin normalizar)
betweenness_pond <- betweenness(g_dir, weights = E(g_dir)$dist, directed = TRUE)

# Cercanía armónica (entrante y saliente)
harm_ent <- harmonic_centrality(g_dir, mode = "in",  weights = E(g_dir)$dist, normalized = TRUE)
harm_sal <- harmonic_centrality(g_dir, mode = "out", weights = E(g_dir)$dist, normalized = TRUE)

# Coreness
coreness_k <- coreness(g_top)

# Clustering local
trans_local <- transitivity(g_top, type = "local", isolates = "zero")
names(trans_local) <- V(g_top)$name


# ============================================================
# MÉTRICAS DE RED
# ============================================================

# Eficiencia global (dirigida, longitud = 1/peso, sin camino = 0)
eficiencia_global <- function(g) {
  n <- vcount(g)
  d <- distances(g, mode = "out", weights = 1 / E(g)$weight)
  diag(d) <- NA
  d[!is.finite(d)] <- NA
  sum(1 / d, na.rm = TRUE) / (n * (n - 1))
}
E0 <- eficiencia_global(g_dir)

# Transitividad global
trans_global <- transitivity(g_top, type = "global")

# Modularidad (Louvain) con los dos pesos, en la red completa y en el núcleo
set.seed(2026)
louv_r <- cluster_louvain(g_top, weights = E(g_top)$r)
louv_m <- cluster_louvain(g_top, weights = E(g_top)$m_media)
g_top_nuc  <- induced_subgraph(g_top, nucleo)
louv_r_nuc <- cluster_louvain(g_top_nuc, weights = E(g_top_nuc)$r)
louv_m_nuc <- cluster_louvain(g_top_nuc, weights = E(g_top_nuc)$m_media)

# Reciprocidad ponderada (Squartini et al. 2013)
W_all <- as_adjacency_matrix(g_dir, attr = "weight", sparse = FALSE)
reciproc_pond <- sum(pmin(W_all, t(W_all))) / sum(W_all)

# Puentes
puentes <- as_ids(bridges(g_top))


# ============================================================
# AISLAMIENTO POR DISTANCIA — dos pruebas de Mantel
# ============================================================

dg <- as.dist(dist_geo[pop_order, pop_order])

# 1) Diferenciación: FST pareado vs distancia (IBD -> correlación positiva)
set.seed(123)
mantel_fst <- mantel.randtest(as.dist(fst_matrix[pop_order, pop_order]), dg,
                              nrepet = 9999, alter = "greater")

# 2) Red: correlación parcial vs distancia (IBD -> correlación negativa)
pc <- Pcorr[pop_order, pop_order]; diag(pc) <- 0
set.seed(123)
mantel_red <- mantel.randtest(as.dist(pc), dg, nrepet = 9999, alter = "less")

mantel_tabla <- data.frame(
  Prueba = c("FST vs distancia", "Correlación parcial vs distancia"),
  Hipotesis_alternativa = c("r > 0", "r < 0"),
  r = round(c(mantel_fst$obs, mantel_red$obs), 3),
  p = round(c(mantel_fst$pvalue, mantel_red$pvalue), 4))
print(mantel_tabla)


# ============================================================
# TABLAS Y EXPORTACIÓN
# ============================================================

metricas_nodo <- data.frame(
  Poblacion         = pop_order,
  Componente        = as.integer(comp$membership[pop_order]),
  Fuerza_saliente   = round(fuerza_sal[pop_order], 3),
  Fuerza_entrante   = round(fuerza_ent[pop_order], 3),
  Fuerza_total      = round(fuerza_total[pop_order], 3),
  rho               = round(rho[pop_order], 3),
  Rol               = rol[pop_order],
  PageRank_global   = round(pagerank[pop_order], 4),
  PageRank_nucleo   = round(pagerank_nucleo[pop_order], 4),   # NA fuera del núcleo
  pi_d_nucleo       = round(pi_markov[pop_order], 4),         # NA fuera del núcleo
  Betweenness       = betweenness_pond[pop_order],
  Armonica_entrante = round(harm_ent[pop_order], 3),
  Armonica_saliente = round(harm_sal[pop_order], 3),
  Coreness          = coreness_k[pop_order],
  Clustering_local  = round(trans_local[pop_order], 3),
  Comunidad_r       = membership(louv_r)[pop_order]
)
print(metricas_nodo)

metricas_red <- data.frame(
  Metrica = c("Nodos", "Aristas (no dirigidas)", "Componentes conexos",
              "Eficiencia global", "Transitividad global", "Reciprocidad ponderada",
              "Modularidad Q (red completa, r)", "Modularidad Q (red completa, m)",
              "Modularidad Q (núcleo, r)", "Modularidad Q (núcleo, m)", "Puentes"),
  Valor = c(vcount(g_top), ecount(g_top), comp$no,
            round(E0, 4), round(trans_global, 4), round(reciproc_pond, 4),
            round(modularity(louv_r), 4), round(modularity(louv_m), 4),
            round(modularity(louv_r_nuc), 4), round(modularity(louv_m_nuc), 4),
            paste(puentes, collapse = ", ")))
print(metricas_red)

if (!dir.exists("outputs")) dir.create("outputs")
write_xlsx(list(
  "Metricas_nodo" = metricas_nodo,
  "Metricas_red"  = metricas_red,
  "Mantel"        = mantel_tabla
), "outputs/metricas_red_10loci.xlsx")


# =============================================================================
# BLOQUE: Modelo dinámico de Markov
# (antes: modelo_dinámico.R)
# =============================================================================
rm(list = ls())   # memoria limpia para este bloque
# ============================================================
#  MODELO DINÁMICO — Ceroxylon quindiuense
#  Script 05 — Cadena de Markov + Laplaciano (tiempo continuo)
#
#  REQUIERE haber corrido antes: 03_red_conectividad_10loci.R
#  (debe existir: "edges_df_10loci.rds")
# ============================================================


# ── PASO 0 — PAQUETES ─────────────────────────────────────────

library(expm)
library(igraph)
library(openxlsx)


# ── PASO 1 — CONSTRUIR LA MATRIZ DE PESOS W (8x8) ─────────────

edges_df <- readRDS("edges_df_10loci.rds")
pop_order <- sort(unique(c(edges_df$from, edges_df$to)))
# Si alguna población quedó totalmente aislada (p. ej. CX) y no aparece
# en edges_df, agrégala manualmente para no perderla del espacio de estados:
pop_order <- sort(unique(c(pop_order, "AN","CO","CX","LI","PI","SF","TE","TO")))
n <- length(pop_order)

W <- matrix(0, nrow = n, ncol = n, dimnames = list(pop_order, pop_order))
for (i in seq_len(nrow(edges_df))) {
  W[edges_df$from[i], edges_df$to[i]] <- edges_df$weight[i]
}
cat("=== Matriz de pesos W ===\n")
print(round(W, 4))


# ============================================================
# PASO 2 — CADENA DE MARKOV: MATRIZ DE TRANSICIÓN
# ============================================================

P <- matrix(0, nrow = n, ncol = n, dimnames = list(pop_order, pop_order))
for (i in seq_len(n)) {
  suma_fila <- sum(W[i, ])
  if (suma_fila > 0) {
    P[i, ] <- W[i, ] / suma_fila
  } else {
    P[i, i] <- 1   # estado sin salidas -> autoabsorbente
  }
}
cat("\n=== Matriz de transición P ===\n")
print(round(P, 4))
cat("\nVerificación (cada fila debe sumar 1):\n")
print(round(rowSums(P), 4))


# ── PASO 3 — CLASES DE COMUNICACIÓN ───────────────────────────

g_markov <- graph_from_adjacency_matrix(P > 0, mode = "directed")
clases <- components(g_markov, mode = "strong")

cat("\n=== Clases de comunicación (fuertemente conexas) ===\n")
for (k in seq_len(clases$no)) {
  miembros <- pop_order[clases$membership == k]
  cat("Clase", k, ":", paste(miembros, collapse = ", "), "\n")
}

# Identificar cuáles son CERRADAS (ninguna transición sale de la clase)
clases_cerradas <- c()
for (k in seq_len(clases$no)) {
  miembros <- pop_order[clases$membership == k]
  suma_saliente <- sum(P[miembros, setdiff(pop_order, miembros), drop = FALSE])
  if (suma_saliente == 0) clases_cerradas <- c(clases_cerradas, k)
}
cat("\nClases CERRADAS:", paste(clases_cerradas, collapse = ", "), "\n")
cat("=> La cadena es", ifelse(length(clases_cerradas) > 1, "REDUCIBLE", "IRREDUCIBLE"), "\n")


# ── PASO 4 — DISTRIBUCIÓN ESTACIONARIA (por clase cerrada) ───

distribuciones_estacionarias <- list()

for (k in clases_cerradas) {
  miembros <- pop_order[clases$membership == k]
  if (length(miembros) == 1) {
    pi_k <- setNames(1, miembros)
  } else {
    P_sub <- P[miembros, miembros]
    ev <- eigen(t(P_sub))
    idx <- which.min(Mod(ev$values - 1))
    vec <- Re(ev$vectors[, idx])
    vec <- vec / sum(vec)
    pi_k <- setNames(vec, miembros)
  }
  distribuciones_estacionarias[[paste0("Clase_", k)]] <- pi_k
  cat("\nDistribución estacionaria — Clase", k, "(", paste(miembros, collapse=","), "):\n")
  print(round(sort(pi_k, decreasing = TRUE), 4))
}


# ============================================================
# PASO 5 — EXTENSIÓN A TIEMPO CONTINUO: LAPLACIANO
# ============================================================

D_out <- diag(rowSums(W))
Q <- W - D_out
rownames(Q) <- colnames(Q) <- pop_order

cat("\n=== Matriz generadora Q (filas deben sumar 0) ===\n")
print(round(rowSums(Q), 6))

# --- Difusión: p(t) para distintos t, partiendo de un pulso en CO (el nodo mas central) ---
nodo_inicial <- "CO"
p0 <- setNames(rep(0, n), pop_order)
p0[nodo_inicial] <- 1

tiempos <- c(0, 0.5, 1, 2, 5, 10, 20, 50, 100)
difusion <- sapply(tiempos, function(t) as.numeric(p0 %*% expm(t * Q)))
rownames(difusion) <- pop_order
colnames(difusion) <- paste0("t=", tiempos)
cat("\n=== Difusión p(t) partiendo de un pulso en", nodo_inicial, "===\n")
print(round(difusion, 4))


# --- Distribución estacionaria continua: vector propio izquierdo de Q, autovalor 0 ---
# (se calcula sobre la clase conectada; CX permanece aislada con masa 1 en si misma)
distribuciones_continuas <- list()

for (k in clases_cerradas) {
  miembros <- pop_order[clases$membership == k]
  if (length(miembros) == 1) {
    pi_continua_k <- setNames(1, miembros)
  } else {
    Q_sub <- Q[miembros, miembros]
    ev_q <- eigen(t(Q_sub))
    idx_q <- which.min(Mod(ev_q$values))
    vec <- Re(ev_q$vectors[, idx_q])
    vec <- vec / sum(vec)
    pi_continua_k <- setNames(vec, miembros)
  }
  distribuciones_continuas[[paste0("Clase_", k)]] <- pi_continua_k
  cat("\nDistribución estacionaria CONTINUA — Clase", k, "(", paste(miembros, collapse=","), "):\n")
  print(round(sort(pi_continua_k, decreasing = TRUE), 4))
}


# ── PASO 6 — COMPARACIÓN: DISCRETA vs. CONTINUA ───────────────

pi_discreta_conectada <- distribuciones_estacionarias[[
  names(which.max(sapply(distribuciones_estacionarias, length)))
]]

for (nombre_clase in names(distribuciones_estacionarias)) {
  pi_d <- distribuciones_estacionarias[[nombre_clase]]
  if (length(pi_d) > 1) {
    pi_c <- distribuciones_continuas[[nombre_clase]]
    comp <- data.frame(Poblacion = names(pi_d),
                       Discreta  = round(pi_d, 4),
                       Continua  = round(pi_c[names(pi_d)], 4))
    comp <- comp[order(-comp$Discreta), ]
    cat("\n=== Comparación discreta vs. continua —", nombre_clase, "===\n")
    print(comp)
  }
}


# ── PASO 7 — TASA DE SALIDA TOTAL POR NODO (explica diferencias) ──

tasa_salida <- data.frame(Poblacion = pop_order, Tasa_salida_total = round(rowSums(W), 4))
tasa_salida <- tasa_salida[order(tasa_salida$Tasa_salida_total), ]
cat("\n=== Tasa de salida total por nodo (menor = 'retiene' más tiempo la señal) ===\n")
print(tasa_salida)


# ── PASO 8 — EXPORTAR ──────────────────────────────────────────

wb <- createWorkbook()
addWorksheet(wb, "Matriz_W");        writeData(wb, "Matriz_W", as.data.frame(W), rowNames = TRUE)
addWorksheet(wb, "Matriz_P");        writeData(wb, "Matriz_P", as.data.frame(round(P,4)), rowNames = TRUE)
addWorksheet(wb, "Difusion_pt");     writeData(wb, "Difusion_pt", as.data.frame(round(difusion,4)), rowNames = TRUE)
addWorksheet(wb, "Comparacion_estacionarias"); writeData(wb, "Comparacion_estacionarias", comparacion)
addWorksheet(wb, "Tasa_salida");     writeData(wb, "Tasa_salida", tasa_salida)

if (!dir.exists("outputs")) dir.create("outputs")
saveWorkbook(wb, "outputs/modelo_dinamico_10loci.xlsx", overwrite = TRUE)
cat("\n✅ Excel guardado: outputs/modelo_dinamico_10loci.xlsx\n")

saveRDS(P, "P_markov_10loci.rds")
saveRDS(Q, "Q_laplaciano_10loci.rds")


# =============================================================================
# BLOQUE: Figuras de métricas y de Markov
# (antes: 18_figuras_metricas.R)
# =============================================================================
rm(list = ls())   # memoria limpia para este bloque
# =============================================================================
# 18_figuras_metricas.R
# Figuras de las métricas de red y del aislamiento por distancia.
#
#   1. metricas_perfil_nodos : seis métricas por población (paneles)
#   2. metricas_pagerank_pid : PageRank vs distribución estacionaria pi^(d)
#                              en el núcleo
#   3. mantel_ibd            : F_ST y correlación parcial vs distancia (28 pares)
#
# Entradas: red_10loci.rds (03a) y outputs/metricas_red_10loci.xlsx (métricas)
# Correr desde Tesis_Ceroxylon.
# =============================================================================

source("scripts/00_estilo_figuras.R")
library(readxl); library(dplyr); library(tidyr)

red  <- readRDS("red_10loci.rds")
nodo <- read_excel("outputs/metricas_red_10loci.xlsx", sheet = "Metricas_nodo")
mant <- read_excel("outputs/metricas_red_10loci.xlsx", sheet = "Mantel")

# Orden de las poblaciones: de mayor a menor fuerza total (arriba la mayor)
orden_pob <- nodo$Poblacion[order(nodo$Fuerza_total)]

# ---------------------------------------------------- 1. Perfil de nodos -----
perfil <- nodo %>%
  transmute(Poblacion,
            `Fuerza saliente`            = Fuerza_saliente,
            `Fuerza entrante`            = Fuerza_entrante,
            `PageRank`                   = PageRank_global,
            `Intermediación`             = Betweenness,
            `Cercanía armónica entrante` = Armonica_entrante,
            `Cercanía armónica saliente` = Armonica_saliente) %>%
  pivot_longer(-Poblacion, names_to = "Metrica", values_to = "valor") %>%
  mutate(Metrica   = factor(Metrica, levels = unique(Metrica)),
         Poblacion = factor(Poblacion, levels = orden_pob))

fig1 <- ggplot(perfil, aes(x = valor, y = Poblacion)) +
  geom_segment(aes(x = 0, xend = valor, yend = Poblacion), colour = "grey70", linewidth = 0.5) +
  geom_point(colour = AZUL, size = 2) +
  facet_wrap(~ Metrica, ncol = 3, scales = "free_x") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = NULL) +
  theme(panel.grid.major.y = element_blank(),
        panel.border = element_rect(fill = NA, colour = "grey55", linewidth = 0.4),
        panel.spacing.x = unit(0.5, "cm"),
        strip.text = element_text(size = TAM_PEQ))
guardar(fig1, "metricas_perfil_nodos", ANCHO_COMPLETO, 9)

# ------------------------- 2. Markov: visitas pi^(d) vs permanencia pi^(c) ---
# W: red dirigida podada (las 7 aristas en los dos sentidos), fila = origen
library(igraph)
g_dir <- graph_from_data_frame(red$edges_df, directed = TRUE, vertices = red$pop_order)
W <- as_adjacency_matrix(g_dir, attr = "weight", sparse = FALSE)
cat("\nTasas q_pq = w_pq (para el diagrama de tiempo continuo):\n"); print(round(W, 3))

clases <- list("Núcleo" = c("AN", "CO", "PI", "SF", "TE"), "LI–TO" = c("LI", "TO"))
pi_clase <- function(C) {
  Wc <- W[C, C]; d <- rowSums(Wc)
  ev <- eigen(t(Wc / d))                                 # P = D^-1 W
  pd <- Re(ev$vectors[, which.min(abs(ev$values - 1))]); pd <- pd / sum(pd)
  pc <- pd / d; pc <- pc / sum(pc)                       # pi^(c) ∝ pi^(d) / d_out
  data.frame(Poblacion = C, pi_d = pd, pi_c = pc)
}
markov <- bind_rows(lapply(names(clases), function(k)
  mutate(pi_clase(clases[[k]]), Clase = k)))
cat("\nDistribuciones estacionarias por clase:\n"); print(markov, digits = 3)

mk <- markov %>%
  pivot_longer(c(pi_d, pi_c), names_to = "Medida", values_to = "valor") %>%
  mutate(Medida = factor(Medida, levels = c("pi_d", "pi_c")),
         Clase  = factor(Clase, levels = names(clases)),
         Poblacion = factor(Poblacion, levels = markov$Poblacion[order(-markov$pi_d)]))

fig2 <- ggplot(mk, aes(x = Poblacion, y = valor, fill = Medida)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  facet_grid(~ Clase, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c(pi_d = AZUL, pi_c = "grey60"),
                    labels = c(pi_d = expression(pi^(d) ~ "(visitas)"),
                               pi_c = expression(pi^(c) ~ "(permanencia)")),
                    name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Probabilidad estacionaria") +
  theme(legend.position = "bottom", panel.grid.major.x = element_blank(),
        panel.border = element_rect(fill = NA, colour = "grey55", linewidth = 0.4))
guardar(fig2, "markov_visitas_permanencia", ANCHO_MEDIO, 7)

# --------------------------------------------- 3. Aislamiento por distancia --
po  <- red$pop_order
dg  <- red$dist_geo[po, po]
if (max(dg) > 1000) dg <- dg / 1000            # por si está en metros
fst <- red$fst_matrix[po, po]
pc  <- red$Pcorr[po, po]

pares <- t(combn(po, 2))
md <- data.frame(dist = dg[pares], fst = fst[pares], pc = pc[pares]) %>%
  pivot_longer(c(fst, pc), names_to = "panel", values_to = "valor") %>%
  mutate(panel = factor(panel, levels = c("fst", "pc"),
                        labels = c("F[ST]~pareado", "'Correlación parcial'")))

etiq <- data.frame(
  panel = factor(c("F[ST]~pareado", "'Correlación parcial'"),
                 levels = c("F[ST]~pareado", "'Correlación parcial'")),
  lab   = gsub("-", "−", sprintf("r = %.3f,  P = %.3f", mant$r, mant$p)))

fig3 <- ggplot(md, aes(x = dist, y = valor)) +
  geom_point(shape = 21, size = 2, fill = AZUL, colour = "white", stroke = 0.4) +
  geom_label(data = etiq, aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE,
             hjust = -0.05, vjust = 1.3, family = FUENTE, size = TAM_PEQ / .pt,
             colour = TINTA, fill = "white", label.size = 0.3,
             label.r = unit(0, "pt"), label.padding = unit(3, "pt")) +
  facet_wrap(~ panel, scales = "free_y", labeller = label_parsed) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.28))) +   # franja libre arriba
  labs(x = "Distancia geográfica (km)", y = NULL) +
  theme(panel.border = element_rect(fill = NA, colour = "grey55", linewidth = 0.4),
        panel.spacing.x = unit(0.5, "cm"))
guardar(fig3, "mantel_ibd", ANCHO_COMPLETO, 7)

