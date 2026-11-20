# =============================================================================
# 04_red.R
# Construcción de la red (Etapas 1 a 3), estabilidad de las aristas, roles y mapas de la red.
#
# Sección de la tesis: 3.3 Red de conectividad (construcción, estabilidad, direccionalidad)
# Lee:      datos_10loci_final.rds, fst_pareado_10loci_final.rds
# Produce:  red_10loci.rds, edges_df_10loci.rds, outputs/construccion_red_10loci.xlsx, outputs/13_*.csv, red_etapa1_mapa, red_etapa2_mapa, red_etapa3_mapa
# Tiempo:   el bloque de estabilidad tarda (1000 réplicas de divMigrate); para probar, B_mig <- 50
#
# Correr desde la carpeta Tesis_Ceroxylon, en una sesión limpia de R
# (Session > Restart R), después de los scripts con número menor.
# =============================================================================

source("scripts/00_configuracion.R")


# =============================================================================
# BLOQUE: Construcción de la red
# (antes: 03a construccion red 10loci.R)
# Cambio: set.seed antes de divMigrate.
# =============================================================================
# ============================================================
#  CONSTRUCCIÓN DE LA RED DE CONECTIVIDAD — Ceroxylon quindiuense
#  Script 03a — Etapas 1, 2 y 3 + roles de nodo (10 loci corregidos)
#
#  REQUIERE haber corrido antes: 02_estructura_genetica_10loci.R
#  (deben existir: "datos_10loci_final.rds" y "fst_pareado_10loci_final.rds")
#
#  GUARDA (para 03b y 06): "red_10loci.rds"  y  "edges_df_10loci.rds"
# ============================================================


# ── PASO 0 — PAQUETES ─────────────────────────────────────────


library(writexl); library(diveRsity); library(readxl); library(dplyr)
library(tidyr); library(MASS); library(igraph); library(geosphere); library(ggplot2)

if (!dir.exists("outputs")) dir.create("outputs")


# ── PASO 1 — CARGAR DATOS CORREGIDOS (10 loci) ───────────────

df <- readRDS("datos_10loci_final.rds")
loci <- c("Cs03","Cs24","Ca17","Ca07","Ca16","Cs05","Ca30","Ca20","Ca25","Cs30")
pop_order <- sort(unique(df$Localidad))   # "AN" "CO" "CX" "LI" "PI" "SF" "TE" "TO"
print(pop_order)
table(df$Localidad)


# ── PASO 2 — DE GENOTIPOS A FRECUENCIAS ALÉLICAS ─────────────

alelos_largo <- df %>%
  pivot_longer(cols = all_of(loci), names_to = "locus", values_to = "genotipo") %>%
  separate(genotipo, into = c("a1", "a2"), sep = "/") %>%
  pivot_longer(cols = c(a1, a2), names_to = "copia", values_to = "alelo") %>%
  filter(alelo != "0")

freq_tabla <- alelos_largo %>%
  group_by(Localidad, locus, alelo) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Localidad, locus) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()


# ── PASO 3 — MATRIZ POBLACIONES × ALELOS ─────────────────────

freq_wide <- freq_tabla %>%
  mutate(columna = paste(locus, alelo, sep = "_")) %>%
  dplyr::select(Localidad, columna, freq) %>%
  pivot_wider(names_from = columna, values_from = freq, values_fill = 0) %>%
  arrange(Localidad)

pop_names <- freq_wide$Localidad
freq_matrix <- as.matrix(freq_wide[, -1])
rownames(freq_matrix) <- pop_names
dim(freq_matrix)


# ============================================================
# ETAPA 1 — GRAFO COMPLETO DE SIMILITUD (1 - FST), descriptivo
# ============================================================

fst_matrix <- readRDS("fst_pareado_10loci_final.rds")
fst_matrix <- fst_matrix[pop_order, pop_order]
diag(fst_matrix) <- 0

w0_matrix <- 1 - fst_matrix
round(w0_matrix, 4)


# ============================================================
# ETAPA 2 — PODA POR INDEPENDENCIA CONDICIONAL
# ============================================================

# Covarianza genética entre poblaciones
Z <- scale(freq_matrix, center = TRUE, scale = FALSE)
m <- ncol(Z)
Sigma <- (Z %*% t(Z)) / (m - 1)
rownames(Sigma) <- colnames(Sigma) <- pop_names
eigen(Sigma)$values   # un autovalor ~ 0 (Sigma es singular por el centrado)

# Matriz de precisión (pseudoinversa) y correlación parcial
Omega <- ginv(Sigma)
d <- sqrt(diag(Omega))
Pcorr <- -Omega / outer(d, d)
diag(Pcorr) <- 1
rownames(Pcorr) <- colnames(Pcorr) <- pop_names
round(Pcorr, 4)

# Significancia: z de Fisher con corrección de Bonferroni
pares <- combn(pop_names, 2)
k <- length(pop_names) - 2
n <- m
resultados <- data.frame()
for (i in 1:ncol(pares)) {
  p1 <- pares[1, i]; p2 <- pares[2, i]
  r <- Pcorr[p1, p2]
  z <- 0.5 * log((1 + r) / (1 - r))
  se <- 1 / sqrt(n - k - 3)
  pval <- 2 * (1 - pnorm(abs(z / se)))
  resultados <- rbind(resultados, data.frame(Pop1 = p1, Pop2 = p2, r = r, p = pval))
}
alpha_bonf <- 0.05 / nrow(resultados)
resultados$significativa <- resultados$p < alpha_bonf
resultados <- resultados[order(resultados$p), ]
print(resultados)
cat("\nNúmero de aristas significativas (E_f):", sum(resultados$significativa), "\n")

aristas_finales <- resultados[resultados$significativa, ]
all(aristas_finales$r > 0)   # TRUE: requisito de Louvain en 03b


# ── COORDENADAS Y DISTANCIA GEODÉSICA (HAVERSINE) ─────────────

dms_a_decimal <- function(g, m, s) g + m/60 + s/3600
coords <- data.frame(
  Localidad = c("SF","CX","CO","TO","LI","PI","AN","TE"),
  lat = c(dms_a_decimal(5,15,18), dms_a_decimal(4,38,15), dms_a_decimal(4,38,15),
          dms_a_decimal(4,31,0),  dms_a_decimal(4,30,47), dms_a_decimal(4,20,2),
          dms_a_decimal(4,19,6),  dms_a_decimal(3,43,50)),
  lon = -c(dms_a_decimal(75,24,45), dms_a_decimal(75,29,18), dms_a_decimal(75,31,40),
           dms_a_decimal(75,28,30), dms_a_decimal(75,35,37), dms_a_decimal(75,42,35),
           dms_a_decimal(75,29,48), dms_a_decimal(76,5,20))
)
mat_coords <- as.matrix(coords[, c("lon","lat")])
rownames(mat_coords) <- coords$Localidad
dist_geo <- distm(mat_coords, fun = distHaversine) / 1000   # km
rownames(dist_geo) <- colnames(dist_geo) <- coords$Localidad
dist_geo <- dist_geo[pop_order, pop_order]


# ============================================================
# ETAPA 3 — DIRECCIONALIDAD (divMigrate, Sundqvist et al. 2016)
# ============================================================

genepop_geno <- function(genotipo) {
  ale <- strsplit(genotipo, "/")[[1]]
  ale_fmt <- ifelse(ale == "0", "000", formatC(as.numeric(ale), width = 3, format = "d", flag = "0"))
  paste0(ale_fmt[1], ale_fmt[2])
}
lineas <- c("Genepop file - Ceroxylon quindiuense - 10 loci corregidos", loci)
for (pop in pop_order) {
  lineas <- c(lineas, "Pop")
  sub <- df %>% filter(Localidad == pop)
  for (i in 1:nrow(sub)) {
    genos <- sapply(loci, function(l) genepop_geno(sub[[l]][i]))
    lineas <- c(lineas, paste0(sub$ID[i], " ,  ", paste(genos, collapse = " ")))
  }
}
writeLines(lineas, "ceroxylon_genepop_10loci.txt")

set.seed(2026)   # para que los intervalos bootstrap sean reproducibles
resultado_migracion <- divMigrate(
  infile = "ceroxylon_genepop_10loci.txt",
  outfile = "divmigrate_resultados_10loci",
  boots = 1000, stat = "gst", filter_threshold = 0
)
gmat <- resultado_migracion$gRelMig      # filas = origen, columnas = destino
dimnames(gmat) <- list(pop_order, pop_order)
gsig <- resultado_migracion$gRelMigSig   # 0 donde la asimetría NO es significativa
dimnames(gsig) <- list(pop_order, pop_order)
print(round(gmat, 3))

# Dirección de las aristas de E_f
direccion_final <- data.frame()
for (i in 1:nrow(aristas_finales)) {
  p1 <- aristas_finales$Pop1[i]; p2 <- aristas_finales$Pop2[i]
  direccion_final <- rbind(direccion_final, data.frame(
    desde = p1, hacia = p2, m_12 = gmat[p1, p2], m_21 = gmat[p2, p1],
    dominante = ifelse(gmat[p1, p2] > gmat[p2, p1], paste(p1, "->", p2), paste(p2, "->", p1)),
    sig_1_a_2 = gsig[p1, p2] > 0,
    sig_2_a_1 = gsig[p2, p1] > 0))
}
direccion_final$asimetria_significativa <- direccion_final$sig_1_a_2 | direccion_final$sig_2_a_1
print(direccion_final)

# Grafo dirigido final
edges_df <- data.frame()
for (i in 1:nrow(aristas_finales)) {
  p1 <- aristas_finales$Pop1[i]; p2 <- aristas_finales$Pop2[i]
  edges_df <- rbind(edges_df,
                    data.frame(from = p1, to = p2, weight = gmat[p1, p2]),
                    data.frame(from = p2, to = p1, weight = gmat[p2, p1]))
}
g_final <- graph_from_data_frame(edges_df, directed = TRUE, vertices = pop_order)


# ============================================================
# ROLES DE NODO (regla de significancia, Sección 4.7)
# ============================================================

fuerza_sal <- sapply(pop_order, function(p) sum(edges_df$weight[edges_df$from == p]))
fuerza_ent <- sapply(pop_order, function(p) sum(edges_df$weight[edges_df$to == p]))
names(fuerza_sal) <- names(fuerza_ent) <- pop_order

asim_sig <- mapply(function(a, b) gsig[a, b] > 0 | gsig[b, a] > 0, edges_df$from, edges_df$to)
pob_sig  <- sapply(pop_order, function(p) any(asim_sig[edges_df$from == p | edges_df$to == p]))

rho <- ifelse(fuerza_sal == 0 & fuerza_ent == 0, NA, fuerza_sal / fuerza_ent)
rol <- ifelse(is.na(rho), "Aislada",
              ifelse(!pob_sig, "Equilibrada",
                     ifelse(rho > 1, "Emisora", "Receptora")))
names(rol) <- pop_order

roles <- data.frame(Poblacion = pop_order, Fuerza_saliente = round(fuerza_sal, 3),
                    Fuerza_entrante = round(fuerza_ent, 3), rho = round(rho, 3), Rol = rol)
print(roles)


# ============================================================
# FIGURA — red dirigida sobre coordenadas reales
# ============================================================

col_rol <- c(Emisora = "#2ca25f", Receptora = "#de2d26",
             Equilibrada = "#fdae6b", Aislada = "#969696")
V(g_final)$color <- col_rol[rol[V(g_final)$name]]
layout_geo_plot <- as.matrix(coords[match(pop_order, coords$Localidad), c("lon","lat")])

edges_igraph <- as_ids(E(g_final))
from_nodes <- sapply(strsplit(edges_igraph, "\\|"), `[`, 1)
to_nodes   <- sapply(strsplit(edges_igraph, "\\|"), `[`, 2)
pesos_aristas <- edge_attr(g_final, "weight")
tipo_linea <- rep(1, length(pesos_aristas))
colores_arista <- rep("#08306b", length(pesos_aristas))
for (i in seq_along(edges_igraph)) {
  peso_inverso <- pesos_aristas[from_nodes == to_nodes[i] & to_nodes == from_nodes[i]]
  if (length(peso_inverso) > 0 && pesos_aristas[i] < peso_inverso) {
    tipo_linea[i] <- 2; colores_arista[i] <- "#9ecae1"
  }
}
plot(g_final, layout = layout_geo_plot,
     edge.curved = 0.35, edge.width = pesos_aristas * 6,
     edge.lty = tipo_linea, edge.color = colores_arista,
     edge.arrow.size = 0.8, edge.arrow.width = 1.0,
     vertex.size = 12, vertex.label.color = "black", vertex.label.font = 2,
     vertex.label.cex = 0.9, asp = 0.8,
     main = "Red de conectividad genética dirigida (10 loci)")
legend("bottomleft",
       legend = c("Emisora", "Receptora", "Equilibrada", "Aislada",
                  "Mayor migración", "Menor migración"),
       col = c(col_rol, "#08306b", "#9ecae1"),
       pch = c(15, 15, 15, 15, NA, NA), lty = c(NA, NA, NA, NA, 1, 2),
       lwd = c(NA, NA, NA, NA, 3, 3), bty = "n", cex = 0.8)


# ============================================================
# GUARDAR OBJETOS Y EXPORTAR
# ============================================================

saveRDS(edges_df, "edges_df_10loci.rds")   # lo usa el script 06
saveRDS(list(pop_order = pop_order, fst_matrix = fst_matrix, w0_matrix = w0_matrix,
             Pcorr = Pcorr, resultados = resultados, aristas_finales = aristas_finales,
             coords = coords, dist_geo = dist_geo, gmat = gmat, gsig = gsig,
             edges_df = edges_df, direccion_final = direccion_final,
             fuerza_sal = fuerza_sal, fuerza_ent = fuerza_ent, rho = rho, rol = rol),
        "red_10loci.rds")

w0_export <- cbind(Poblacion = rownames(w0_matrix), as.data.frame(round(w0_matrix, 4)))
pcorr_export <- cbind(Poblacion = rownames(Pcorr), as.data.frame(round(Pcorr, 4)))
write_xlsx(list(
  "Etapa1_1menosFST"      = w0_export,
  "Etapa2_Pcorr"          = pcorr_export,
  "Etapa2_28_pares"       = resultados,
  "Etapa2_significativas" = aristas_finales,
  "Etapa3_Direccion"      = direccion_final,
  "Roles"                 = roles
), "outputs/construccion_red_10loci.xlsx")


# =============================================================================
# BLOQUE: Estabilidad de las aristas y asimetría pareada
# (antes: 13_estabilidad_red.R)
# =============================================================================
rm(list = ls())   # memoria limpia para este bloque
# =============================================================================
# 13_estabilidad_red.R
# Responde a la revisión de la construcción de la red (Etapas 2 y 3).
#
#   A. Comprueba la identidad  Sigma^+ = (Sigma + J)^(-1) - J,  J = 11'/8,
#      que dice qué calcula la pseudo-inversa cuando el centrado por alelo
#      deja a Sigma con rango 7.
#   B. Sensibilidad de la Etapa 2:
#      B1. Fisher z con m - L observaciones (una restricción por locus).
#      B2. Covarianza regularizada Sigma_l = (1-l) Sigma + l diag(Sigma).
#      B3. Centrado convencional (por población) con inversa usual.
#   C. Estabilidad de las aristas por remuestreo (1000 réplicas):
#      C1. Bootstrap sobre loci.
#      C2. Bootstrap sobre individuos dentro de cada población.
#   D. Etapa 3: bootstrap pareado de Delta_pq = w_pq - w_qp con divMigrate.
#      divMigrate decide la significancia comparando por separado los
#      intervalos de cada sentido y no devuelve las réplicas, así que aquí se
#      remuestrea igual que divMigrate (individuos dentro de población) y en
#      cada réplica se estiman los dos sentidos sobre la misma muestra.
#
# Entradas: datos_10loci_final.rds
# Tiempo: C tarda pocos minutos; D tarda más (una llamada a divMigrate por
# réplica). Para una prueba rápida, poner B_mig <- 50.
# =============================================================================
library(dplyr); library(tidyr); library(MASS); library(diveRsity)
select <- dplyr::select
set.seed(2026)
dir.create("outputs", showWarnings = FALSE)

df   <- readRDS("datos_10loci_final.rds")
loci <- c("Cs03","Cs24","Ca17","Ca07","Ca16","Cs05","Ca30","Ca20","Ca25","Cs30")
pops <- sort(unique(df$Localidad))
L    <- length(loci)
B    <- 1000          # réplicas para C1 y C2
B_mig <- 1000         # réplicas para D

# ---------------------------------------------------------------- funciones --
frecuencias <- function(d) {
  d %>% select(Localidad, all_of(loci)) %>%
    pivot_longer(all_of(loci), names_to = "locus", values_to = "g") %>%
    filter(!is.na(g), g != "0/0") %>%
    separate(g, into = c("a1", "a2"), sep = "/") %>%
    pivot_longer(c(a1, a2), values_to = "alelo") %>%
    count(Localidad, locus, alelo) %>%
    group_by(Localidad, locus) %>% mutate(freq = n / sum(n)) %>% ungroup()
}
a_matriz <- function(fr) {
  fw <- fr %>% mutate(col = paste(locus, alelo, sep = "_")) %>%
    select(Localidad, col, freq) %>%
    pivot_wider(names_from = col, values_from = freq, values_fill = 0) %>%
    arrange(Localidad)
  X <- as.matrix(fw[, -1]); rownames(X) <- fw$Localidad; X
}
# Etapa 2 con opciones. Por defecto reproduce exactamente el script 03a.
etapa2 <- function(X, centrado = "alelo", lambda = 0, descuento = 0) {
  Z <- if (centrado == "alelo") scale(X, center = TRUE, scale = FALSE)
  else t(scale(t(X), center = TRUE, scale = FALSE))
  m <- ncol(Z)
  S <- (Z %*% t(Z)) / (m - 1)
  # Regularización que conserva la restricción del centrado (Sigma 1 = 0):
  # se encoge hacia mu (I - J) y no hacia diag(Sigma)
  if (lambda > 0) {
    n <- nrow(S); Jn <- matrix(1 / n, n, n); mu <- sum(diag(S)) / (n - 1)
    S <- (1 - lambda) * S + lambda * mu * (diag(n) - Jn)
  }
  O <- if (centrado == "alelo") ginv(S) else solve(S)
  d <- sqrt(diag(O)); R <- -O / outer(d, d)
  dimnames(R) <- list(rownames(X), rownames(X))
  pr <- combn(rownames(X), 2)
  r  <- R[cbind(pr[1, ], pr[2, ])]
  z  <- 0.5 * log((1 + r) / (1 - r)) * sqrt(m - descuento - 6 - 3)
  P  <- 2 * (1 - pnorm(abs(z)))
  data.frame(par = paste(pr[1, ], pr[2, ], sep = "-"), r = r, P = P,
             arista = P < 0.05 / 28)
}

fr0 <- frecuencias(df)
X0  <- a_matriz(fr0)
original <- etapa2(X0)
E_f <- original$par[original$arista]
cat("\nControl, aristas originales (deben ser 7):", length(E_f), "\n"); print(E_f)

# --------------------------------------------------------------------- A -----
Z0 <- scale(X0, center = TRUE, scale = FALSE)
S0 <- (Z0 %*% t(Z0)) / (ncol(Z0) - 1)
J  <- matrix(1 / 8, 8, 8)
cat("\n== A. Pseudo-inversa ==\n")
cat("Suma de cada columna de Z (debe ser 0):", signif(max(abs(colSums(Z0))), 3), "\n")
cat("Rango de Sigma:", qr(S0)$rank, " | Sigma %*% 1 (debe ser 0):",
    signif(max(abs(S0 %*% rep(1, 8))), 3), "\n")
cat("max |Sigma^+ - ((Sigma+J)^-1 - J)| :",
    signif(max(abs(ginv(S0) - (solve(S0 + J) - J))), 3), "\n")

# --------------------------------------------------------------------- B -----
# --------------------------------------------------------------------- B -----
variantes <- list(
  original        = etapa2(X0),
  fisher_m_menos_L = etapa2(X0, descuento = L),
  regul_0.05      = etapa2(X0, lambda = 0.05),
  regul_0.10      = etapa2(X0, lambda = 0.10),
  regul_0.20      = etapa2(X0, lambda = 0.20),
  centrado_pob    = etapa2(X0, centrado = "poblacion")
)
tabla_B <- bind_rows(lapply(names(variantes), function(n)
  variantes[[n]] %>% filter(arista) %>% transmute(variante = n, par, x = "x"))) %>%
  pivot_wider(names_from = variante, values_from = x, values_fill = "") %>%
  arrange(par)
cat("\n== B. Aristas en cada variante de la Etapa 2 ==\n"); print(tabla_B, n = Inf)
write.csv(tabla_B, "outputs/13_variantes_etapa2.csv", row.names = FALSE)

# Jackknife sobre loci: quitar un locus a la vez (sin loci duplicados)
jack <- sapply(loci, function(l) etapa2(a_matriz(fr0 %>% filter(locus != l)))$arista)
rownames(jack) <- original$par
tabla_jack <- data.frame(par = original$par, en_red_original = original$arista,
                         veces_de_10 = rowSums(jack)) %>%
  filter(en_red_original | veces_de_10 > 0) %>% arrange(desc(veces_de_10))
cat("\n== Jackknife sobre loci (en cuántas de las 10 versiones aparece cada arista) ==\n")
print(tabla_jack)
cat("\nLocus quitado y aristas que cambian:\n")
for (l in loci) {
  a <- original$par[jack[, l]]
  cat(sprintf("  sin %s: pierde {%s}  gana {%s}\n", l,
              paste(setdiff(E_f, a), collapse = ","), paste(setdiff(a, E_f), collapse = ",")))
}
write.csv(tabla_jack, "outputs/13_jackknife_loci.csv", row.names = FALSE)


# --------------------------------------------------------------------- C -----
cat("\n== C1. Bootstrap sobre loci ==\n")
res_loci <- sapply(1:B, function(b) {
  s <- sample(loci, L, replace = TRUE)
  frb <- bind_rows(lapply(seq_along(s), function(j)
    fr0 %>% filter(locus == s[j]) %>% mutate(locus = paste0(s[j], "_", j))))
  etapa2(a_matriz(frb))$arista
})
cat("== C2. Bootstrap sobre individuos dentro de población ==\n")
filas_pob <- split(seq_len(nrow(df)), df$Localidad)
res_ind <- sapply(1:B, function(b) {
  filas <- unlist(lapply(filas_pob, function(i) i[sample.int(length(i), length(i), replace = TRUE)]))
  etapa2(a_matriz(frecuencias(df[filas, ])))$arista
})
estab <- data.frame(par = original$par,
                    en_red_original = original$arista,
                    frec_boot_loci = rowMeans(res_loci),
                    frec_boot_individuos = rowMeans(res_ind)) %>%
  filter(en_red_original | frec_boot_loci > 0.05 | frec_boot_individuos > 0.05) %>%
  arrange(desc(en_red_original), desc(frec_boot_loci))
print(estab, digits = 3)
write.csv(estab, "outputs/13_estabilidad_aristas.csv", row.names = FALSE)

# --------------------------------------------------------------------- D -----
# Archivo genepop igual al del 03a
genepop_geno <- function(genotipo) {
  ale <- strsplit(genotipo, "/")[[1]]
  ale_fmt <- ifelse(ale == "0", "000", formatC(as.numeric(ale), width = 3, format = "d", flag = "0"))
  paste0(ale_fmt[1], ale_fmt[2])
}
escribir_genepop <- function(d, archivo) {
  lineas <- c("Genepop - bootstrap", loci)
  for (p in pops) {
    lineas <- c(lineas, "Pop")
    sub <- d[d$Localidad == p, ]
    for (i in seq_len(nrow(sub))) {
      genos <- sapply(loci, function(l) genepop_geno(sub[[l]][i]))
      lineas <- c(lineas, paste0(p, "_", i, " ,  ", paste(genos, collapse = " ")))
    }
  }
  writeLines(lineas, archivo)
}
migrar <- function(d) {
  f <- tempfile(fileext = ".txt"); escribir_genepop(d, f)
  g <- divMigrate(infile = f, outfile = NULL, boots = 0, stat = "gst",
                  filter_threshold = 0, plot_network = FALSE)$gRelMig
  unlink(f); dimnames(g) <- list(pops, pops); g
}

cat("\n== D. Bootstrap pareado de divMigrate (", B_mig, "réplicas) ==\n")
g0 <- migrar(df)
cat("Control: compara esta matriz con gRelMig del 03a\n"); print(round(g0, 3))
boot_g <- array(NA, c(8, 8, B_mig), dimnames = list(pops, pops, NULL))
for (b in 1:B_mig) {
  filas <- unlist(lapply(filas_pob, function(i) i[sample.int(length(i), length(i), replace = TRUE)]))
  boot_g[, , b] <- migrar(df[filas, ])
  if (b %% 100 == 0) cat("  réplica", b, "\n")
}
asim <- bind_rows(lapply(E_f, function(e) {
  pq <- strsplit(e, "-")[[1]]; p <- pq[1]; q <- pq[2]
  bp <- boot_g[p, q, ]; bq <- boot_g[q, p, ]; delta <- bp - bq
  ip <- quantile(bp, c(.025, .975)); iq <- quantile(bq, c(.025, .975))
  id <- quantile(delta, c(.025, .975))
  data.frame(arista = e, w_pq = g0[p, q], w_qp = g0[q, p],
             delta = g0[p, q] - g0[q, p],
             delta_IC_inf = id[1], delta_IC_sup = id[2],
             prop_delta_pos = mean(delta > 0),
             asim_pareada = id[1] > 0 | id[2] < 0,
             asim_no_solapamiento = ip[1] > iq[2] | iq[1] > ip[2])
}))
rownames(asim) <- NULL
print(asim, digits = 3)
write.csv(asim, "outputs/13_asimetria_pareada.csv", row.names = FALSE)


# =============================================================================
# BLOQUE: Roles de nodo
# (antes: nuevo)
# =============================================================================
rm(list = ls())   # memoria limpia para este bloque
# Roles de nodo con el bootstrap pareado (método de la tesis). El script de
# construcción asignaba los roles con la regla de no solapamiento de
# divMigrate; aquí se reemplazan por los del bootstrap pareado y se vuelve a
# guardar red_10loci.rds, que usan 05_metricas_markov.R y las figuras.
red  <- readRDS("red_10loci.rds")
asim <- read.csv("outputs/13_asimetria_pareada.csv")
par_sig <- asim$arista[asim$asim_pareada]
sig_pob <- sapply(red$pop_order, function(p)
  any(sapply(strsplit(par_sig, "-"), function(e) p %in% e)))
rol <- ifelse(is.na(red$rho), "Aislada",
              ifelse(!sig_pob, "Equilibrada",
                     ifelse(red$rho > 1, "Emisora", "Receptora")))
names(rol) <- red$pop_order
red$rol <- rol
saveRDS(red, "red_10loci.rds")
roles <- data.frame(Poblacion = red$pop_order, Fuerza_saliente = round(red$fuerza_sal, 3),
                    Fuerza_entrante = round(red$fuerza_ent, 3), rho = round(red$rho, 3), Rol = rol)
print(roles)
write.csv(roles, "outputs/roles_nodos.csv", row.names = FALSE)
# Nota: la hoja "Roles" de construccion_red_10loci.xlsx tiene la regla anterior;
# los roles de la tesis son los de outputs/roles_nodos.csv.


# =============================================================================
# BLOQUE: Mapas de la red (Etapas 1, 2 y 3)
# (antes: 17_mapas_red.R)
# =============================================================================
rm(list = ls())   # memoria limpia para este bloque
# =============================================================================
# 17_mapas_red.R
# Las tres etapas de la red sobre el mapa, con el estilo común.
#
#   Etapa 1: grafo completo, 28 aristas con color según 1 - F_ST
#   Etapa 2: red podada, 7 aristas con color según r
#   Etapa 3: red dirigida, cada arista con color según la migración media de
#            los dos sentidos; solo la asimetría significativa (bootstrap
#            pareado, script 13) lleva flecha
#
# Fondo: departamentos con los mismos colores del mapa del área de estudio.
# Entradas: red_10loci.rds (03a), outputs/13_asimetria_pareada.csv (13),
#           capas/ (GADM, descargado por el script 11)
# Correr desde Tesis_Ceroxylon.
# =============================================================================

source("scripts/00_estilo_figuras.R")
library(sf); library(geodata); library(ggspatial); library(ggrepel); library(dplyr)

red    <- readRDS("red_10loci.rds")
coords <- red$coords                       # Localidad, lat, lon
w0     <- red$w0_matrix                    # 1 - F_ST
ar     <- red$aristas_finales              # Pop1, Pop2, r
gmat   <- red$gmat                         # divMigrate, fila = origen

# Coordenadas por población (búsqueda por nombre, sin depender del orden)
lon_de <- setNames(coords$lon, coords$Localidad)
lat_de <- setNames(coords$lat, coords$Localidad)

# ---------------------------------------------------------------- base -------
col_dep <- st_as_sf(gadm("COL", level = 1, path = "capas"))
deps_estudio <- c("Quindío", "Tolima", "Caldas", "Valle del Cauca")
col_dep$depto <- ifelse(col_dep$NAME_1 %in% deps_estudio, col_dep$NAME_1, "Otros")
colores_dep <- c("Quindío" = "#a8d5c2", "Tolima" = "#d6e2f5", "Caldas" = "#f0d5ec",
                 "Valle del Cauca" = "#fbf1cf", "Otros" = "#ececec")

margen <- c(0.30, 0.20)
lim_x  <- range(coords$lon) + c(-1, 1) * margen[1]
lim_y  <- range(coords$lat) + c(-1, 1) * margen[2]

base_mapa <- function() {
  ggplot() +
    geom_sf(data = col_dep, aes(fill = depto), colour = "grey55", linewidth = 0.3,
            show.legend = FALSE) +
    scale_fill_manual(values = colores_dep) +
    coord_sf(xlim = lim_x, ylim = lim_y, expand = FALSE) +
    annotation_scale(location = "bl", width_hint = 0.25, text_family = FUENTE) +
    annotation_north_arrow(location = "tr", height = unit(0.9, "cm"), width = unit(0.7, "cm"),
                           style = north_arrow_orienteering(text_size = 7, text_family = FUENTE)) +
    labs(x = NULL, y = NULL) +
    theme(panel.grid.major = element_line(colour = "grey80", linewidth = 0.15),
          axis.text = element_text(size = TAM_PEQ - 1),
          legend.position = "bottom")
}
capa_nodos <- list(
  geom_point(data = coords, aes(x = lon, y = lat), shape = 21, size = 2.6,
             fill = "grey15", colour = "white", stroke = 0.7),
  geom_label_repel(data = coords, aes(x = lon, y = lat, label = Localidad),
                   family = FUENTE, size = TAM_PEQ / .pt, colour = TINTA,
                   fill = alpha("white", 0.85), label.size = 0, label.padding = 0.12,
                   box.padding = 0.4, point.padding = 0.3, min.segment.length = 0.2, seed = 1)
)
aristas_xy <- function(d) {
  d %>% mutate(x = lon_de[p], y = lat_de[p], xend = lon_de[q], yend = lat_de[q])
}
# Escala de calor sin el amarillo pálido (se pierde sobre el fondo del Valle)
escala_calor <- function(nombre) {
  scale_colour_gradientn(colours = RColorBrewer::brewer.pal(9, "YlOrRd")[4:9], name = nombre)
}
leyenda_barra <- guides(colour = guide_colourbar(barwidth = unit(4, "cm"),
                                                 barheight = unit(0.3, "cm")),
                        linewidth = "none")

# ---------------------------------------------------------------- Etapa 1 ----
pares <- t(combn(rownames(w0), 2))
e1 <- data.frame(p = pares[, 1], q = pares[, 2], valor = w0[pares]) %>%
  aristas_xy() %>% arrange(valor)                      # las fuertes quedan encima
fig1 <- base_mapa() +
  geom_segment(data = e1, aes(x = x, y = y, xend = xend, yend = yend,
                              colour = valor, linewidth = valor), lineend = "round") +
  escala_calor(expression(1 - F[ST])) +
  scale_linewidth(range = c(0.3, 1.3)) +
  capa_nodos + leyenda_barra
guardar(fig1, "red_etapa1_mapa", ANCHO_MEDIO, 14)

# ---------------------------------------------------------------- Etapa 2 ----
e2 <- ar %>% transmute(p = Pop1, q = Pop2, valor = r) %>% aristas_xy() %>% arrange(valor)
fig2 <- base_mapa() +
  geom_segment(data = e2, aes(x = x, y = y, xend = xend, yend = yend,
                              colour = valor, linewidth = valor), lineend = "round") +
  escala_calor(expression(r[pq])) +
  scale_linewidth(range = c(0.8, 1.8)) +
  capa_nodos + leyenda_barra
guardar(fig2, "red_etapa2_mapa", ANCHO_MEDIO, 14)

# ---------------------------------------------------------------- Etapa 3 ----
# Asimetría significativa con el bootstrap pareado (script 13)
sig <- if (file.exists("outputs/13_asimetria_pareada.csv")) {
  s <- read.csv("outputs/13_asimetria_pareada.csv") %>% filter(asim_pareada)
  pq <- strsplit(s$arista, "-")
  data.frame(p = ifelse(s$delta > 0, sapply(pq, `[`, 1), sapply(pq, `[`, 2)),
             q = ifelse(s$delta > 0, sapply(pq, `[`, 2), sapply(pq, `[`, 1)))
} else data.frame(p = "LI", q = "TO")

# Flecha curva al lado de la arista, de p hacia q, recortada en los extremos
flecha <- function(p, q, k = 0.25, n = 40, recorte = 0.15) {
  x1 <- lon_de[p]; y1 <- lat_de[p]; x2 <- lon_de[q]; y2 <- lat_de[q]
  cx <- (x1 + x2) / 2 + k * (y2 - y1); cy <- (y1 + y2) / 2 - k * (x2 - x1)
  t <- seq(recorte, 1 - recorte, length.out = n)
  data.frame(x = (1 - t)^2 * x1 + 2 * (1 - t) * t * cx + t^2 * x2,
             y = (1 - t)^2 * y1 + 2 * (1 - t) * t * cy + t^2 * y2,
             id = paste(p, q))
}
# Dos flechas por arista (una por sentido), cada una curvada hacia su lado
# y coloreada por su propia migración relativa w_pq (fila = origen)
dir3 <- bind_rows(lapply(seq_len(nrow(ar)), function(i) {
  a <- ar$Pop1[i]; b <- ar$Pop2[i]
  bind_rows(flecha(a, b, k = 0.15) %>% mutate(w = gmat[a, b], p = a, q = b),
            flecha(b, a, k = 0.15) %>% mutate(w = gmat[b, a], p = b, q = a))
})) %>%
  mutate(sig = ifelse(paste(p, q) %in% paste(sig$p, sig$q), "si", "no"))
orden <- dir3 %>% distinct(id, w) %>% arrange(w) %>% pull(id)   # fuertes encima
dir3$id <- factor(dir3$id, levels = orden)

fig3 <- base_mapa() +
  geom_path(data = dir3, aes(x = x, y = y, group = id, colour = w,
                             linewidth = w, linetype = sig),
            lineend = "round", arrow = arrow(length = unit(0.16, "cm"), type = "closed")) +
  escala_calor("Migración relativa") +
  scale_linewidth(range = c(0.4, 1.4)) +
  scale_linetype_manual(values = c(si = "solid", no = "22"), name = NULL,
                        labels = c(si = "Asimetría significativa",
                                   no = "Sin asimetría significativa")) +
  capa_nodos + leyenda_barra +
  guides(linetype = guide_legend(nrow = 2, override.aes = list(linewidth = 0.6)))
guardar(fig3, "red_etapa3_mapa", ANCHO_MEDIO, 14)
