# =============================================================================
# A2_agrupamiento.R
# Agrupamiento de alelos separados por 1 pb (abajo, arriba y aleatorio).
#
# Sección de la tesis: Anexo A (Agrupamiento de alelos); Tabla de estabilidad (columna de agrupamiento)
# Lee:      datos_10loci_final.rds
# Produce:  outputs/08_*.csv, outputs/09_*.csv
#
# Correr desde la carpeta Tesis_Ceroxylon, en una sesión limpia de R
# (Session > Restart R), después de los scripts con número menor.
# =============================================================================

source("scripts/00_configuracion.R")


# =============================================================================
# BLOQUE: Agrupamiento hacia abajo y hacia arriba
# (antes: 08 agrupamiento alelos sensibilidad.R)
# =============================================================================
# ============================================================
#  SENSIBILIDAD AL AGRUPAMIENTO DE ALELOS (binning) — Script 08
#  Responde al comentario A3 de Víctor Hugo (escenarios alternativos
#  de corrección) y a los alelos privados a 1 pb entre fuentes.
#
#  Idea: en la fuente CV casi todos los tamaños son pares, y en GR hay
#  pares e impares en proporciones parecidas. Si los loci son de
#  repetición de 2 pb, un mismo alelo real puede quedar como dos tamaños
#  vecinos (n y n + 1). Aquí se agrupan los tamaños en clases de 2 pb y se
#  repiten los análisis para ver qué cambia.
#
#  Dos escenarios, porque un tamaño impar queda justo entre dos clases:
#    "abajo":  impar -> par inmediatamente inferior   (2 * floor(t / 2))
#    "arriba": impar -> par inmediatamente superior   (2 * ceiling(t / 2))
#
#  REQUIERE: "datos_10loci_final.rds", "fst_pareado_10loci_final.rds"
#  GUARDA:   outputs/08_resumen_escenarios.csv, outputs/08_aristas_escenarios.csv
# ============================================================

library(dplyr); library(tidyr); library(adegenet); library(hierfstat); library(MASS)

df <- readRDS("datos_10loci_final.rds")
loci <- c("Cs03","Cs24","Ca17","Ca07","Ca16","Cs05","Ca30","Ca20","Ca25","Cs30")
pop_order <- sort(unique(df$Localidad))
fuente <- c(AN = "GR", CO = "GR", PI = "GR", SF = "GR", TE = "GR",
            CX = "CV", LI = "CV", TO = "CV")

# ── Funciones ────────────────────────────────────────────────

agrupar <- function(d, modo) {
  f <- if (modo == "abajo") function(t) 2 * floor(t / 2) else function(t) 2 * ceiling(t / 2)
  for (l in loci) {
    d[[l]] <- sapply(d[[l]], function(g) {
      if (is.na(g) || g == "0/0") return(g)
      a <- as.numeric(strsplit(g, "/")[[1]])
      paste(f(a), collapse = "/")
    })
  }
  d
}

a_genind <- function(d) {
  geno <- d[, loci]
  geno[] <- lapply(geno, function(x) ifelse(x == "0/0", NA, x))
  df2genind(geno, sep = "/", ind.names = d$ID, pop = d$Localidad,
            ploidy = 2, type = "codom", NA.char = NA)
}

# Frecuencias alélicas por población (como en el 03a)
frecuencias <- function(d) {
  d %>% dplyr::select(Localidad, all_of(loci)) %>%
    pivot_longer(all_of(loci), names_to = "locus", values_to = "g") %>%
    filter(g != "0/0") %>%
    separate(g, into = c("a1", "a2"), sep = "/") %>%
    pivot_longer(c(a1, a2), values_to = "alelo") %>%
    count(Localidad, locus, alelo) %>%
    group_by(Localidad, locus) %>% mutate(freq = n / sum(n)) %>% ungroup()
}

# Etapa 2 del 03a: correlación parcial con pseudo-inversa y poda de Bonferroni
aristas_poda <- function(fr) {
  fw <- fr %>% mutate(col = paste(locus, alelo, sep = "_")) %>%
    dplyr::select(Localidad, col, freq) %>%
    pivot_wider(names_from = col, values_from = freq, values_fill = 0) %>%
    arrange(Localidad)
  X <- as.matrix(fw[, -1]); rownames(X) <- fw$Localidad
  Z <- scale(X, center = TRUE, scale = FALSE)
  m <- ncol(Z)
  S <- (Z %*% t(Z)) / (m - 1)
  O <- ginv(S)
  d <- sqrt(diag(O))
  R <- -O / outer(d, d); dimnames(R) <- list(rownames(X), rownames(X))
  pares <- combn(rownames(X), 2)
  out <- data.frame(p = pares[1, ], q = pares[2, ])
  out$r <- mapply(function(a, b) R[a, b], out$p, out$q)
  out$z <- 0.5 * log((1 + out$r) / (1 - out$r)) * sqrt(m - 6 - 3)
  out$P <- 2 * (1 - pnorm(abs(out$z)))
  out$arista <- out$P < 0.05 / 28
  out$m <- m
  out
}

resumen_escenario <- function(d, nombre) {
  gen <- a_genind(d)
  fr  <- frecuencias(d)
  # alelos privados
  priv <- fr %>% group_by(locus, alelo) %>% filter(n_distinct(Localidad) == 1) %>%
    ungroup() %>% count(Localidad, name = "privados")
  # diversidad
  bs <- basic.stats(genind2hierfstat(gen))
  div <- data.frame(Localidad = colnames(bs$Ho),
                    Ho = round(colMeans(bs$Ho, na.rm = TRUE), 3),
                    He = round(colMeans(bs$Hs, na.rm = TRUE), 3),
                    Fis = round(colMeans(bs$Fis, na.rm = TRUE), 3))
  res <- div %>% left_join(priv, by = "Localidad") %>%
    mutate(escenario = nombre, alelos_totales = sum(gen@loc.n.all))
  # red
  ar <- aristas_poda(fr) %>% mutate(escenario = nombre)
  # F_ST pareado
  fst <- as.matrix(pairwise.WCfst(genind2hierfstat(gen)))
  list(res = res, aristas = ar, fst = fst)
}

# ── Escenarios ───────────────────────────────────────────────

esc <- list(
  original = resumen_escenario(df, "original"),
  abajo    = resumen_escenario(agrupar(df, "abajo"), "abajo"),
  arriba   = resumen_escenario(agrupar(df, "arriba"), "arriba")
)

resumen <- bind_rows(lapply(esc, `[[`, "res"))
aristas <- bind_rows(lapply(esc, `[[`, "aristas"))
write.csv(resumen, "outputs/08_resumen_escenarios.csv", row.names = FALSE)
write.csv(aristas, "outputs/08_aristas_escenarios.csv", row.names = FALSE)

cat("\nAlelos totales por escenario:\n")
print(resumen %>% distinct(escenario, alelos_totales))

cat("\nDiversidad y alelos privados por población y escenario:\n")
print(resumen %>% dplyr::select(escenario, Localidad, Ho, He, Fis, privados) %>%
        arrange(Localidad, escenario))

cat("\nAristas de la poda (P < 0.05/28) en cada escenario:\n")
print(aristas %>% filter(arista) %>%
        mutate(par = paste(p, q, sep = "-"), r = round(r, 3)) %>%
        dplyr::select(escenario, par, r) %>% arrange(par, escenario))

cat("\nCorrelación de los F_ST pareados con el escenario original:\n")
tri <- function(M) M[upper.tri(M)]
fst_orig <- tri(esc$original$fst)
cat("abajo :", round(cor(fst_orig, tri(esc$abajo$fst)), 3), "\n")
cat("arriba:", round(cor(fst_orig, tri(esc$arriba$fst)), 3), "\n")

# Control: la red original de este script debe dar las mismas 7 aristas del 03a
cat("\nControl, aristas del escenario original:",
    sum(esc$original$aristas$arista), "(deben ser 7)\n")


# =============================================================================
# BLOQUE: Agrupamiento aleatorio (500 réplicas)
# (antes: 09 agrupamiento aleatorio.R)
# =============================================================================
rm(list = ls())   # memoria limpia para este bloque
# ============================================================
#  INCERTIDUMBRE POR AGRUPAMIENTO DE ALELOS — Script 09
#
#  Todos los loci son de repetición de 2 pb (Chacón 2019, tabla A1 de la
#  tesis), así que los alelos reales tienen tamaños de una sola paridad.
#  Los tamaños impares no dicen si el alelo real era el par de abajo o el
#  de arriba. En vez de escoger un solo redondeo, aquí se repite muchas
#  veces un agrupamiento al azar: cada tamaño impar de cada locus se manda
#  a t - 1 o a t + 1 con probabilidad 1/2 (el mismo destino para todos los
#  individuos que tienen ese tamaño). En cada réplica se rehace la poda
#  de la Etapa 2 y se cuentan los alelos privados.
#
#  Resultado: con qué frecuencia aparece cada arista y cuántos alelos
#  privados tiene cada población, considerando la incertidumbre del
#  agrupamiento.
#
#  REQUIERE: "datos_10loci_final.rds"
#  GUARDA:   outputs/09_frecuencia_aristas.csv, outputs/09_privados.csv
# ============================================================

library(dplyr); library(tidyr); library(MASS)

df <- readRDS("datos_10loci_final.rds")
loci <- c("Cs03","Cs24","Ca17","Ca07","Ca16","Cs05","Ca30","Ca20","Ca25","Cs30")
pops <- df$Localidad

# Alelos como números (NA para los faltantes "0/0")
partir <- function(x, k) {
  v <- suppressWarnings(as.numeric(sapply(strsplit(x, "/"), `[`, k)))
  v[!is.na(v) & v == 0] <- NA
  v
}
A1 <- lapply(loci, function(l) partir(df[[l]], 1)); names(A1) <- loci
A2 <- lapply(loci, function(l) partir(df[[l]], 2)); names(A2) <- loci

impares <- lapply(loci, function(l) {
  t <- unique(c(A1[[l]], A2[[l]])); sort(t[!is.na(t) & t %% 2 == 1])
}); names(impares) <- loci
cat("Tamaños impares por locus:\n"); print(sapply(impares, length))

# Frecuencias alélicas por población a partir de los alelos numéricos
frecuencias <- function(B1, B2) {
  bind_rows(lapply(loci, function(l)
    data.frame(Localidad = c(pops, pops), locus = l, alelo = c(B1[[l]], B2[[l]])))) %>%
    filter(!is.na(alelo)) %>%
    count(Localidad, locus, alelo) %>%
    group_by(Localidad, locus) %>% mutate(freq = n / sum(n)) %>% ungroup()
}

# Etapa 2 del 03a (idéntica a la del script 08)
aristas_poda <- function(fr) {
  fw <- fr %>% mutate(col = paste(locus, alelo, sep = "_")) %>%
    dplyr::select(Localidad, col, freq) %>%
    pivot_wider(names_from = col, values_from = freq, values_fill = 0) %>%
    arrange(Localidad)
  X <- as.matrix(fw[, -1]); rownames(X) <- fw$Localidad
  Z <- scale(X, center = TRUE, scale = FALSE)
  m <- ncol(Z)
  O <- ginv((Z %*% t(Z)) / (m - 1))
  d <- sqrt(diag(O))
  R <- -O / outer(d, d); dimnames(R) <- list(rownames(X), rownames(X))
  pares <- combn(rownames(X), 2)
  r <- R[cbind(pares[1, ], pares[2, ])]
  z <- 0.5 * log((1 + r) / (1 - r)) * sqrt(m - 6 - 3)
  P <- 2 * (1 - pnorm(abs(z)))
  data.frame(par = paste(pares[1, ], pares[2, ], sep = "-"), r = r,
             arista = P < 0.05 / 28, m = m)
}

privados_por_pob <- function(fr) {
  fr %>% group_by(locus, alelo) %>% filter(n_distinct(Localidad) == 1) %>%
    ungroup() %>% count(Localidad, name = "privados") %>%
    right_join(data.frame(Localidad = sort(unique(pops))), by = "Localidad") %>%
    mutate(privados = ifelse(is.na(privados), 0, privados))
}

# ── Réplicas ─────────────────────────────────────────────────
set.seed(2026)
n_rep <- 500
res_aristas <- vector("list", n_rep)
res_priv <- vector("list", n_rep)
alelos_tot <- numeric(n_rep)

for (b in 1:n_rep) {
  B1 <- A1; B2 <- A2
  for (l in loci) {
    imp <- impares[[l]]
    if (length(imp) == 0) next
    destino <- imp + sample(c(-1, 1), length(imp), replace = TRUE)
    mapa <- setNames(destino, imp)
    i1 <- !is.na(B1[[l]]) & B1[[l]] %% 2 == 1
    i2 <- !is.na(B2[[l]]) & B2[[l]] %% 2 == 1
    B1[[l]][i1] <- mapa[as.character(B1[[l]][i1])]
    B2[[l]][i2] <- mapa[as.character(B2[[l]][i2])]
  }
  fr <- frecuencias(B1, B2)
  res_aristas[[b]] <- aristas_poda(fr) %>% mutate(rep = b)
  res_priv[[b]] <- privados_por_pob(fr) %>% mutate(rep = b)
  alelos_tot[b] <- nrow(distinct(fr, locus, alelo))
  if (b %% 50 == 0) cat("Réplica", b, "\n")
}

# ── Resúmenes ────────────────────────────────────────────────
aristas <- bind_rows(res_aristas)
frec_aristas <- aristas %>% group_by(par) %>%
  summarise(frecuencia = mean(arista), r_mediana = round(median(r), 3), .groups = "drop") %>%
  arrange(desc(frecuencia))
write.csv(frec_aristas, "outputs/09_frecuencia_aristas.csv", row.names = FALSE)
cat("\nFrecuencia con que aparece cada arista (solo las que aparecen al menos una vez):\n")
print(frec_aristas %>% filter(frecuencia > 0), n = 28)

cat("\nNúmero de aristas por réplica:\n")
print(table(aristas %>% group_by(rep) %>% summarise(k = sum(arista)) %>% pull(k)))

priv <- bind_rows(res_priv) %>% group_by(Localidad) %>%
  summarise(mediana = median(privados), p2.5 = quantile(privados, 0.025),
            p97.5 = quantile(privados, 0.975), .groups = "drop")
write.csv(priv, "outputs/09_privados.csv", row.names = FALSE)
cat("\nAlelos privados por población (mediana e intervalo del 95 % entre réplicas):\n")
print(priv)

cat("\nAlelos totales: mediana", median(alelos_tot), "; rango", range(alelos_tot), "\n")
