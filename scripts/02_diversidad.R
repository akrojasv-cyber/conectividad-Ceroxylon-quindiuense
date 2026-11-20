# =============================================================================
# 02_diversidad.R
# Diversidad genética por población y por locus, alelos nulos, alelos privados, F_IS y figuras.
#
# Sección de la tesis: 3.1 Diversidad genética; Anexo B (tabla por locus)
# Lee:      datos_10loci_final.rds
# Produce:  genind_10loci_final.rds, outputs/Resultados_Diversidad_Genetica_10loci_v2.xlsx, outputs/07_*.csv, figuras de diversidad
#
# Correr desde la carpeta Tesis_Ceroxylon, en una sesión limpia de R
# (Session > Restart R), después de los scripts con número menor.
# =============================================================================

source("scripts/00_configuracion.R")


# =============================================================================
# BLOQUE: Diversidad genética
# (antes: Diversidad_v2_final.R)
# =============================================================================
# ============================================================
#  DIVERSIDAD GENÉTICA — Ceroxylon quindiuense
#  Script 01 — VERSIÓN FINAL v2 (10 loci corregidos)
#
#  REQUIERE: "datos_10loci_final.rds" (de 00_estandarizacion.R)
#
#  CAMBIOS RESPECTO A LA VERSIÓN ANTERIOR (buscar "[v2]"):
#   1. El alelo "0" (sin amplificación) se trata como dato FALTANTE (NA),
#      no como un alelo más.
#   2. FST, FIS y FIT globales con Weir & Cockerham (wc) + IC 95% (boot.vc).
#   3. Significancia del FIS por población: IC 95% por bootstrap sobre loci
#      (boot.ppfis), en vez de permutar individuos entre poblaciones.
#   4. HWE: se agrega p ajustado por FDR (Benjamini-Hochberg).
#   5. Dunn: p bilateral directo (altp = TRUE).
#   6. Alelos privados y PIC: NA.method = "asis" (no imputar faltantes).
# ============================================================


# ── PASO 0 — PAQUETES ─────────────────────────────────────────


library(adegenet); library(hierfstat); library(poppr); library(mmod)
library(vegan); library(pegas); library(openxlsx); library(PopGenReport)
library(dunn.test); library(ggplot2)

sig_label <- function(p) {
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))
}


# ── PASO 1 — IMPORTAR DATOS Y CONSTRUIR genind ───────────────

df_10 <- readRDS("datos_10loci_final.rds")
loci_10 <- c("Cs03","Cs24","Ca17","Ca07","Ca16","Cs05","Ca30","Ca20","Ca25","Cs30")

# [v2] Cualquier genotipo con un alelo "0" se convierte en NA (dato faltante)
X <- df_10[, loci_10]
X[] <- lapply(X, function(col) {
  col <- as.character(col)
  col[grepl("(^|/)0($|/)", col)] <- NA
  col
})

# [v2] Reporte de datos faltantes (para Métodos)
cat("=== Genotipos faltantes ===\n")
cat("Total:", sum(is.na(X)), "de", length(as.matrix(X)),
    "(", round(100 * mean(is.na(X)), 2), "% )\n")
cat("Por población y locus:\n")
print(aggregate(is.na(X), by = list(Poblacion = df_10$Localidad), FUN = sum))

genind_obj <- df2genind(
  X,
  ploidy    = 2,
  ind.names = df_10$ID,
  pop       = as.factor(df_10$Localidad),
  sep       = "/",
  type      = "codom"
)

hf <- genind2hierfstat(genind_obj)
pops_names <- levels(pop(genind_obj))
loci_names <- locNames(genind_obj)
pops_sep   <- seppop(genind_obj)


# ── PASO 2 — ESTADÍSTICOS BÁSICOS ────────────────────────────

bs <- basic.stats(hf)
Ho <- bs$Ho
He <- bs$Hs                       # diversidad génica de Nei (1987), insesgada
Ne <- 1 / (1 - He)
Na <- nAll(genind_obj)
Na_locus_pop <- sapply(pops_sep, function(p) nAll(p))

cat("Ho por población:\n"); print(round(colMeans(Ho, na.rm = TRUE), 4))
cat("He por población:\n"); print(round(colMeans(He, na.rm = TRUE), 4))


# ── PASO 2b — [v2] F-ESTADÍSTICOS GLOBALES (Weir & Cockerham) ─

wc_global <- wc(hf)
FST_wc <- wc_global$FST
FIS_wc <- wc_global$FIS
FIT_wc <- 1 - (1 - FST_wc) * (1 - FIS_wc)

set.seed(123)
boot_global <- boot.vc(levels = data.frame(pop = hf[, 1]), loci = hf[, -1],
                       diploid = TRUE, nboot = 1000)
# Columnas de boot.vc: F-Pop/Total = FST ; F-Ind/Total = FIT ; F-Ind/Pop = FIS
cat("\n=== F-estadísticos globales (Weir & Cockerham) ===\n")
cat("FST =", round(FST_wc, 4), "| FIS =", round(FIS_wc, 4), "| FIT =", round(FIT_wc, 4), "\n")
cat("IC 95% (bootstrap sobre loci, 1000 réplicas):\n"); print(round(boot_global$ci, 4))

# Solo para comparar con la versión anterior (estimadores de Nei, NO usar en el artículo)
cat("\n=== basic.stats()$overall (Nei) — solo para comparar ===\n")
print(round(bs$overall, 4))

# F-estadísticos por locus (Tabla 3)
wc_locus <- wc_global$per.loc
Tabla_F_locus <- data.frame(
  Locus = loci_names,
  FIS = round(wc_locus$FIS, 4),
  FST = round(wc_locus$FST, 4),
  FIT = round(1 - (1 - wc_locus$FST) * (1 - wc_locus$FIS), 4)
)
print(Tabla_F_locus)

# Prueba de de Meeûs (2018): correlación entre FIS y FST a través de loci
# Nulos -> correlación positiva fuerte; Wahlund -> débil o nula
cor_meeus <- cor.test(wc_locus$FIS, wc_locus$FST, method = "spearman")
print(cor_meeus)
# ── PASO 3 — RIQUEZA ALÉLICA Y PIC ───────────────────────────

Ar <- allelic.richness(hf)
cat("\n[v2] Tamaño de rarefacción de Ar (número mínimo de alelos):",
    unique(Ar$min.all), "\n")

calc_PIC <- function(obj_genind) {
  freq_matrix <- tab(obj_genind, freq = TRUE, NA.method = "asis")   # [v2]
  loci_factor <- locFac(obj_genind)
  nombres_loci <- locNames(obj_genind)
  pic_values <- numeric(length(nombres_loci)); names(pic_values) <- nombres_loci
  for (i in seq_along(nombres_loci)) {
    cols_locus <- which(loci_factor == nombres_loci[i])
    p <- colMeans(freq_matrix[, cols_locus], na.rm = TRUE)
    p <- p[p > 0]
    pic_values[i] <- 1 - sum(p^2) - sum(outer(p^2, p^2)[lower.tri(outer(p, p))])
  }
  return(pic_values)
}
PIC <- calc_PIC(genind_obj)


# ── PASO 4 — [v2] FIS POR POBLACIÓN CON IC 95% (bootstrap sobre loci) ──
# Mantiene a cada individuo en su población y evalúa si FIS > 0.
# FIS significativo si el IC 95% NO incluye el cero.

FIS_obs_pop <- colMeans(bs$Fis, na.rm = TRUE)

set.seed(123)
boot_fis <- boot.ppfis(hf, nboot = 1000)
ci_fis <- boot_fis$fis.ci

Tabla_FIS_poblacion <- data.frame(
  Poblacion = pops_names,
  FIS       = round(FIS_obs_pop, 4),
  IC_inf    = round(ci_fis$ll, 4),
  IC_sup    = round(ci_fis$hl, 4),
  Sig       = ifelse(ci_fis$ll > 0 | ci_fis$hl < 0, "*", "ns")
)
cat("\n=== FIS por población con IC 95% ===\n"); print(Tabla_FIS_poblacion)

# FIS por locus: solo descriptivo (su respaldo estadístico son las pruebas HWE)
FIS_obs_locus <- bs$perloc[, "Fis"]
Tabla_FIS_locus <- data.frame(Locus = loci_names, FIS = round(FIS_obs_locus, 4))
print(Tabla_FIS_locus)


# ── PASO 5 — HWE DETALLADO ───────────────────────────────────

set.seed(123)
hwe_detail <- do.call(rbind, lapply(pops_names, function(p) {
  hwe_p <- hw.test(pops_sep[[p]], B = 1000)
  data.frame(Poblacion = p, Locus = rownames(hwe_p),
             Chi2 = round(hwe_p[, "chi^2"], 4), df = hwe_p[, "df"],
             p_exacto = round(hwe_p[, "Pr.exact"], 4))
}))
# TO-Ca30 es todo dato faltante y su prueba sale NA; se excluye (79 pruebas)
hwe_detail <- hwe_detail[!is.na(hwe_detail$p_exacto), ]

# [v2] Corrección por comparaciones múltiples (79 pruebas)
hwe_detail$p_FDR  <- round(p.adjust(hwe_detail$p_exacto, method = "BH"), 4)
hwe_detail$p_FDR  <- round(p.adjust(hwe_detail$p_exacto, method = "BH"), 4)
hwe_detail$Sig     <- sig_label(hwe_detail$p_exacto)
hwe_detail$Sig_FDR <- sig_label(hwe_detail$p_FDR)
hwe_detail$En_HWE     <- ifelse(hwe_detail$p_exacto >= 0.05, "SI", "NO")
hwe_detail$En_HWE_FDR <- ifelse(hwe_detail$p_FDR    >= 0.05, "SI", "NO")

HWE_resumen_poblacion <- do.call(rbind, lapply(pops_names, function(p) {
  sub_p <- hwe_detail[hwe_detail$Poblacion == p, ]
  data.frame(Poblacion = p, Loci_total = nrow(sub_p),
             Loci_en_HWE = sum(sub_p$En_HWE == "SI"),
             Pct_en_HWE = round(mean(sub_p$En_HWE == "SI") * 100, 1),
             Loci_en_HWE_FDR = sum(sub_p$En_HWE_FDR == "SI"),
             Pct_en_HWE_FDR = round(mean(sub_p$En_HWE_FDR == "SI") * 100, 1))
}))
HWE_resumen_locus <- do.call(rbind, lapply(loci_names, function(l) {
  sub_l <- hwe_detail[hwe_detail$Locus == l, ]
  data.frame(Locus = l, Pops_fuera_HWE = sum(sub_l$En_HWE == "NO"),
             Pops_fuera_HWE_FDR = sum(sub_l$En_HWE_FDR == "NO"),
             Pops_total = nrow(sub_l))
}))
print(HWE_resumen_poblacion); print(HWE_resumen_locus)


# ── PASO 6 — ALELOS NULOS ────────────────────────────────────
# [v2] Se guardan los dos estimadores; revisar ?null.all para saber
# a qué método corresponde cada resumen y citar ese en Métodos.

alelos_nulos <- PopGenReport::null.all(genind_obj)
tabla_alelos_nulos   <- as.data.frame(alelos_nulos$null.allele.freq$summary1)
tabla_alelos_nulos_2 <- as.data.frame(alelos_nulos$null.allele.freq$summary2)
cat("\n=== Alelos nulos: summary1 ===\n"); print(round(tabla_alelos_nulos, 3))
cat("\n=== Alelos nulos: summary2 ===\n"); print(round(tabla_alelos_nulos_2, 3))


# ── PASO 7 — ALELOS PRIVADOS ──────────────────────────────────

freq_list <- lapply(seppop(genind_obj), function(x) {
  colMeans(tab(x, freq = TRUE, NA.method = "asis"), na.rm = TRUE)   # [v2]
})
all_alleles <- names(freq_list[[1]])
allele_presence <- sapply(all_alleles, function(al) {
  sum(sapply(freq_list, function(f) f[al] > 0), na.rm = TRUE)
})
private_alleles <- names(allele_presence[allele_presence == 1])
private_by_pop <- sapply(pops_names, function(p) {
  sum(sapply(private_alleles, function(al) freq_list[[p]][al] > 0), na.rm = TRUE)
})
Tabla_Privados <- data.frame(
  Poblacion = pops_names, N = as.numeric(table(pop(genind_obj))),
  Alelos_privados = as.numeric(private_by_pop),
  Pct_del_total = round(as.numeric(private_by_pop) / length(private_alleles) * 100, 1)
)
print(Tabla_Privados)

# [v2] Verificación: ya no debe aparecer ningún alelo ".0"
cat("Alelos '.0' entre los privados (debe ser 0):",
    sum(grepl("\\.0$", private_alleles)), "\n")


# ── PASO 8 — % POLIMORFISMO ──────────────────────────────────

loci_polimorficos <- sapply(seppop(genind_obj), function(x) sum(isPoly(x)))
Tabla_Polimorfismo <- data.frame(
  Poblacion = names(loci_polimorficos), Loci_Polimorficos = as.numeric(loci_polimorficos),
  Porcentaje_P = round((as.numeric(loci_polimorficos) / length(loci_10)) * 100, 1)
)


# ── PASO 9 — CURVA DE ACUMULACIÓN DE LOCI ────────────────────

set.seed(123)
matriz_curva <- poppr::genotype_curve(genind_obj, sample = 1000, quiet = TRUE, plot = FALSE)
df_curva <- data.frame(
  Numero_Loci = 1:ncol(matriz_curva),
  Media_Genotipos = colMeans(matriz_curva),
  Desviacion_Est = apply(matriz_curva, 2, sd)
)


# ── PASO 10 — COMPARACIÓN DE He ENTRE POBLACIONES ────────────

df_stats_dunn <- stack(as.data.frame(He))
colnames(df_stats_dunn) <- c("He", "Poblacion")
prueba_kruskal <- kruskal.test(He ~ Poblacion, data = df_stats_dunn)
print(prueba_kruskal)

df_dunn_excel <- NULL
if (prueba_kruskal$p.value < 0.05) {
  dunn_res <- dunn.test::dunn.test(df_stats_dunn$He, df_stats_dunn$Poblacion,
                                   method = "bonferroni", altp = TRUE)   # [v2] p bilateral
  df_dunn_excel <- data.frame(Comparacion = dunn_res$comparisons, Z = dunn_res$Z,
                              P_bilateral_ajustado = dunn_res$altP.adjusted)
}


# ── PASO 11 — TABLAS RESUMEN ──────────────────────────────────

Tabla_Locus_Global <- data.frame(
  Locus = loci_names, Na = Na, Ho = round(rowMeans(Ho, na.rm = TRUE), 4),
  He = round(rowMeans(He, na.rm = TRUE), 4), Ne = round(rowMeans(Ne, na.rm = TRUE), 4),
  PIC = round(PIC, 4), FIS = round(FIS_obs_locus, 4),
  Nulos = round(as.numeric(tabla_alelos_nulos[1, loci_names]), 3)
)

Tabla_Poblaciones <- data.frame(
  Poblacion = pops_names, N = as.numeric(table(pop(genind_obj))),
  Ho = round(colMeans(Ho, na.rm = TRUE), 4), He = round(colMeans(He, na.rm = TRUE), 4),
  Ne = round(colMeans(Ne, na.rm = TRUE), 4), Ar = round(colMeans(Ar$Ar, na.rm = TRUE), 4),
  FIS = round(FIS_obs_pop, 4), FIS_IC_inf = round(ci_fis$ll, 4), FIS_IC_sup = round(ci_fis$hl, 4),
  HWE_pct_en_HWE = HWE_resumen_poblacion$Pct_en_HWE[match(pops_names, HWE_resumen_poblacion$Poblacion)],
  Alelos_privados = as.numeric(private_by_pop),
  Pct_Polimorfismo = Tabla_Polimorfismo$Porcentaje_P[match(pops_names, Tabla_Polimorfismo$Poblacion)]
)
print(Tabla_Poblaciones)

Tabla_Locus_Poblacion <- data.frame(
  Poblacion = rep(colnames(Ho), each = nrow(Ho)), Locus = rep(rownames(Ho), times = ncol(Ho)),
  Na = as.vector(Na_locus_pop), Ho = as.vector(Ho), He = as.vector(He),
  Ne = as.vector(Ne), FIS = as.vector((He - Ho) / He)
)

Resumen_Global <- data.frame(
  Parametro = c("Genotipos faltantes (%)", "Ho global", "He global", "Na total",
                "Na promedio por locus", "Ne promedio", "Ar promedio",
                "Tamano de rarefaccion (alelos)", "PIC medio",
                "FST (W&C)", "FIS (W&C)", "FIT (W&C)",
                "Total alelos privados", "Poblaciones con alelos privados"),
  Valor = round(c(100 * mean(is.na(X)),
                  mean(colMeans(Ho, na.rm = TRUE)), mean(colMeans(He, na.rm = TRUE)),
                  sum(Na), mean(Na), mean(colMeans(Ne, na.rm = TRUE)),
                  mean(colMeans(Ar$Ar, na.rm = TRUE)), unique(Ar$min.all)[1], mean(PIC),
                  FST_wc, FIS_wc, FIT_wc,
                  length(private_alleles), sum(private_by_pop > 0)), 4)
)
print(Resumen_Global)


# ── PASO 12 — EXPORTAR TODO A EXCEL ──────────────────────────

wb <- createWorkbook()
addWorksheet(wb, "Locus_Global");          writeData(wb, "Locus_Global", Tabla_Locus_Global)
addWorksheet(wb, "Resumen_Poblaciones");   writeData(wb, "Resumen_Poblaciones", Tabla_Poblaciones)
addWorksheet(wb, "Resumen_Global");        writeData(wb, "Resumen_Global", Resumen_Global)
addWorksheet(wb, "F_globales_IC");         writeData(wb, "F_globales_IC", as.data.frame(boot_global$ci), rowNames = TRUE)
addWorksheet(wb, "F_por_locus");           writeData(wb, "F_por_locus", Tabla_F_locus)
addWorksheet(wb, "Locus_x_Poblacion");     writeData(wb, "Locus_x_Poblacion", Tabla_Locus_Poblacion)
addWorksheet(wb, "FIS_poblacion_IC");      writeData(wb, "FIS_poblacion_IC", Tabla_FIS_poblacion)
addWorksheet(wb, "HWE_detalle");           writeData(wb, "HWE_detalle", hwe_detail)
addWorksheet(wb, "HWE_resumen_poblacion"); writeData(wb, "HWE_resumen_poblacion", HWE_resumen_poblacion)
addWorksheet(wb, "HWE_resumen_locus");     writeData(wb, "HWE_resumen_locus", HWE_resumen_locus)
addWorksheet(wb, "Alelos_Nulos_summary1"); writeData(wb, "Alelos_Nulos_summary1", tabla_alelos_nulos, rowNames = TRUE)
addWorksheet(wb, "Alelos_Nulos_summary2"); writeData(wb, "Alelos_Nulos_summary2", tabla_alelos_nulos_2, rowNames = TRUE)
addWorksheet(wb, "Alelos_privados");       writeData(wb, "Alelos_privados", Tabla_Privados)
addWorksheet(wb, "Polimorfismo");          writeData(wb, "Polimorfismo", Tabla_Polimorfismo)
addWorksheet(wb, "Curva_Acumulacion");     writeData(wb, "Curva_Acumulacion", df_curva)
addWorksheet(wb, "Kruskal_He");            writeData(wb, "Kruskal_He", data.frame(
  Estadistico = prueba_kruskal$statistic,
  df = prueba_kruskal$parameter, p_valor = prueba_kruskal$p.value))
if (!is.null(df_dunn_excel)) { addWorksheet(wb, "Dunn_post_hoc"); writeData(wb, "Dunn_post_hoc", df_dunn_excel) }

if (!dir.exists("outputs")) dir.create("outputs")
saveWorkbook(wb, "outputs/Resultados_Diversidad_Genetica_10loci_v2.xlsx", overwrite = TRUE)
cat("\n✅ Excel guardado: outputs/Resultados_Diversidad_Genetica_10loci_v2.xlsx\n")

# El script 02 lee este archivo: ahora contendrá el genind con faltantes como NA
saveRDS(genind_obj, "genind_10loci_final.rds")


# =============================================================================
# BLOQUE: Verificación de alelos nulos, privados y F_IS
# (antes: 07 verificacion nulos privados fis.R)
# =============================================================================
rm(list = ls())   # memoria limpia para este bloque
# ============================================================
#  VERIFICACIÓN DE ALELOS NULOS, ALELOS PRIVADOS Y F_IS
#  Script 07 — responde a los comentarios de Víctor Hugo (A4, A5, A6)
#
#  REQUIERE: "datos_10loci_final.rds" (el mismo que usa el 03a)
#  GUARDA en outputs/:
#    07_nulos_global.csv, 07_nulos_por_poblacion.csv,
#    07_faltantes_vs_esperados.csv, 07_privados_detalle.csv,
#    07_privados_rarefaccion.csv, 07_fis_comparacion.csv
# ============================================================


# ── PASO 0 — PAQUETES ─────────────────────────────────────────

library(dplyr); library(tidyr); library(adegenet); library(hierfstat); library(PopGenReport)

if (!dir.exists("outputs")) dir.create("outputs")


# ── PASO 1 — DATOS ───────────────────────────────────────────

df <- readRDS("datos_10loci_final.rds")
loci <- c("Cs03","Cs24","Ca17","Ca07","Ca16","Cs05","Ca30","Ca20","Ca25","Cs30")
pop_order <- sort(unique(df$Localidad))
fuente <- c(AN = "GR", CO = "GR", PI = "GR", SF = "GR", TE = "GR",
            CX = "CV", LI = "CV", TO = "CV")

# Revisar cómo están codificados los faltantes antes de seguir
cat("\nValores que contienen '0' en los genotipos:\n")
print(table(unlist(df[loci])[grepl("(^0/)|(/0$)", unlist(df[loci]))]))

# Genotipos con algún alelo "0" se tratan como faltantes
geno <- df[, loci]
geno[] <- lapply(geno, function(x) ifelse(grepl("(^0/)|(/0$)", x), NA, x))

gen <- df2genind(geno, sep = "/", ind.names = df$ID, pop = df$Localidad,
                 ploidy = 2, type = "codom", NA.char = NA)
gen


# ============================================================
# PARTE A — ALELOS NULOS
# ============================================================

# A1. Estimación con todos los individuos juntos (como en la tesis)
nul_global <- null.all(gen)
str(nul_global$null.allele.freq, max.level = 1)
print(nul_global$null.allele.freq)
saveRDS(nul_global, "outputs/07_nulos_global.rds")

# A2. La misma estimación dentro de cada población.
#     Con todas las poblaciones juntas, el déficit de heterocigotos que
#     produce la estructura (efecto Wahlund) se atribuye a alelos nulos.
#     Dentro de cada población ese componente desaparece.
nul_pob <- lapply(pop_order, function(p) {
  sub <- gen[pop(gen) == p]
  r <- tryCatch(null.all(sub)$null.allele.freq, error = function(e) NULL)
  if (is.null(r)) return(NULL)
  # summary1 = Chakraborty et al. (1992); summary2 = Brookfield (1996)
  data.frame(Localidad = p,
             locus = colnames(r$summary1),
             chakraborty = as.numeric(r$summary1[1, ]),
             brookfield1 = as.numeric(r$summary2[1, ]))
})
nul_pob <- bind_rows(nul_pob)
write.csv(nul_pob, "outputs/07_nulos_por_poblacion.csv", row.names = FALSE)

cat("\nNulos (Chakraborty) dentro de cada población, promedio por locus:\n")
print(nul_pob %>% group_by(locus) %>%
        summarise(media = round(mean(chakraborty, na.rm = TRUE), 3),
                  min = round(min(chakraborty, na.rm = TRUE), 3),
                  max = round(max(chakraborty, na.rm = TRUE), 3)))

# A3. Prueba de consistencia: si un alelo nulo tiene frecuencia r, en
#     equilibrio una fracción r^2 de los individuos es homocigota para el
#     nulo y no amplifica (aparece como faltante). Se compara esa
#     expectativa con los faltantes observados.
falt <- df %>%
  dplyr::select(Localidad, all_of(loci)) %>%
  pivot_longer(all_of(loci), names_to = "locus", values_to = "g") %>%
  mutate(faltante = grepl("(^0/)|(/0$)", g) | is.na(g)) %>%
  group_by(locus) %>%
  summarise(n = n(), faltantes = sum(faltante),
            prop_faltante = round(mean(faltante), 3), .groups = "drop")

r_global <- as.numeric(nul_global$null.allele.freq$summary1[1, ])
names(r_global) <- colnames(nul_global$null.allele.freq$summary1)
falt$r_chakraborty <- round(r_global[falt$locus], 3)
falt$faltantes_esperados <- round(falt$n * falt$r_chakraborty^2, 1)
write.csv(falt, "outputs/07_faltantes_vs_esperados.csv", row.names = FALSE)
cat("\nFaltantes observados frente a los esperados si los nulos fueran reales:\n")
print(falt)


# ============================================================
# PARTE B — ALELOS PRIVADOS
# ============================================================

# Conteo de copias de cada alelo por población (sin faltantes)
copias <- df %>%
  dplyr::select(Localidad, all_of(loci)) %>%
  pivot_longer(all_of(loci), names_to = "locus", values_to = "g") %>%
  filter(!grepl("(^0/)|(/0$)", g), !is.na(g)) %>%
  separate(g, into = c("a1", "a2"), sep = "/") %>%
  pivot_longer(c(a1, a2), values_to = "alelo") %>%
  count(locus, alelo, Localidad, name = "n_copias")

n_por_pob <- copias %>% group_by(locus, Localidad) %>%
  summarise(N = sum(n_copias), .groups = "drop")

# B1. Lista de alelos privados con su frecuencia dentro de la población
privados <- copias %>%
  group_by(locus, alelo) %>% filter(n() == 1) %>% ungroup() %>%
  left_join(n_por_pob, by = c("locus", "Localidad")) %>%
  mutate(frecuencia = round(n_copias / N, 3),
         fuente = fuente[Localidad],
         una_sola_copia = n_copias == 1)

# B2. Distancia (pb) al alelo más cercano de la OTRA fuente en el mismo
#     locus. Un privado a 1 o 2 pb de un alelo de la otra fuente podría
#     ser el mismo alelo mal corregido.
tamanos_fuente <- copias %>% mutate(fuente = fuente[Localidad], tam = as.numeric(alelo)) %>%
  distinct(locus, fuente, tam)
privados$dist_otra_fuente <- mapply(function(l, a, f) {
  otros <- tamanos_fuente$tam[tamanos_fuente$locus == l & tamanos_fuente$fuente != f]
  if (length(otros) == 0) NA else min(abs(otros - as.numeric(a)))
}, privados$locus, privados$alelo, privados$fuente)
write.csv(privados, "outputs/07_privados_detalle.csv", row.names = FALSE)

cat("\nAlelos privados por población (comparar con la tabla: total 117, CX 31, TE 26):\n")
print(privados %>% group_by(Localidad, fuente) %>%
        summarise(privados = n(),
                  una_copia = sum(una_sola_copia),
                  frec_menor_0.05 = sum(frecuencia < 0.05),
                  a_1o2pb_otra_fuente = sum(dist_otra_fuente <= 2, na.rm = TRUE),
                  .groups = "drop"))
cat("Total:", nrow(privados), "\n")

# B3. Riqueza de alelos privados por rarefacción (Kalinowski 2004).
#     Para un tamaño de submuestra g, la probabilidad de que el alelo i NO
#     aparezca en g copias tomadas de la población j es
#        Q_ij = C(N_j - N_ij, g) / C(N_j, g),
#     y la riqueza privada de j es la suma, sobre alelos y loci, de
#        (1 - Q_ij) * prod_{k != j} Q_ik.
g <- min(n_por_pob$N)
cat("\nTamaño de rarefacción g (mínimo de copias en un locus y población):", g, "\n")

Q <- function(Nj, Nij, g) ifelse(Nj - Nij < g, 0, exp(lchoose(Nj - Nij, g) - lchoose(Nj, g)))

tabla_q <- copias %>%
  dplyr::select(locus, alelo, Localidad, n_copias) %>%
  complete(nesting(locus, alelo), Localidad, fill = list(n_copias = 0)) %>%
  left_join(n_por_pob, by = c("locus", "Localidad")) %>%
  mutate(Qij = Q(N, n_copias, g))

rare_priv <- tabla_q %>%
  group_by(locus, alelo) %>%
  mutate(prod_todos = prod(Qij)) %>%
  ungroup() %>%
  mutate(aporte = (1 - Qij) * ifelse(Qij > 0, prod_todos / Qij, 0))
# Cuando Qij = 0 (el alelo aparece seguro en j) el producto de los demás
# se calcula aparte para no dividir por cero
idx0 <- which(rare_priv$Qij == 0)
for (i in idx0) {
  otros <- tabla_q %>% filter(locus == rare_priv$locus[i], alelo == rare_priv$alelo[i],
                              Localidad != rare_priv$Localidad[i])
  rare_priv$aporte[i] <- prod(otros$Qij)
}

priv_rare <- rare_priv %>% group_by(Localidad) %>%
  summarise(riqueza_privada_total = round(sum(aporte), 2),
            riqueza_privada_por_locus = round(sum(aporte) / length(loci), 3),
            .groups = "drop") %>%
  mutate(fuente = fuente[Localidad]) %>%
  left_join(privados %>% count(Localidad, name = "privados_crudos"), by = "Localidad")
write.csv(priv_rare, "outputs/07_privados_rarefaccion.csv", row.names = FALSE)
cat("\nRiqueza de alelos privados por rarefacción:\n")
print(priv_rare)

privados %>%
  mutate(dist = ifelse(dist_otra_fuente == 1, "1 pb", ifelse(dist_otra_fuente == 2, "2 pb", ">2 pb"))) %>%
  count(Localidad, dist) %>%
  tidyr::pivot_wider(names_from = dist, values_from = n, values_fill = 0)

tam <- copias %>% mutate(fuente = fuente[Localidad], tam = as.numeric(alelo)) %>%
  distinct(locus, fuente, tam)

# 1. Paridad de los tamaños por locus y fuente: si una fuente queda casi toda
#    en pares y la otra en impares, hay un desfase de 1 pb en la corrección
paridad <- tam %>% mutate(par = ifelse(tam %% 2 == 0, "par", "impar")) %>%
  count(locus, fuente, par) %>%
  tidyr::pivot_wider(names_from = c(fuente, par), values_from = n, values_fill = 0)
print(paridad)

# 2. Pares de alelos a 1 pb DENTRO de la misma fuente frente a ENTRE fuentes
vecinos <- tam %>% inner_join(tam, by = "locus", suffix = c("_a", "_b"),
                              relationship = "many-to-many") %>%
  filter(tam_b - tam_a == 1) %>%
  mutate(tipo = ifelse(fuente_a == fuente_b, "misma fuente", "entre fuentes")) %>%
  count(locus, tipo) %>%
  tidyr::pivot_wider(names_from = tipo, values_from = n, values_fill = 0)
print(vecinos)

# 3. ¿En qué loci se concentran los privados a 1 pb?
print(privados %>% filter(dist_otra_fuente == 1) %>% count(locus, sort = TRUE))

# 1. Paridad por POBLACIÓN: si CX se comporta distinto de LI y TO,
#    es probable que venga de otro genotipado (Sanín 2013)
copias %>% mutate(tam = as.numeric(alelo)) %>%
  distinct(Localidad, locus, tam) %>%
  mutate(par = ifelse(tam %% 2 == 0, "par", "impar")) %>%
  count(Localidad, par) %>%
  tidyr::pivot_wider(names_from = par, values_from = n, values_fill = 0) %>%
  mutate(prop_impar = round(impar / (par + impar), 2))

# 2. Número de alelos de la fuente GR en tu matriz. El artículo reporta
#    237 alelos en 126 individuos. Si da 237, tu matriz es la misma que se publicó.
copias %>% filter(fuente[Localidad] == "GR") %>% distinct(locus, alelo) %>% nrow()
# ============================================================
# PARTE C — ¿QUÉ ESTIMADOR DE F_IS TIENE LA TABLA?
# ============================================================

hf <- genind2hierfstat(gen)
bs <- basic.stats(hf)

fis_comp <- data.frame(
  Localidad = colnames(bs$Ho),
  Ho = round(colMeans(bs$Ho, na.rm = TRUE), 3),
  Hs = round(colMeans(bs$Hs, na.rm = TRUE), 3),
  # (a) la fórmula de la metodología con los promedios: 1 - Ho/He
  fis_formula = round(1 - colMeans(bs$Ho, na.rm = TRUE) / colMeans(bs$Hs, na.rm = TRUE), 3),
  # (b) promedio de los F_IS por locus que da basic.stats
  fis_promedio_loci = round(colMeans(bs$Fis, na.rm = TRUE), 3)
)
# (c) intervalos de boot.ppfis (Weir y Cockerham por población, bootstrap sobre loci)
set.seed(123)
bp <- boot.ppfis(hf, nboot = 1000)
fis_comp$ic_inf <- round(bp$fis.ci$ll, 3)
fis_comp$ic_sup <- round(bp$fis.ci$hl, 3)
fis_comp$tabla_tesis <- c(AN = 0.622, CO = 0.712, CX = 0.364, LI = 0.000,
                          PI = 0.587, SF = 0.520, TE = 0.778, TO = 0.066)[fis_comp$Localidad]
write.csv(fis_comp, "outputs/07_fis_comparacion.csv", row.names = FALSE)
cat("\nComparación de estimadores de F_IS con la tabla de la tesis:\n")
print(fis_comp)

library(dplyr); library(tidyr)
df   <- readRDS("datos_10loci_final.rds")
loci <- c("Cs03","Cs24","Ca17","Ca07","Ca16","Cs05","Ca30","Ca20","Ca25","Cs30")
fr <- df %>% select(Localidad, all_of(loci)) %>%
  pivot_longer(all_of(loci), names_to = "locus", values_to = "g") %>%
  filter(g != "0/0") %>% separate(g, c("a1", "a2"), sep = "/") %>%
  pivot_longer(c(a1, a2), values_to = "alelo") %>%
  count(Localidad, locus, alelo) %>%
  group_by(Localidad, locus) %>% mutate(freq = n / sum(n)) %>% ungroup()
privados <- fr %>% group_by(locus, alelo) %>% filter(n_distinct(Localidad) == 1) %>% ungroup()
resumen_priv <- privados %>% group_by(Localidad) %>%
  summarise(privados = n(),
            una_copia = sum(n == 1),
            frec_menor_0.05 = sum(freq < 0.05),
            frec_media = round(mean(freq), 3),
            frec_max = round(max(freq), 3))
print(resumen_priv)
cat("Total:", nrow(privados), "| una sola copia:", sum(privados$n == 1),
    "| frecuencia < 0.05:", sum(privados$freq < 0.05), "\n")
privados %>% filter(freq > 0.3)
# para el locus que salga, por ejemplo "Ca30":
fr %>% filter(locus == "Ca30") %>% arrange(as.numeric(alelo), Localidad) %>% print(n = Inf)

rangos <- fr %>% mutate(alelo = as.numeric(alelo)) %>%
  group_by(locus, Localidad) %>%
  summarise(min = min(alelo), max = max(alelo),
            moda = alelo[which.max(n)], n_alelos = n(), .groups = "drop")

# Para cada población y locus: ¿cuántas de sus copias caen fuera del rango
# que tienen las otras siete poblaciones en ese locus?
fuera <- fr %>% mutate(alelo = as.numeric(alelo)) %>%
  group_by(locus) %>%
  group_modify(function(d, k) {
    bind_rows(lapply(unique(d$Localidad), function(p) {
      otros <- d$alelo[d$Localidad != p]
      mio <- d[d$Localidad == p, ]
      fuera_rango <- mio$alelo < min(otros) | mio$alelo > max(otros)
      data.frame(Localidad = p,
                 prop_fuera = round(sum(mio$n[fuera_rango]) / sum(mio$n), 3))
    }))
  }) %>% ungroup()

print(fuera %>% filter(prop_fuera > 0.2) %>% arrange(desc(prop_fuera)), n = Inf)

# Heterocigosis fija: proporción de heterocigotos por población y locus
het <- df %>% select(Localidad, all_of(loci)) %>%
  pivot_longer(all_of(loci), names_to = "locus", values_to = "g") %>%
  filter(g != "0/0") %>% separate(g, c("a1", "a2"), sep = "/") %>%
  group_by(Localidad, locus) %>%
  summarise(prop_het = mean(a1 != a2),
            genotipo_mas_comun = names(which.max(table(paste(pmin(a1, a2), pmax(a1, a2), sep = "/")))),
            .groups = "drop")
print(het %>% filter(prop_het > 0.9), n = Inf)


# =============================================================================
# BLOQUE: Figuras de diversidad
# (antes: 14_figuras_diversidad.R)
# =============================================================================
rm(list = ls())   # memoria limpia para este bloque
# =============================================================================
# 14_figuras_diversidad.R
# Dos figuras para la subsección "Diversidad genética".
#
#   Figura 1. (a) Ho y He por población, unidas por una línea cuyo largo es el
#             déficit de heterocigotos. (b) F_IS con su IC del 95 %.
#             Valores tomados de la Tabla tab:diversidad-pob.
#   Figura 2. Individuos sin amplificación por locus, observados (en la
#             matriz) frente a esperados si los nulos fueran reales (N r^2).
#
# Colores: azul = González-Rivillas (GR), naranja = Chacón-Vargas (CV).
# Paleta validada para daltonismo (ΔE 24.7 en protanopía).
# =============================================================================

library(dplyr); library(tidyr); library(patchwork)
source("scripts/00_estilo_figuras.R")
col_fuente <- COL_FUENTE
tinta <- TINTA
tema <- theme(panel.grid.major.y = element_blank())   # solo ajustes propios

# --- Figura 1 -----------------------------------------------------------------
div <- data.frame(
  pob    = c("AN","CO","CX","LI","PI","SF","TE","TO"),
  fuente = c("GR","GR","CV","CV","GR","GR","GR","CV"),
  Ho  = c(0.310, 0.238, 0.471, 0.590, 0.330, 0.390, 0.150, 0.616),
  He  = c(0.863, 0.857, 0.697, 0.606, 0.838, 0.824, 0.793, 0.642),
  Fis = c(0.622, 0.712, 0.364, 0.000, 0.587, 0.520, 0.778, 0.066),
  inf = c(0.457, 0.561, 0.149, -0.233, 0.392, 0.368, 0.683, -0.250),
  sup = c(0.806, 0.856, 0.485, 0.261, 0.793, 0.699, 0.902, 0.310)
) %>%
  mutate(fuente = factor(fuente, c("GR", "CV"))) %>%
  arrange(desc(fuente), Fis) %>%             # abajo CV, arriba GR; dentro, por F_IS
  mutate(pob = factor(pob, levels = pob))

largo <- div %>% pivot_longer(c(Ho, He), names_to = "medida", values_to = "valor")

pa <- ggplot(div, aes(y = pob)) +
  geom_segment(aes(x = Ho, xend = He, yend = pob, colour = fuente), linewidth = 0.8) +
  geom_point(data = largo, aes(x = valor, shape = medida, colour = fuente),
             size = 2.6, fill = "white", stroke = 0.9) +
  scale_shape_manual(values = c(He = 16, Ho = 21), breaks = c("He", "Ho"), name = NULL,
                     labels = c(expression(H[e]), expression(H[o]))) +
  scale_colour_manual(values = col_fuente, name = "Fuente", labels = ETIQ_FUENTE) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "Heterocigosidad", y = NULL, tag = "(a)") + tema

pb <- ggplot(div, aes(y = pob, x = Fis, colour = fuente)) +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.4, linetype = "dashed") +
  geom_linerange(aes(xmin = inf, xmax = sup), linewidth = 0.8) +
  geom_point(size = 2.6) +
  scale_colour_manual(values = col_fuente, guide = "none") +
  scale_x_continuous(limits = c(-0.3, 1), breaks = seq(-0.2, 1, 0.2)) +
  labs(x = expression(F[IS]~"(IC 95 %)"), y = NULL, tag = "(b)") + tema +
  theme(axis.text.y = element_blank())

fig1 <- (pa | pb) + plot_layout(guides = "collect", widths = c(1.15, 1)) &
  theme(legend.position = "bottom", legend.box = "vertical")
guardar(fig1, "diversidad_poblaciones", ANCHO_COMPLETO, 8.5)

# --- Figura 2 -----------------------------------------------------------------
df   <- readRDS("datos_10loci_final.rds")
nulo <- c(Cs03 = 0.414, Cs24 = 0.174, Ca17 = 0.523, Ca07 = 0.270, Ca16 = 0.280,
          Cs05 = 0.471, Ca30 = 0.459, Ca20 = 0.524, Ca25 = 0.602, Cs30 = 0.357)
N <- nrow(df)
amp <- data.frame(locus = names(nulo), r = unname(nulo)) %>%
  mutate(esperados = N * r^2,
         observados = sapply(locus, function(l) sum(df[[l]] == "0/0" | is.na(df[[l]]))),
         locus = factor(locus, levels = locus[order(esperados)]))
print(amp)   # control: esperados 6-77, observados 0-4

fig2 <- ggplot(amp, aes(y = locus)) +
  geom_segment(aes(x = observados, xend = esperados, yend = locus),
               colour = "grey70", linewidth = 0.8) +
  geom_point(aes(x = esperados, shape = "Esperados si los nulos fueran reales"),
             size = 2.6, colour = AZUL) +
  geom_point(aes(x = observados, shape = "Observados en la matriz"),
             size = 2.6, colour = AZUL, fill = "white", stroke = 0.9) +
  scale_shape_manual(values = c("Esperados si los nulos fueran reales" = 16,
                                "Observados en la matriz" = 21), name = NULL) +
  geom_text(aes(x = esperados, label = round(esperados)), hjust = -0.6) +
  scale_x_continuous(limits = c(0, 85), breaks = seq(0, 80, 20)) +
  labs(x = "Individuos sin amplificación", y = NULL) + tema
guardar(fig2, "nulos_esperados_observados", ANCHO_MEDIO, 7)

# --- Figura 3. Curva de acumulación de genotipos ------------------------------
# Para k = 1, ..., 9 loci, poppr toma 1000 veces k loci al azar y cuenta cuántos
# genotipos multilocus distintos se distinguen. Si la curva se aplana antes de
# llegar al total de loci, el panel basta para separar a todos los individuos.
library(adegenet); library(poppr)
loci <- c("Cs03","Cs24","Ca17","Ca07","Ca16","Cs05","Ca30","Ca20","Ca25","Cs30")
geno <- df[, loci]; geno[] <- lapply(geno, function(x) ifelse(x == "0/0", NA, x))
gen  <- df2genind(geno, sep = "/", ind.names = df$ID, pop = df$Localidad,
                  ploidy = 2, type = "codom", NA.char = NA)
set.seed(2026)
gc <- genotype_curve(gen, sample = 1000, quiet = TRUE, plot = FALSE)
curva <- data.frame(loci = as.integer(colnames(gc)),
                    mediana = apply(gc, 2, median),
                    inf = apply(gc, 2, quantile, 0.025),
                    sup = apply(gc, 2, quantile, 0.975))
total_mlg <- mlg(gen, quiet = TRUE)
print(curva); cat("Genotipos multilocus con los 10 loci:", total_mlg, "\n")

fig3 <- ggplot(curva, aes(x = loci)) +
  geom_hline(yintercept = total_mlg, colour = "grey50", linewidth = 0.4, linetype = "dashed") +
  annotate("text", x = 1, y = total_mlg, label = paste(total_mlg, "genotipos"),
           hjust = 0, vjust = -0.5, family = FUENTE, size = TAM_PEQ / .pt, colour = tinta) +
  geom_ribbon(aes(ymin = inf, ymax = sup), fill = AZUL, alpha = 0.18) +
  geom_line(aes(y = mediana), colour = AZUL, linewidth = 0.8) +
  geom_point(aes(y = mediana), colour = AZUL, size = 2.2) +
  scale_x_continuous(breaks = curva$loci) +
  scale_y_continuous(limits = c(0, total_mlg * 1.08)) +
  labs(x = "Número de loci", y = "Genotipos multilocus distintos") + tema +
  theme(panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3))
guardar(fig3, "curva_genotipos", ANCHO_MEDIO, 7)
