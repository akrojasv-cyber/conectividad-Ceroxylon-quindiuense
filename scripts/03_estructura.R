# =============================================================================
# 03_estructura.R
# Diferenciación (F_ST, D_est), DAPC, k-means, sNMF, AMOVA y mapa de calor.
#
# Sección de la tesis: 3.2 Estructura genética
# Lee:      genind_10loci_final.rds
# Produce:  fst_pareado_10loci_final.rds, dest_pareado_10loci_final.rds, outputs/Resultados_Estructura_Genetica_10loci_v2.xlsx, Figura_DAPC, Figura_snmf, mapa_calor_fst_dest
# Tiempo:   las permutaciones pareadas tardan 15-30 min
#
# El Mantel de F_ST contra distancia está en 05_metricas_markov.R (distancia haversine).
#
# Correr desde la carpeta Tesis_Ceroxylon, en una sesión limpia de R
# (Session > Restart R), después de los scripts con número menor.
# =============================================================================

source("scripts/00_configuracion.R")


# =============================================================================
# BLOQUE: Estructura genética
# (antes: estructura_v2.R)
# Cambios: se quitó un bloque viejo (grupos de sNMF K = 2), las corridas repetidas
# de la AMOVA y la hoja "Mantel", que tomaba un objeto que este script no calcula.
# =============================================================================
# ============================================================
#  ESTRUCTURA GENÉTICA — Ceroxylon quindiuense
#  Script 02 — VERSIÓN FINAL v2 (10 loci corregidos)
#
#  REQUIERE haber corrido antes: 01_diversidad_genetica_10loci_v2.R
#  (debe existir "genind_10loci_final.rds" con los faltantes como NA)
#
#  CAMBIOS RESPECTO A LA VERSIÓN ANTERIOR (buscar "[v2]"):
#   1. Verifica que el genind ya no tenga el alelo "0".
#   2. Dest global y por locus (D_Jost) con IC 95% por bootstrap.
#   3. Permutaciones pareadas para FST Y para Dest, con corrección FDR.
#   4. DAPC no interactivo: número de PCs fijado con criterios explícitos.
#   5. snmf: faltantes = 9 (NA.method = "asis") y semilla fija.
#   6. Tablas de entropía cruzada, AMOVA y DAPC exportadas a Excel.
#   7. Figuras guardadas en TIFF a 600 dpi (requisito de Biotropica: >300 dpi).
#
#  OJO: el PASO 2 (permutaciones de FST y Dest) puede tardar 15–30 min.
# ============================================================


# ── PASO 0 — PAQUETES ─────────────────────────────────────────

library(adegenet); library(hierfstat); library(poppr); library(mmod)
library(LEA); library(openxlsx); library(ade4); library(geosphere)
library(ggplot2); library(ggrepel); library(dplyr); library(tidyr); library(patchwork)

if (!dir.exists("outputs")) dir.create("outputs")

Q <- LEA::Q
cross.entropy <- LEA::cross.entropy
# ── PASO 1 — IMPORTAR EL genind CORREGIDO ─────────────────────

genind_10 <- readRDS("genind_10loci_final.rds")
hf_10 <- genind2hierfstat(genind_10)
pops_names <- levels(pop(genind_10))

# [v2] Verificación: no debe existir ningún alelo "0"
alelos_cero <- grep("\\.0$", colnames(tab(genind_10)), value = TRUE)
if (length(alelos_cero) > 0) {
  stop("El genind todavía tiene alelos '0': ", paste(alelos_cero, collapse = ", "),
       "\nCorre primero 01_diversidad_genetica_10loci_v2.R")
} else cat("✅ Sin alelos '0': el genind está corregido.\n")


# ============================================================
# PASO 2 — DIFERENCIACIÓN: FST (Weir & Cockerham) y Dest (Jost)
# ============================================================

# --- 2a. FST pareado + IC 95% por bootstrap sobre loci ---
fst_pareado <- as.matrix(pairwise.WCfst(hf_10))
set.seed(123)
fst_boot <- boot.ppfst(hf_10, nboot = 1000, quant = c(0.025, 0.975))

# --- 2b. [v2] Dest global y por locus, con IC 95% (bootstrap) ---
dest_obs <- D_Jost(genind_10)
cat("\n=== Dest (Jost 2008) ===\n"); print(dest_obs)
set.seed(123)
dest_bs <- chao_bootstrap(genind_10, nreps = 1000)
dest_ic <- summarise_bootstrap(dest_bs, D_Jost)
cat("\n=== Dest con IC 95% (bootstrap) ===\n"); print(dest_ic)

dest_pareado <- as.matrix(pairwise_D(genind_10, linearized = FALSE))

# --- 2c. [v2] Permutaciones pareadas para FST y Dest (999) + FDR ---
cat("\nPermutaciones pareadas de FST y Dest (999 réplicas; puede tardar)...\n")
set.seed(123)
nperm <- 999
pares <- combn(pops_names, 2, simplify = FALSE)

resultados_pareados <- do.call(rbind, lapply(pares, function(par) {
  sub <- genind_10[pop(genind_10) %in% par]
  pop(sub) <- droplevels(pop(sub))
  fst_o  <- pairwise.WCfst(genind2hierfstat(sub))[1, 2]
  dest_o <- as.numeric(pairwise_D(sub, linearized = FALSE))
  perm <- replicate(nperm, {
    s <- sub; pop(s) <- sample(pop(sub))      # reparte individuos entre las dos poblaciones
    c(fst  = tryCatch(pairwise.WCfst(genind2hierfstat(s))[1, 2], error = function(e) NA),
      dest = tryCatch(as.numeric(pairwise_D(s, linearized = FALSE)), error = function(e) NA))
  })
  cat("  ", par[1], "-", par[2], "listo\n")
  data.frame(
    Pop1 = par[1], Pop2 = par[2],
    FST = round(fst_o, 4),
    FST_IC_inf = round(fst_boot$ll[par[1], par[2]], 4),
    FST_IC_sup = round(fst_boot$ul[par[1], par[2]], 4),
    FST_p  = (sum(perm["fst", ]  >= fst_o,  na.rm = TRUE) + 1) / (nperm + 1),
    Dest = round(dest_o, 4),
    Dest_p = (sum(perm["dest", ] >= dest_o, na.rm = TRUE) + 1) / (nperm + 1)
  )
}))
resultados_pareados$FST_p_FDR  <- round(p.adjust(resultados_pareados$FST_p,  "BH"), 4)
resultados_pareados$Dest_p_FDR <- round(p.adjust(resultados_pareados$Dest_p, "BH"), 4)
cat("\n=== FST y Dest pareados ===\n"); print(resultados_pareados)


# ============================================================
# PASO 3 — [v2] DAPC (no interactivo)
# ============================================================
# Criterios:
#  - find.clusters: se retienen prácticamente todos los PCs (n.pca = 200),
#    como recomienda Jombart para k-means (no hay riesgo de sobreajuste en
#    este paso), y K se elige con BIC (criterio "diffNgroup").
#  - DAPC a priori: el número de PCs se elige con optim.a.score, que busca
#    el número que maximiza la reasignación correcta sin sobreajustar.

N_PCA_MAX <- min(200, nInd(genind_10) - 1)

set.seed(123)
grp <- find.clusters(genind_10, max.n.clust = 10, n.pca = N_PCA_MAX,
                     choose.n.clust = FALSE, criterion = "diffNgroup")
K_kmeans <- length(levels(grp$grp))
bic_tabla <- data.frame(K = 1:length(grp$Kstat), BIC = round(as.numeric(grp$Kstat), 2))
cat("\n=== find.clusters: K elegido por BIC =", K_kmeans, "===\n"); print(bic_tabla)

n_da <- length(pops_names) - 1
dapc_full <- dapc(genind_10, pop(genind_10), n.pca = N_PCA_MAX, n.da = n_da)
set.seed(123)
oa <- optim.a.score(dapc_full, plot = FALSE)
n_pca_opt <- oa$best
dapc_apriori <- dapc(genind_10, pop(genind_10), n.pca = n_pca_opt, n.da = n_da)

var_explicada <- round(100 * dapc_apriori$eig / sum(dapc_apriori$eig), 1)
reasignacion  <- mean(as.character(dapc_apriori$assign) == as.character(pop(genind_10)))
cat("\n=== DAPC a priori ===\n")
cat("PCs retenidos (optim.a.score):", n_pca_opt, "\n")
cat("Varianza discriminante LD1 =", var_explicada[1], "% ; LD2 =", var_explicada[2], "%\n")
cat("Reasignación correcta a la población de origen:", round(100 * reasignacion, 1), "%\n")

tabla_clusters <- data.frame(ID = indNames(genind_10), Poblacion = as.character(pop(genind_10)),
                             Cluster_kmeans = grp$grp)
cat("\nPoblación x grupo de k-means:\n"); print(table(tabla_clusters$Poblacion, tabla_clusters$Cluster_kmeans))

# --- Figura DAPC ---
coords_ind <- as.data.frame(dapc_apriori$ind.coord[, 1:2])
colnames(coords_ind) <- c("LD1", "LD2")
coords_ind$Codigo <- as.character(pop(genind_10))
centroides <- coords_ind %>% group_by(Codigo) %>% summarise(LD1 = mean(LD1), LD2 = mean(LD2))
nombres_completos <- c(AN = "Anaime", CO = "Cocora occidental", CX = "Cocora oriental",
                       LI = "La Línea", PI = "Pijao", SF = "San Félix",
                       TE = "Tenerife", TO = "Toche")

fig_dapc <- ggplot(coords_ind, aes(x = LD1, y = LD2, color = Codigo)) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  geom_point(size = 2, alpha = 0.7) +
  stat_ellipse(aes(group = Codigo), linetype = "dotted", linewidth = 0.6) +
  geom_point(data = centroides, aes(x = LD1, y = LD2, fill = Codigo),
             shape = 21, size = 5, color = "black", stroke = 1) +
  geom_text_repel(data = centroides, aes(x = LD1, y = LD2, label = Codigo),
                  color = "black", size = 3.5, box.padding = 0.6, max.overlaps = Inf) +
  scale_color_discrete(labels = nombres_completos) +
  labs(x = paste0("LD1 (", var_explicada[1], "%)"), y = paste0("LD2 (", var_explicada[2], "%)"),
       color = "Población") +
  theme_classic(base_size = 12, base_family = "sans") +   # Biotropica: Arial, sin negrita
  guides(fill = "none")
print(fig_dapc)
ggsave("outputs/Figura_DAPC.tiff", fig_dapc, width = 17, height = 12, units = "cm", dpi = 600)

ggsave("figuras/Figura_DAPC.png", fig_dapc, width = 17, height = 12, units = "cm", dpi = 400)
# ============================================================
# PASO 4 — [v2] snmf
# ============================================================

mat_geno <- tab(genind_10, NA.method = "asis")    # [v2] faltantes como NA, no como promedio
mat_geno[is.na(mat_geno)] <- 9                     # 9 = dato faltante en formato .geno
write.geno(mat_geno, "ceroxylon_10loci.geno")

proyecto_snmf <- snmf("ceroxylon_10loci.geno", K = 1:10, repetitions = 10,
                      entropy = TRUE, project = "new", CPU = 2, seed = 123)   # [v2] semilla

ce_tabla <- data.frame(
  K = 1:10,
  CE_media = sapply(1:10, function(k) mean(cross.entropy(proyecto_snmf, K = k))),
  CE_min   = sapply(1:10, function(k) min(cross.entropy(proyecto_snmf, K = k))),
  CE_de    = sapply(1:10, function(k) sd(cross.entropy(proyecto_snmf, K = k)))
)
cat("\n=== Entropía cruzada por K ===\n"); print(round(ce_tabla, 4))
cat("K con menor entropía cruzada media:", ce_tabla$K[which.min(ce_tabla$CE_media)], "\n")

ancestria_K <- function(proyecto, K_val, genind_obj) {
  run <- which.min(LEA::cross.entropy(proyecto, K = K_val))
  Q_mat <- LEA::Q(proyecto, K = K_val, run = run)
  tabla <- data.frame(ID = indNames(genind_obj), Poblacion = as.character(pop(genind_obj)), Q_mat)
  colnames(tabla)[3:ncol(tabla)] <- paste0("C", 1:K_val)
  tabla
}

armar_panel_K <- function(proyecto, K_val, genind_obj) {
  tabla <- ancestria_K(proyecto, K_val, genind_obj)
  tabla_larga <- tabla %>% pivot_longer(cols = starts_with("C"), names_to = "Cluster", values_to = "Ancestria")
  orden_pops <- c("AN","CO","PI","SF","TE","CX","TO","LI")
  tabla_larga$Poblacion <- factor(tabla_larga$Poblacion, levels = orden_pops)
  tabla_larga <- tabla_larga %>% arrange(Poblacion)
  tabla_larga$ID <- factor(tabla_larga$ID, levels = unique(tabla_larga$ID))
  ggplot(tabla_larga, aes(x = ID, y = Ancestria, fill = Cluster)) +
    geom_bar(stat = "identity", width = 1) +
    facet_grid(~ Poblacion, scales = "free_x", space = "free_x") +
    theme_minimal(base_size = 10, base_family = "sans") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.spacing = unit(0.05, "lines"),
          strip.text = if (K_val == 2) element_text() else element_blank(),
          legend.position = "none", axis.title.x = element_blank()) +
    labs(y = paste0("K = ", K_val)) +
    scale_fill_brewer(palette = "Set1")
}
# Mismo estilo que las demás figuras de la tesis (fuente, tamaños, guardar())
source("scripts/00_estilo_figuras.R")
# Ocho colores distinguibles con daltonismo (Okabe e Ito)
col_pob <- c(AN = "#E69F00", CO = "#009E73", CX = "#56B", LI = "#D55E00",
             PI = "#CC79A7", SF = "#0072B2", TE = "#8C6D31", TO = "#000000")

fig_dapc <- ggplot(coords_ind, aes(x = LD1, y = LD2, colour = Codigo)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  stat_ellipse(aes(group = Codigo), linetype = "dashed", linewidth = 0.4, level = 0.95) +
  geom_point(size = 1.3, alpha = 0.6) +
  geom_point(data = centroides, aes(fill = Codigo), shape = 21, size = 3.2,
             colour = TINTA, stroke = 0.5) +
  geom_text_repel(data = centroides, aes(label = Codigo), colour = TINTA,
                  box.padding = 0.5, min.segment.length = 0, seed = 1) +
  scale_colour_manual(values = col_pob, labels = nombres_completos, name = "Población") +
  scale_fill_manual(values = col_pob, guide = "none") +
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 2.2, alpha = 1))) +
  labs(x = paste0("LD1 (", var_explicada[1], " %)"), y = paste0("LD2 (", var_explicada[2], " %)"))
print(fig_dapc)
guardar(fig_dapc, "Figura_DAPC", ANCHO_COMPLETO, 11)
ggsave("outputs/Figura_DAPC.tiff", fig_dapc, width = ANCHO_COMPLETO, height = 11,
       units = "cm", dpi = 600)


# Ancestría media por población (para Resultados)
Q_medias <- lapply(c(2, 4, 6), function(k) {
  t <- ancestria_K(proyecto_snmf, k, genind_10)
  cbind(K = k, aggregate(. ~ Poblacion, data = t[, -1], FUN = function(x) round(mean(x), 3)))
})
names(Q_medias) <- paste0("K", c(2, 4, 6))
print(Q_medias)


ggsave("figuras/Figura_snmf.png", fig_snmf, width = 17, height = 14, units = "cm", dpi = 400)

# ============================================================
# PASO 5 — AMOVA (10 loci)
# ============================================================
# poppr.amova elimina por defecto los loci con más del 5 % de datos faltantes.
# Ca30 tiene 12.3 % (todos los genotipos de TO), así que se usa cutoff = 0.15
# para conservar los 10 loci, como en el resto de los análisis.

genclone_10 <- as.genclone(genind_10)
grupo_pob <- ifelse(pop(genind_10) %in% c("CX","LI","TO"), "CV", "GR")
strata(genclone_10) <- data.frame(Grupo = factor(grupo_pob), Pop = pop(genind_10))

# Modelo 1 — sin jerarquía
amova_1 <- poppr.amova(genclone_10, ~Pop, within = TRUE, cutoff = 0.15)
set.seed(123); amova_1_test <- randtest(amova_1, nrepet = 999)

# Modelo 2 — jerárquico, poblaciones anidadas en la fuente de los datos
amova_2 <- poppr.amova(genclone_10, ~Grupo/Pop, within = TRUE, cutoff = 0.15)
set.seed(123); amova_2_test <- randtest(amova_2, nrepet = 999)

tabla_amova <- function(a, t, modelo) {
  cc <- a$componentsofcovariance
  out <- data.frame(Modelo = modelo, Componente = rownames(cc),
                    Sigma = round(cc$Sigma, 4), Porcentaje = round(cc[["%"]], 2))
  attr_p <- data.frame(Modelo = modelo, Prueba = t$names, p = t$pvalue)
  list(componentes = out, p = attr_p)
}
am1 <- tabla_amova(amova_1, amova_1_test, "Sin jerarquia")
am2 <- tabla_amova(amova_2, amova_2_test, "Jerarquico (fuente)")
cat("\n=== AMOVA (10 loci) ===\n"); print(am1); print(am2)
print(amova_2$statphi)

# ============================================================
# PASO 7 — EXPORTAR A EXCEL
# ============================================================

wb <- createWorkbook()
addWorksheet(wb, "FST_pareado");       writeData(wb, "FST_pareado", round(fst_pareado, 4), rowNames = TRUE)
addWorksheet(wb, "Dest_pareado");      writeData(wb, "Dest_pareado", round(dest_pareado, 4), rowNames = TRUE)
addWorksheet(wb, "Pareados_IC_p_FDR"); writeData(wb, "Pareados_IC_p_FDR", resultados_pareados)
addWorksheet(wb, "Dest_global_locus"); writeData(wb, "Dest_global_locus",
                                                 data.frame(Salida = capture.output(print(dest_obs))))
addWorksheet(wb, "Dest_IC_bootstrap"); writeData(wb, "Dest_IC_bootstrap",
                                                 data.frame(Salida = capture.output(print(dest_ic))))
addWorksheet(wb, "DAPC_resumen");      writeData(wb, "DAPC_resumen", data.frame(
  Parametro = c("PCs find.clusters", "K (BIC, diffNgroup)", "PCs DAPC a priori (optim.a.score)",
                "LD1 (%)", "LD2 (%)", "Reasignacion correcta (%)"),
  Valor = c(N_PCA_MAX, K_kmeans, n_pca_opt, var_explicada[1], var_explicada[2], round(100 * reasignacion, 1))))
addWorksheet(wb, "BIC_kmeans");        writeData(wb, "BIC_kmeans", bic_tabla)
addWorksheet(wb, "Pob_x_kmeans");      writeData(wb, "Pob_x_kmeans",
                                                 as.data.frame.matrix(table(tabla_clusters$Poblacion, tabla_clusters$Cluster_kmeans)),
                                                 rowNames = TRUE)
addWorksheet(wb, "snmf_entropia");     writeData(wb, "snmf_entropia", round(ce_tabla, 5))
for (k in names(Q_medias)) { addWorksheet(wb, paste0("Q_medio_", k)); writeData(wb, paste0("Q_medio_", k), Q_medias[[k]]) }
addWorksheet(wb, "AMOVA");             writeData(wb, "AMOVA", rbind(am1$componentes, am2$componentes))
addWorksheet(wb, "AMOVA_p");           writeData(wb, "AMOVA_p", rbind(am1$p, am2$p))

saveWorkbook(wb, "outputs/Resultados_Estructura_Genetica_10loci_v2.xlsx", overwrite = TRUE)
cat("\n✅ Excel guardado: outputs/Resultados_Estructura_Genetica_10loci_v2.xlsx\n")

saveRDS(fst_pareado,  "fst_pareado_10loci_final.rds")
saveRDS(dest_pareado, "dest_pareado_10loci_final.rds")


# =============================================================================
# BLOQUE: Mapa de calor de F_ST y D_est
# (antes: 16_mapa_calor_fst.R)
# =============================================================================
# =============================================================================
# 16_mapa_calor_fst.R
# Mapas de calor de F_ST y D_est pareados, ordenados por agrupamiento
# jerárquico UPGMA sobre F_ST, con el dendrograma encima.
#
# Usa el objeto resultados_pareados del bloque anterior.
# =============================================================================

source("scripts/00_estilo_figuras.R")
library(dplyr); library(tidyr); library(patchwork); library(ggdendro)
pp   <- resultados_pareados
pops <- sort(unique(c(pp$Pop1, pp$Pop2)))

a_matriz <- function(columna) {
  M <- matrix(0, length(pops), length(pops), dimnames = list(pops, pops))
  for (i in seq_len(nrow(pp))) {
    M[pp$Pop1[i], pp$Pop2[i]] <- pp[[columna]][i]
    M[pp$Pop2[i], pp$Pop1[i]] <- pp[[columna]][i]
  }
  M
}
M_fst  <- a_matriz("FST")
M_dest <- a_matriz("Dest")

# UPGMA sobre F_ST (los valores negativos, si los hubiera, se toman como 0)
hc    <- hclust(as.dist(pmax(M_fst, 0)), method = "average")
orden <- hc$labels[hc$order]
cat("Orden del dendrograma:", orden, "\n")
cat("Alturas de unión:\n"); print(round(hc$height, 3))

a_largo <- function(M) {
  as.data.frame(as.table(M)) %>%
    setNames(c("p", "q", "valor")) %>%
    mutate(p = factor(p, rev(orden)), q = factor(q, orden)) %>%
    filter(p != q)
}

mapa <- function(M, etiqueta, letra) {
  ggplot(a_largo(M), aes(x = q, y = p, fill = valor)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    scale_fill_gradient(low = "#e6eefa", high = "#0b3d91", name = etiqueta,
                        limits = c(0, max(M))) +
    labs(x = NULL, y = NULL, tag = letra) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "bottom", legend.key.width = unit(0.8, "cm"))
}
dd <- dendro_data(as.dendrogram(hc))
dendro <- ggplot(segment(dd)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend),
               colour = TINTA, linewidth = 0.4) +
  scale_x_continuous(limits = c(0.5, length(pops) + 0.5), expand = c(0, 0)) +
  theme_void()

panel_a <- dendro / mapa(M_fst,  expression(F[ST]),  "(a)") + plot_layout(heights = c(1, 4))
panel_b <- dendro / mapa(M_dest, expression(D[est]), "(b)") + plot_layout(heights = c(1, 4))
fig <- (panel_a | panel_b)

guardar(fig, "mapa_calor_fst_dest", ANCHO_COMPLETO, 10)
