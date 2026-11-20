# =============================================================================
# A3_sensibilidad_loci.R
# Sensibilidad a la exclusión de Ca25 y Cs30 (10 frente a 8 loci).
#
# Sección de la tesis: Anexo A (Sensibilidad sin Ca25 ni Cs30)
# Lee:      datos_10loci_final.rds
# Produce:  outputs/ con el resumen de 10 frente a 8 loci
#
# Correr desde la carpeta Tesis_Ceroxylon, en una sesión limpia de R
# (Session > Restart R), después de los scripts con número menor.
# =============================================================================

source("scripts/00_configuracion.R")


# =============================================================================
# BLOQUE: Sensibilidad 10 frente a 8 loci
# (antes: sensibilidad_10L_vs_8L.R)
# =============================================================================
# ============================================================
# Análisis de sensibilidad: 10 loci (versión final) vs 8 loci (sin Ca25 ni Cs30)
# Artículo Biotropica, Ceroxylon quindiuense
# ============================================================
# PASO 0 (OBLIGATORIO): antes de mirar la columna de 8 loci, verifica que la
# columna de 10 loci reproduzca tu manuscrito:
#   FST global = 0.171 | Tabla 4 | K DAPC = 6 | K snmf = 4 |
#   7 aristas: AN-PI 0.510, CO-SF 0.372, LI-TO 0.313, CO-TE 0.309,
#              AN-CO 0.300, PI-SF 0.294, CO-PI 0.207
# Si un paso NO reproduce tu resultado, reemplaza esa función por tu código original.
# ============================================================

library(adegenet)
library(hierfstat)
library(ade4)
library(LEA)
library(MASS)      # ginv(): pseudo-inversa de Moore-Penrose
library(writexl)

set.seed(2026)

# ---- Parámetros: deben ser LOS MISMOS que usaste en el análisis original ----
NPCA_FC   <- 200            # PCs retenidos en find.clusters  [ajusta al tuyo]
CRIT_FC   <- "diffNgroup"   # criterio de selección de K con BIC [ajusta al tuyo]
K_SNMF    <- 1:10
REP_SNMF  <- 10
ALFA_BONF <- 0.05 / choose(8, 2)   # 0.00179
N_MANTEL  <- 9999

# ---- 0. Datos ----------------------------------------------------------------
df <- readRDS("datos_10loci_final.rds")   # columnas: ID, Localidad, 10 loci "a/b"
loci_10 <- c("Cs03","Cs24","Ca17","Ca07","Ca16","Cs05","Ca30","Ca20","Ca25","Cs30")
loci_8  <- setdiff(loci_10, c("Ca25", "Cs30"))
pops    <- c("AN","CO","CX","LI","PI","SF","TE","TO")

# Convierte la matriz "a/b" en objeto genind; genotipos con alelo "0" = dato faltante
hacer_genind <- function(df, loci) {
  X <- df[, loci]
  X[] <- lapply(X, function(col) {
    col <- as.character(col)
    col[grepl("(^|/)0($|/)", col)] <- NA
    col
  })
  df2genind(X, sep = "/", ind.names = df$ID,
            pop = factor(df$Localidad, levels = pops),
            ploidy = 2, type = "codom")
}

# ---- 1. FST global y pareado (Weir & Cockerham 1984) --------------------------
fst_global <- function(g) wc(genind2hierfstat(g))$FST

fst_pareado <- function(g) {
  m <- as.matrix(pairwise.WCfst(genind2hierfstat(g)))
  m[pops, pops]
}

# ---- 2. K del DAPC (find.clusters + BIC, sin poblaciones a priori) ------------
k_dapc <- function(g) {
  fc <- find.clusters(g, max.n.clust = 10, n.pca = NPCA_FC,
                      choose.n.clust = FALSE, criterion = CRIT_FC)
  length(levels(fc$grp))
}

# ---- 3. K de snmf (mínima entropía cruzada) -----------------------------------
# Cada alelo es una columna con conteos 0/1/2; 9 = dato faltante (formato .geno)
k_snmf <- function(g, etiqueta) {
  G <- tab(g, NA.method = "asis")
  G[is.na(G)] <- 9
  archivo <- paste0("snmf_", etiqueta, ".geno")
  write.geno(G, archivo)
  proj <- snmf(archivo, K = K_SNMF, repetitions = REP_SNMF, entropy = TRUE,
               ploidy = 2, project = "new", seed = 2026)
  ce <- sapply(K_SNMF, function(k) mean(cross.entropy(proj, K = k)))
  list(K = K_SNMF[which.min(ce)], ce = ce)
}

# ---- 4. Poda por independencia condicional (Dyer & Nason 2004) ----------------
podar_red <- function(g) {
  gp <- genind2genpop(g, quiet = TRUE)
  Fq <- tab(gp, freq = TRUE)[pops, ]                 # 8 poblaciones x alelos
  Fc <- scale(Fq, center = TRUE, scale = FALSE)      # centra cada alelo entre poblaciones
  C  <- tcrossprod(Fc) / (ncol(Fq) - 1)              # covarianza 8 x 8 entre poblaciones
  P  <- ginv(C)                                      # matriz de precisión (pseudo-inversa)
  R  <- -P / sqrt(outer(diag(P), diag(P)))           # correlaciones parciales
  dimnames(R) <- list(pops, pops)
  n <- ncol(Fq)                  # n usado en la prueba z [verifica que coincida con el tuyo]
  k <- length(pops) - 2          # variables de control
  pares <- t(combn(pops, 2))
  r <- R[pares]
  z <- atanh(r) * sqrt(n - k - 3)
  p <- 2 * pnorm(-abs(z))
  data.frame(par = paste(pares[, 1], pares[, 2], sep = "-"),
             r_parcial = round(r, 3), p = signif(p, 3),
             retenida = p < ALFA_BONF)
}

# ---- 5. Correr ambas versiones -------------------------------------------------
g10 <- hacer_genind(df, loci_10)
g8  <- hacer_genind(df, loci_8)

fst10 <- fst_pareado(g10)
fst8  <- fst_pareado(g8)
mt    <- mantel.rtest(as.dist(fst10), as.dist(fst8), nrepet = N_MANTEL)

red10 <- podar_red(g10)
red8  <- podar_red(g8)

snmf10 <- k_snmf(g10, "10L")
snmf8  <- k_snmf(g8,  "8L")

resumen <- data.frame(
  Medida = c("FST global (Weir & Cockerham)",
             "K DAPC (find.clusters, BIC)",
             "K snmf (mínima entropía cruzada)",
             "Aristas retenidas (de 28)"),
  Loci_10 = c(round(fst_global(g10), 3), k_dapc(g10), snmf10$K, sum(red10$retenida)),
  Loci_8  = c(round(fst_global(g8), 3),  k_dapc(g8),  snmf8$K,  sum(red8$retenida))
)

# ---- 6. Salidas ----------------------------------------------------------------
cat("\n===== RESUMEN 10 vs 8 loci =====\n"); print(resumen)

cat("\nMantel entre matrices de FST pareado (10L vs 8L): r =",
    round(mt$obs, 3), "; p =", mt$pvalue, "\n")
cat("Diferencia absoluta máxima de FST pareado:",
    round(max(abs(fst10 - fst8), na.rm = TRUE), 4), "\n")

comp_red <- merge(red10, red8, by = "par", suffixes = c("_10L", "_8L"))
cat("\n===== Aristas retenidas en al menos una versión =====\n")
print(comp_red[comp_red$retenida_10L | comp_red$retenida_8L, ])

cat("\nEntropía cruzada media por K (10L):\n"); print(round(snmf10$ce, 4))
cat("Entropía cruzada media por K (8L):\n");   print(round(snmf8$ce, 4))

a_tabla <- function(m) cbind(Poblacion = rownames(m), as.data.frame(round(m, 4)))
write_xlsx(list(Resumen    = resumen,
                Mantel     = data.frame(r = mt$obs, p = mt$pvalue),
                FST_10L    = a_tabla(fst10),
                FST_8L     = a_tabla(fst8),
                Red_10L_8L = comp_red,
                snmf_CE    = data.frame(K = K_SNMF, CE_10L = snmf10$ce, CE_8L = snmf8$ce)),
           "sensibilidad_10L_vs_8L.xlsx")
cat("\nGuardado: sensibilidad_10L_vs_8L.xlsx\n")
