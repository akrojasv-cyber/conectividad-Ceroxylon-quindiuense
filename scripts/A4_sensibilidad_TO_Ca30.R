# =============================================================================
# A4_sensibilidad_TO_Ca30.R
# Efecto de los alelos de TO en Ca30: original, corregido (-25 pb) y faltante.
#
# Sección de la tesis: Metodología (decisión de dejar TO en Ca30 como faltante)
# Lee:      datos_10loci_preTO.rds
# Produce:  outputs/ con el resumen de las tres versiones
#
# Cambio: lee la matriz anterior a la decisión, para que la versión "original"
# siga siendo la original al volver a correrlo.
#
# Correr desde la carpeta Tesis_Ceroxylon, en una sesión limpia de R
# (Session > Restart R), después de los scripts con número menor.
# =============================================================================

source("scripts/00_configuracion.R")


# =============================================================================
# BLOQUE: Sensibilidad TO en Ca30
# (antes: 15_sensibilidad_TO_Ca30.R)
# =============================================================================
# =============================================================================
# 15_sensibilidad_TO_Ca30.R
# ¿Cambian los resultados por los alelos de TO en Ca30?
#
# En Ca30, las otras siete poblaciones tienen alelos entre 96 y 132 pb, y TO
# solo tiene 150 y 152, exactamente 25 pb por encima de 125 y 127. Varios geles
# de Toche y La Línea se leyeron con un marcador de 25 pb (Chacón-Vargas,
# tesis), así que lo más probable es un error de un paso del marcador.
#
# Versiones:
#   original    : matriz de la tesis
#   corregido   : alelos de TO en Ca30 menos 25 pb
#   faltante    : genotipos de TO en Ca30 como datos faltantes ("0/0")
#
# En cada versión: privados, Ho/He/F_IS de TO, F_ST pareado (pares de TO y
# media general), aristas de la Etapa 2 y divMigrate (estimación puntual) en
# las aristas de TO.
#
# Usa datos_10loci_preTO.rds (lo guarda 01_datos.R).
# =============================================================================

library(dplyr); library(tidyr); library(MASS)
library(adegenet); library(hierfstat); library(diveRsity)
select <- dplyr::select

df   <- readRDS("datos_10loci_preTO.rds")   # matriz antes de dejar TO en Ca30 como faltante
loci <- c("Cs03","Cs24","Ca17","Ca07","Ca16","Cs05","Ca30","Ca20","Ca25","Cs30")
pops <- sort(unique(df$Localidad))

# ---------------------------------------------------------------- versiones --
corregir_TO <- function(d, modo) {
  i <- d$Localidad == "TO" & d$Ca30 != "0/0"
  if (modo == "corregido") {
    d$Ca30[i] <- sapply(strsplit(d$Ca30[i], "/"),
                        function(a) paste(as.numeric(a) - 25, collapse = "/"))
  } else if (modo == "faltante") {
    d$Ca30[i] <- "0/0"
  }
  d
}
versiones <- list(original  = df,
                  corregido = corregir_TO(df, "corregido"),
                  faltante  = corregir_TO(df, "faltante"))

cat("\n== Control: alelos de TO en Ca30 en cada versión ==\n")
for (v in names(versiones)) {
  cat(v, ":", paste(names(table(versiones[[v]]$Ca30[versiones[[v]]$Localidad == "TO"])),
                    collapse = " "), "\n")
}

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
a_genind <- function(d) {
  geno <- d[, loci]
  geno[] <- lapply(geno, function(x) ifelse(x == "0/0", NA, x))
  df2genind(geno, sep = "/", ind.names = d$ID, pop = d$Localidad,
            ploidy = 2, type = "codom", NA.char = NA)
}
aristas_poda <- function(fr) {
  fw <- fr %>% mutate(col = paste(locus, alelo, sep = "_")) %>%
    select(Localidad, col, freq) %>%
    pivot_wider(names_from = col, values_from = freq, values_fill = 0) %>%
    arrange(Localidad)
  X <- as.matrix(fw[, -1]); rownames(X) <- fw$Localidad
  Z <- scale(X, center = TRUE, scale = FALSE); m <- ncol(Z)
  O <- ginv((Z %*% t(Z)) / (m - 1)); dd <- sqrt(diag(O))
  R <- -O / outer(dd, dd); dimnames(R) <- list(rownames(X), rownames(X))
  pr <- combn(rownames(X), 2)
  r <- R[cbind(pr[1, ], pr[2, ])]
  z <- 0.5 * log((1 + r) / (1 - r)) * sqrt(m - 6 - 3)
  P <- 2 * (1 - pnorm(abs(z)))
  data.frame(par = paste(pr[1, ], pr[2, ], sep = "-"), r = round(r, 3),
             arista = P < 0.05 / 28)
}
genepop_geno <- function(genotipo) {
  ale <- strsplit(genotipo, "/")[[1]]
  ale_fmt <- ifelse(ale == "0", "000", formatC(as.numeric(ale), width = 3, format = "d", flag = "0"))
  paste0(ale_fmt[1], ale_fmt[2])
}
migrar <- function(d) {
  f <- tempfile(fileext = ".txt")
  lineas <- c("Genepop - sensibilidad", loci)
  for (p in pops) {
    lineas <- c(lineas, "Pop")
    sub <- d[d$Localidad == p, ]
    for (i in seq_len(nrow(sub))) {
      genos <- sapply(loci, function(l) genepop_geno(sub[[l]][i]))
      lineas <- c(lineas, paste0(p, "_", i, " ,  ", paste(genos, collapse = " ")))
    }
  }
  writeLines(lineas, f)
  g <- divMigrate(infile = f, outfile = NULL, boots = 0, stat = "gst",
                  filter_threshold = 0, plot_network = FALSE)$gRelMig
  unlink(f); dimnames(g) <- list(pops, pops); g
}

# ---------------------------------------------------------------- análisis --
res <- lapply(names(versiones), function(v) {
  d   <- versiones[[v]]
  fr  <- frecuencias(d)
  gen <- a_genind(d)
  hf  <- genind2hierfstat(gen)
  bs  <- basic.stats(hf)
  fst <- as.matrix(pairwise.WCfst(hf))
  priv <- fr %>% group_by(locus, alelo) %>% filter(n_distinct(Localidad) == 1) %>%
    ungroup() %>% count(Localidad, name = "privados")
  ar <- aristas_poda(fr)
  g  <- migrar(d)
  list(
    resumen = data.frame(
      version = v,
      alelos_totales = nrow(distinct(fr, locus, alelo)),
      privados_total = sum(priv$privados),
      privados_TO = sum(priv$privados[priv$Localidad == "TO"]),
      Ho_TO  = round(mean(bs$Ho[, "TO"], na.rm = TRUE), 3),
      He_TO  = round(mean(bs$Hs[, "TO"], na.rm = TRUE), 3),
      Fis_TO = round(mean(bs$Fis[, "TO"], na.rm = TRUE), 3),
      Fst_Ca30 = round(bs$perloc["Ca30", "Fst"], 3),
      Fst_medio_TO = round(mean(fst["TO", colnames(fst) != "TO"]), 3),
      Fst_medio_todos = round(mean(fst[upper.tri(fst)]), 3),
      aristas = sum(ar$arista),
      w_LI_TO = round(g["LI", "TO"], 3),
      w_TO_LI = round(g["TO", "LI"], 3)),
    fst_TO = round(fst["TO", ], 3),
    aristas = ar %>% filter(arista) %>% pull(par)
  )
})
names(res) <- names(versiones)

cat("\n== Resumen por versión ==\n")
print(bind_rows(lapply(res, `[[`, "resumen")), row.names = FALSE)

cat("\n== F_ST pareado de TO con cada población ==\n")
print(do.call(rbind, lapply(res, `[[`, "fst_TO")))

cat("\n== Aristas de la Etapa 2 ==\n")
for (v in names(res)) cat(sprintf("  %-10s %s\n", v, paste(res[[v]]$aristas, collapse = ", ")))
cat("\n¿Iguales a la original?\n")
for (v in names(res)[-1]) cat(sprintf("  %-10s %s\n", v,
                                      identical(sort(res[[v]]$aristas), sort(res$original$aristas))))

dir.create("outputs", showWarnings = FALSE)
write.csv(bind_rows(lapply(res, `[[`, "resumen")),
          "outputs/15_sensibilidad_TO_Ca30.csv", row.names = FALSE)
