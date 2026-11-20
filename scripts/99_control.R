# =============================================================================
# 99_control.R
# Compara los números clave de la tesis con lo que dan los objetos guardados.
# Correr al final, en una sesión limpia. Cada línea dice OK o REVISAR.
# Si algo sale REVISAR, lo más probable es que un script se haya corrido con
# objetos viejos en memoria o que falte volver a correr un script posterior.
#
# Lee: genind_10loci_final.rds, fst_pareado_10loci_final.rds, red_10loci.rds
# =============================================================================

source("scripts/00_configuracion.R")
suppressPackageStartupMessages({ library(adegenet); library(hierfstat); library(poppr) })

revisar <- function(que, obtenido, tesis, tol = 0) {
  ok <- all(abs(obtenido - tesis) <= tol + 1e-9)
  cat(sprintf("%-8s %-38s obtenido: %-22s tesis: %s\n", if (ok) "OK" else "REVISAR",
              que, paste(round(obtenido, 3), collapse = " "), paste(tesis, collapse = " ")))
  invisible(ok)
}
revisar_txt <- function(que, obtenido, tesis) {
  ok <- identical(as.character(obtenido), as.character(tesis))
  cat(sprintf("%-8s %-38s obtenido: %s\n", if (ok) "OK" else "REVISAR", que,
              paste(obtenido, collapse = " ")))
}

# --- Fechas de los archivos intermedios (deben ir en orden de los scripts) ---
cat("\n== Fechas de los archivos intermedios ==\n")
arch <- c("datos_10loci_final.rds", "genind_10loci_final.rds", "fst_pareado_10loci_final.rds",
          "red_10loci.rds", "edges_df_10loci.rds", "outputs/metricas_red_10loci.xlsx",
          "outputs/modelo_dinamico_10loci.xlsx", "outputs/simulaciones_conservacion_10loci.xlsx")
print(data.frame(archivo = arch, modificado = format(file.mtime(arch), "%Y-%m-%d %H:%M")))

# --- 3.1 Diversidad ---
cat("\n== Diversidad ==\n")
g  <- readRDS("genind_10loci_final.rds")
revisar("Número de loci", nLoc(g), 10)
revisar("Alelos en total", sum(nAll(g)), 293)
conteo <- apply(tab(genind2genpop(g, quiet = TRUE)), 2, function(x) sum(x > 0))
revisar("Alelos privados", sum(conteo == 1), 115)
hf <- genind2hierfstat(g)
wcg <- wc(hf)
revisar("F_ST global (Weir y Cockerham)", wcg$FST, 0.177, tol = 0.0005)
revisar("F_IS global (Weir y Cockerham)", wcg$FIS, 0.484, tol = 0.0005)
revisar("He media (Nei)", mean(basic.stats(hf)$Hs, na.rm = TRUE), 0.772, tol = 0.0005)

# --- 3.2 Estructura ---
cat("\n== F_ST pareado de TO (1 - F_ST, como en la red) ==\n")
fst <- readRDS("fst_pareado_10loci_final.rds")
otras <- c("AN", "CO", "CX", "LI", "PI", "SF", "TE")
revisar("1 - F_ST con TO", 1 - fst[otras, "TO"],
        c(0.854, 0.858, 0.752, 0.775, 0.855, 0.835, 0.807), tol = 0.0005)

cat("\n== AMOVA (cutoff = 0.15, porcentajes) ==\n")
gc <- as.genclone(g)
strata(gc) <- data.frame(Grupo = factor(ifelse(pop(g) %in% c("CX", "LI", "TO"), "CV", "GR")),
                         Pop = pop(g))
a1 <- poppr.amova(gc, ~Pop, within = TRUE, cutoff = 0.15, quiet = TRUE)
a2 <- poppr.amova(gc, ~Grupo/Pop, within = TRUE, cutoff = 0.15, quiet = TRUE)
revisar("Sin jerarquía", a1$componentsofcovariance[["%"]][1:3], c(17.22, 39.64, 43.14), tol = 0.005)
revisar("Jerárquico (fuente)", a2$componentsofcovariance[["%"]][1:4],
        c(4.65, 14.29, 38.81, 42.24), tol = 0.005)

# --- 3.3 Red ---
cat("\n== Red ==\n")
red <- readRDS("red_10loci.rds")
revisar("Aristas de la Etapa 2", nrow(red$aristas_finales), 7)
revisar_txt("Roles (AN CO CX LI PI SF TE TO)", red$rol[c("AN","CO","CX","LI","PI","SF","TE","TO")],
            c("Equilibrada","Equilibrada","Aislada","Emisora","Equilibrada","Equilibrada","Equilibrada","Receptora"))
cat("\nSi los roles salen REVISAR, falta correr el bloque 'Roles de nodo' de 04_red.R",
    "\n(o 05 y 06 se corrieron antes de ese bloque).\n")
