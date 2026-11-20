# =============================================================================
# 01_datos.R
# Corrección de tamaños alélicos entre fuentes y matriz final de 10 loci.
#
# Sección de la tesis: Metodología (fuentes de datos, marcadores y estandarización)
# Lee:      Matriz_datos_8poblaciones.xlsx
# Produce:  datos_10loci_final.rds / .xlsx, datos_10loci_preTO.rds
#
# Correr desde la carpeta Tesis_Ceroxylon, en una sesión limpia de R
# (Session > Restart R), después de los scripts con número menor.
# =============================================================================

source("scripts/00_configuracion.R")


# =============================================================================
# BLOQUE: Estandarización
# (antes: Estandarización_10loci.R)
# =============================================================================
library(readxl)
library(dplyr)
library(writexl)

# ============================================================
# ESTANDARIZACIÓN DE LA MATRIZ MULTI-FUENTE — VERSIÓN FINAL (10 loci)
# ============================================================

df <- as.data.frame(read_excel("Matriz_datos_8poblaciones.xlsx"))

poblaciones_gonzalez <- c("AN","CO","PI","SF","TE")  # fuente: González-Rivillas et al. 2018
poblaciones_chacon   <- c("CX","LI","TO")             # fuente: Chacón-Vargas et al. 2019 (referencia, sin corregir)

loci_10 <- c("Cs03","Cs24","Ca17","Ca07","Ca16","Cs05","Ca30","Ca20","Ca25","Cs30")

# Offsets estimados para los 10 loci (diferencia de mediana de tamaño alélico, González - Chacón)
offset_todos <- c(Ca07=13, Ca16=17, Ca17=11, Ca20=5, Ca25=-7, Ca30=17,
                  Cs03=23, Cs05=22, Cs24=18, Cs30=-16)

corregir_genotipo <- function(genotipo, offset) {
  ale <- strsplit(genotipo, "/")[[1]]
  ale_corr <- ifelse(ale == "0", "0", as.character(as.numeric(ale) - offset))
  paste(ale_corr, collapse = "/")
}

df_10 <- df
for (locus in loci_10) {
  idx <- df_10$Localidad %in% poblaciones_gonzalez
  df_10[idx, locus] <- sapply(df_10[idx, locus], corregir_genotipo, offset = offset_todos[[locus]])
}

# ============================================================
# VERIFICACIÓN 1 — rango de valores (ningún alelo corregido debería ser negativo)
# ============================================================
cat("=== Verificación de rangos por locus ===\n")
for (l in loci_10) {
  vals <- unlist(strsplit(df_10[[l]], "/"))
  vals <- as.numeric(vals[vals != "0"])
  cat(l, ": min =", min(vals), " max =", max(vals), "\n")
}

# ============================================================
# VERIFICACIÓN 2 — variación real dentro de cada población de González
# ============================================================
cat("\n=== Verificación de variación intrapoblacional (debe ser > 1 en todas) ===\n")
alerta_encontrada <- FALSE
for (pop in poblaciones_gonzalez) {
  for (l in loci_10) {
    n_unicos <- length(unique(df_10[[l]][df_10$Localidad == pop]))
    if (n_unicos <= 1) {
      cat("⚠️  ALERTA:", pop, "-", l, "solo tiene", n_unicos, "valor(es) único(s)\n")
      alerta_encontrada <- TRUE
    }
  }
}
if (!alerta_encontrada) cat("Sin alertas — la corrección quedó bien.\n")

# ============================================================
# Matriz final de trabajo (10 loci, versión definitiva)
# ============================================================

df_10_final <- df_10[, c("ID", "Localidad", loci_10)]

# Versión antes de dejar TO en Ca30 como faltante (la usa A4_sensibilidad_TO_Ca30.R)
saveRDS(df_10_final, "datos_10loci_preTO.rds")

# TO en Ca30: todos sus alelos quedan 25 pb por encima de los de las demás
# poblaciones (marcador de 25 pb en los geles de Chacón-Vargas). Se dejan como faltantes.
df_10_final$Ca30[df_10_final$Localidad == "TO"] <- "0/0"
saveRDS(df_10_final, "datos_10loci_final.rds")
write_xlsx(df_10_final, "datos_10loci_final.xlsx")

print(table(df_10_final$Localidad))
cat("\n✅ Estandarización de 10 loci guardada: datos_10loci_final.rds / .xlsx\n")

