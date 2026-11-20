# =============================================================================
# A1_correccion.R
# Evaluación de la corrección entre fuentes: cuantiles, equivalencias, pares y cuatro escenarios.
#
# Sección de la tesis: Anexo A (Evaluación de la corrección)
# Lee:      datos_10loci_final.rds
# Produce:  outputs/12_*.csv, figuras/12_distribuciones_alelicas
#
# Correr desde la carpeta Tesis_Ceroxylon, en una sesión limpia de R
# (Session > Restart R), después de los scripts con número menor.
# =============================================================================

source("scripts/00_configuracion.R")


# =============================================================================
# BLOQUE: Evaluación de la corrección
# (antes: 12_evaluacion_correccion_delta.R)
# =============================================================================
# =============================================================================
# 12_evaluacion_correccion_delta.R
# Responde a las seis sugerencias de VH sobre la corrección a'_ik = a_ik - delta_k
#
#   1. Distribuciones alélicas completas por locus y fuente, antes y después.
#   2. Si el desplazamiento es constante a lo largo de la distribución
#      (cuantiles 10, 25, 50, 75 y 90).
#   3. Alelos que coinciden entre fuentes después de corregir, los ambiguos
#      (a 1 pb) y los que no tienen par.
#   4. Qué cambia y qué no con la corrección (dentro de fuente vs entre fuentes).
#   5. Escenarios alternativos de corrección y su efecto en F_ST, privados y red.
#   6. Acotar el efecto biológico: delta por par de poblaciones, CO frente a CX
#      y (si hay datos) etapa de vida.
#
# Idea central: delta_k = tau_k (técnico) + beta_k (biológico).
#   - Un corrimiento fijo solo cambia el nombre de los alelos de GR, no sus
#     frecuencias, así que lo que está dentro de una fuente no depende de delta.
#   - Lo que compara fuentes sí depende de delta, y eso es lo que se evalúa.
#   - Si las medianas de poblaciones de una misma fuente difieren poco entre sí,
#     la parte biológica entre fuentes probablemente también es pequeña.
#
# Entradas: datos_10loci_final.rds (matriz corregida de la tesis).
# Los tamaños crudos de GR se recuperan como tam + delta_k (verificado con el
# script 10: rds - xlsx = -delta_k en todos los individuos).
# =============================================================================

library(dplyr); library(tidyr); library(ggplot2)
library(adegenet); library(hierfstat); library(MASS)
select <- dplyr::select

dir.create("outputs", showWarnings = FALSE)
dir.create("figuras", showWarnings = FALSE)

df   <- readRDS("datos_10loci_final.rds")
loci <- c("Cs03","Cs24","Ca17","Ca07","Ca16","Cs05","Ca30","Ca20","Ca25","Cs30")
fuente <- c(AN = "GR", CO = "GR", PI = "GR", SF = "GR", TE = "GR",
            CX = "CV", LI = "CV", TO = "CV")
delta <- c(Cs03 = 23, Cs24 = 18, Ca17 = 11, Ca07 = 13, Ca16 = 17,
           Cs05 = 22, Ca30 = 17, Ca20 = 5, Ca25 = -7, Cs30 = -16)

# Copias alélicas en formato largo: tam = corregido (tesis), crudo = original
copias <- df %>%
  select(ID, Localidad, all_of(loci)) %>%
  pivot_longer(all_of(loci), names_to = "locus", values_to = "g") %>%
  filter(!is.na(g), g != "0/0") %>%
  separate(g, into = c("a1", "a2"), sep = "/", convert = TRUE) %>%
  pivot_longer(c(a1, a2), values_to = "tam") %>% select(-name) %>%
  mutate(fuente = unname(fuente[Localidad]),
         crudo  = ifelse(fuente == "GR", tam + delta[locus], tam))

cat("\n== Control: mediana(GR crudo) - mediana(CV) debe ser delta_k ==\n")
print(copias %>% group_by(locus) %>%
        summarise(dif_medianas = median(crudo[fuente == "GR"]) - median(crudo[fuente == "CV"]),
                  delta = delta[first(locus)]))

# -----------------------------------------------------------------------------
# 1. Distribuciones completas, antes y después
# -----------------------------------------------------------------------------
fig_dat <- bind_rows(
  copias %>% mutate(momento = "Antes de la corrección", x = crudo),
  copias %>% mutate(momento = "Después de la corrección", x = tam)) %>%
  mutate(momento = factor(momento, c("Antes de la corrección", "Después de la corrección")),
         fuente  = factor(fuente, c("GR", "CV")))

p1 <- ggplot(fig_dat, aes(x, y = after_stat(density), fill = fuente)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.55) +
  facet_wrap(vars(locus, momento), ncol = 4, scales = "free") +
  scale_fill_manual(values = c(GR = "#2c3e9e", CV = "#c0392b"),
                    labels = c(GR = "González-Rivillas et al. (2018)",
                               CV = "Chacón-Vargas et al. (2019)"), name = NULL) +
  labs(x = "Tamaño alélico (pb)", y = "Proporción de copias") +
  theme_bw(base_size = 8) +
  theme(legend.position = "top", strip.text = element_text(size = 6.5))
ggsave("figuras/12_distribuciones_alelicas.png", p1, width = 19, height = 26, units = "cm", dpi = 400)
ggsave("figuras/12_distribuciones_alelicas.pdf", p1, width = 19, height = 26, units = "cm")

# -----------------------------------------------------------------------------
# 2. ¿El desplazamiento es constante? Diferencia GR - CV en varios cuantiles
#    (datos crudos). Si es un corrimiento puro, todas las columnas q valen
#    lo mismo y el rango es cercano a cero.
# -----------------------------------------------------------------------------
qs <- c(0.10, 0.25, 0.50, 0.75, 0.90)
cuant <- copias %>%
  group_by(locus, fuente) %>%
  reframe(q = qs, valor = quantile(crudo, qs, type = 1)) %>%
  pivot_wider(names_from = fuente, values_from = valor) %>%
  mutate(dif = GR - CV)
desplazamiento <- cuant %>%
  select(locus, q, dif) %>%
  pivot_wider(names_from = q, names_prefix = "q", values_from = dif) %>%
  left_join(cuant %>% group_by(locus) %>%
              summarise(rango = max(dif) - min(dif)), by = "locus") %>%
  mutate(delta = delta[locus])
cat("\n== Diferencia GR - CV por cuantil (pb, datos crudos) ==\n")
print(desplazamiento, n = Inf)
write.csv(desplazamiento, "outputs/12_desplazamiento_cuantiles.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. Coincidencias y ambigüedades después de corregir
#    Para cada alelo de una fuente, distancia al alelo más cercano de la otra.
#    0 pb = coincide; 1 pb = ambiguo; >= 2 pb = sin par.
# -----------------------------------------------------------------------------
clasificar <- function(desde, hacia) {
  ref <- copias %>% filter(fuente == hacia) %>% distinct(locus, tam)
  copias %>% filter(fuente == desde) %>% count(locus, tam, name = "copias") %>%
    rowwise() %>%
    mutate(dist = min(abs(ref$tam[ref$locus == locus] - tam))) %>%
    ungroup() %>%
    mutate(clase = cut(dist, c(-Inf, 0, 1, Inf),
                       labels = c("Coincide", "Ambiguo (1 pb)", "Sin par (>= 2 pb)")),
           fuente = desde)
}
equiv <- bind_rows(clasificar("GR", "CV"), clasificar("CV", "GR"))
tabla_equiv <- equiv %>%
  group_by(fuente, locus, clase) %>%
  summarise(alelos = n(), copias = sum(copias), .groups = "drop") %>%
  group_by(fuente, locus) %>%
  mutate(prop_copias = round(copias / sum(copias), 3)) %>% ungroup()
cat("\n== Equivalencias entre fuentes después de corregir (resumen) ==\n")
print(equiv %>% group_by(fuente, clase) %>%
        summarise(alelos = n(), copias = sum(copias), .groups = "drop") %>%
        group_by(fuente) %>% mutate(prop_copias = round(copias / sum(copias), 3)))
write.csv(tabla_equiv, "outputs/12_equivalencias_por_locus.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 4 y 5. Escenarios de corrección
#    sin_correccion   : no se toca GR (s = 0)
#    m13_18pb         : 18 pb en los ocho loci con desplazamiento positivo,
#                       0 en Ca25 y Cs30 (solo la causa técnica conocida)
#    mediana          : delta_k (el usado en la tesis)
#    max_solapamiento : el corrimiento entero que maximiza el solapamiento de
#                       las distribuciones de frecuencia de las dos fuentes
#    En todos, alelo GR = crudo - s_k.
# -----------------------------------------------------------------------------
solapamiento <- function(l) {
  g <- copias %>% filter(locus == l, fuente == "GR") %>% count(crudo) %>% mutate(p = n / sum(n))
  v <- copias %>% filter(locus == l, fuente == "CV") %>% count(tam)   %>% mutate(p = n / sum(n))
  ss <- (delta[l] - 30):(delta[l] + 30)
  ov <- sapply(ss, function(s) {
    m <- match(g$crudo - s, v$tam)
    sum(pmin(g$p, ifelse(is.na(m), 0, v$p[m])))
  })
  mejores <- ss[ov == max(ov)]
  mejores[which.min(abs(mejores - delta[l]))]
}
escenarios <- list(
  sin_correccion   = setNames(rep(0, 10), loci),
  m13_18pb         = setNames(ifelse(delta > 0, 18, 0), loci),
  mediana          = delta,
  max_solapamiento = setNames(sapply(loci, solapamiento), loci)
)
cat("\n== Corrimiento s_k de cada escenario ==\n")
print(do.call(rbind, escenarios))

aplicar <- function(s) {
  d <- df
  for (l in loci) {
    i <- fuente[d$Localidad] == "GR" & !is.na(d[[l]]) & d[[l]] != "0/0"
    d[[l]][i] <- sapply(strsplit(d[[l]][i], "/"),
                        function(a) paste(as.numeric(a) + delta[l] - s[l], collapse = "/"))
  }
  d
}

# Funciones del script 08 (idénticas)
a_genind <- function(d) {
  geno <- d[, loci]
  geno[] <- lapply(geno, function(x) ifelse(x == "0/0", NA, x))
  df2genind(geno, sep = "/", ind.names = d$ID, pop = d$Localidad,
            ploidy = 2, type = "codom", NA.char = NA)
}
frecuencias <- function(d) {
  d %>% select(Localidad, all_of(loci)) %>%
    pivot_longer(all_of(loci), names_to = "locus", values_to = "g") %>%
    filter(g != "0/0") %>%
    separate(g, into = c("a1", "a2"), sep = "/") %>%
    pivot_longer(c(a1, a2), values_to = "alelo") %>%
    count(Localidad, locus, alelo) %>%
    group_by(Localidad, locus) %>% mutate(freq = n / sum(n)) %>% ungroup()
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
  pares <- combn(rownames(X), 2)
  r <- R[cbind(pares[1, ], pares[2, ])]
  z <- 0.5 * log((1 + r) / (1 - r)) * sqrt(m - 6 - 3)
  P <- 2 * (1 - pnorm(abs(z)))
  data.frame(par = paste(pares[1, ], pares[2, ], sep = "-"), r = r,
             arista = P < 0.05 / 28, m = m)
}

evaluar <- function(d, nombre) {
  fr  <- frecuencias(d)
  gen <- a_genind(d)
  hf  <- genind2hierfstat(gen)
  bs  <- basic.stats(hf)
  fst <- as.matrix(pairwise.WCfst(hf))
  pr  <- rownames(fst)
  tipo <- outer(fuente[pr], fuente[pr],
                function(a, b) ifelse(a == b, paste0("dentro_", a), "entre"))
  ut <- upper.tri(fst)
  fst_tipo <- tapply(fst[ut], tipo[ut], mean)
  compartidos <- fr %>% mutate(f = fuente[Localidad]) %>%
    distinct(locus, alelo, f) %>% count(locus, alelo) %>% filter(n == 2) %>% nrow()
  priv <- fr %>% group_by(locus, alelo) %>% filter(n_distinct(Localidad) == 1) %>%
    ungroup() %>% count(Localidad, name = "privados") %>% mutate(escenario = nombre)
  ar <- aristas_poda(fr) %>% mutate(escenario = nombre)
  resumen <- data.frame(
    escenario = nombre,
    alelos_totales = nrow(distinct(fr, locus, alelo)),
    alelos_compartidos_GR_CV = compartidos,
    privados_totales = sum(priv$privados),
    He_media_GR = round(mean(colMeans(bs$Hs, na.rm = TRUE)[names(fuente)[fuente == "GR"]]), 4),
    Fst_dentro_GR = round(fst_tipo["dentro_GR"], 4),
    Fst_dentro_CV = round(fst_tipo["dentro_CV"], 4),
    Fst_entre_fuentes = round(fst_tipo["entre"], 4),
    Fst_CO_CX = round(fst["CO", "CX"], 4),
    aristas = sum(ar$arista))
  list(resumen = resumen, priv = priv, aristas = ar)
}

res <- lapply(names(escenarios), function(n) evaluar(aplicar(escenarios[[n]]), n))
resumen_esc <- bind_rows(lapply(res, `[[`, "resumen"))
cat("\n== Escenarios de corrección ==\n")
cat("   He y F_ST dentro de cada fuente deben ser iguales en todos (control).\n")
print(resumen_esc, row.names = FALSE)

aristas_esc <- bind_rows(lapply(res, `[[`, "aristas"))
tabla_aristas <- aristas_esc %>% filter(arista) %>%
  select(escenario, par) %>% mutate(presente = "x") %>%
  pivot_wider(names_from = escenario, values_from = presente, values_fill = "") %>%
  arrange(par)
cat("\n== Aristas presentes en cada escenario ==\n")
print(tabla_aristas, n = Inf)

priv_esc <- bind_rows(lapply(res, `[[`, "priv")) %>%
  pivot_wider(names_from = escenario, values_from = privados, values_fill = 0)
cat("\n== Alelos privados por población y escenario ==\n")
print(priv_esc)

write.csv(resumen_esc,   "outputs/12_escenarios_resumen.csv", row.names = FALSE)
write.csv(tabla_aristas, "outputs/12_escenarios_aristas.csv", row.names = FALSE)
write.csv(priv_esc,      "outputs/12_escenarios_privados.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 6a. Delta por par de poblaciones (datos crudos)
#     - Pares entre fuentes (5 x 3 = 15 por locus): GR - CV con signo.
#       Si el desplazamiento fuera solo técnico, los 15 valores serían
#       parecidos entre sí. Su dispersión indica el tamaño de beta_k.
#     - Pares dentro de una misma fuente: |diferencia| de medianas.
#       Es la escala de la variación biológica entre poblaciones.
# -----------------------------------------------------------------------------
med_pob <- copias %>% group_by(locus, Localidad, fuente) %>%
  summarise(med = median(crudo), .groups = "drop")
pares_med <- inner_join(med_pob, med_pob, by = "locus", suffix = c("1", "2"),
                        relationship = "many-to-many") %>%
  filter(Localidad1 < Localidad2) %>%
  mutate(tipo = ifelse(fuente1 == fuente2, "dentro", "entre"),
         dif = case_when(fuente1 == "GR" & fuente2 == "CV" ~ med1 - med2,
                         fuente1 == "CV" & fuente2 == "GR" ~ med2 - med1,
                         TRUE ~ abs(med1 - med2)))
delta_pares <- pares_med %>% group_by(locus) %>%
  summarise(delta = delta[first(locus)],
            entre_mediana = median(dif[tipo == "entre"]),
            entre_min     = min(dif[tipo == "entre"]),
            entre_max     = max(dif[tipo == "entre"]),
            dentro_mediana_abs = median(dif[tipo == "dentro"]),
            dentro_max_abs     = max(dif[tipo == "dentro"]))
cat("\n== Delta por pares de poblaciones (pb) ==\n")
cat("   entre_*: 15 pares GR-CV. dentro_*: pares de la misma fuente (escala biológica).\n")
print(delta_pares, n = Inf)
write.csv(delta_pares, "outputs/12_delta_por_pares.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 6b. CO frente a CX: mismo valle, distinta fuente
# -----------------------------------------------------------------------------
co_cx <- med_pob %>% filter(Localidad %in% c("CO", "CX")) %>%
  select(locus, Localidad, med) %>%
  pivot_wider(names_from = Localidad, values_from = med) %>%
  mutate(dif_CO_CX = CO - CX, delta = delta[locus], residuo = dif_CO_CX - delta)
cat("\n== CO frente a CX (medianas crudas) ==\n")
print(co_cx, n = Inf)
write.csv(co_cx, "outputs/12_CO_vs_CX.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 6c. Etapa de vida (solo si hay una columna 'etapa' por individuo en GR)
#     Si VH envía la etapa de cada individuo, unirla por ID y descomentar.
# -----------------------------------------------------------------------------
# etapas <- read.csv("etapas_GR.csv")            # columnas: ID, etapa
# d_et <- df %>% inner_join(etapas, by = "ID")
# for (p in c("AN","CO","PI","SF","TE")) {
#   sub <- d_et %>% filter(Localidad == p)
#   g <- a_genind(sub); pop(g) <- sub$etapa
#   cat("\n", p, ": F_ST entre etapas\n"); print(pairwise.WCfst(genind2hierfstat(g)))
# }
