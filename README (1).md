# Conectividad genética de *Ceroxylon quindiuense*

Scripts de R de la tesis de Maestría en Biomatemáticas (Universidad del Quindío):
conectividad genética de la palma de cera con microsatélites y redes poblacionales.

## Cómo correr

1. Abrir R o RStudio con la carpeta principal del proyecto como directorio de trabajo.
2. Poner en esa carpeta la matriz de genotipos `Matriz_datos_8poblaciones.xlsx`.
3. Correr los scripts en el orden de la tabla, cada uno en una sesión limpia de R.
   Cada script carga primero `scripts/00_configuracion.R`, que instala los paquetes que falten.
4. Al final, correr `scripts/99_control.R`, que compara los valores principales con los de la tesis.

## Scripts

| Script | Sección de la tesis | Contenido |
|---|---|---|
| `00_configuracion.R` | (todos) | Paquetes y carpetas |
| `00_estilo_figuras.R` | (figuras) | Tema y colores de las figuras |
| `01_datos.R` | Metodología | Corrección entre fuentes y matriz de 10 loci |
| `01b_mapa_area_estudio.R` | Metodología | Mapa del área de estudio |
| `02_diversidad.R` | Resultados: diversidad; Anexo B | Diversidad, alelos nulos y privados, F_IS |
| `03_estructura.R` | Resultados: estructura | F_ST, D_est, DAPC, sNMF y AMOVA |
| `04_red.R` | Resultados: red | Etapas 1 a 3, estabilidad, roles y mapas de la red |
| `05_metricas_markov.R` | Resultados: red | Métricas del grafo, Mantel y cadenas de Markov |
| `06_simulaciones.R` | Resultados: escenarios | Pérdida de poblaciones y corredores |
| `A1_correccion.R` | Anexo A | Evaluación de la corrección |
| `A2_agrupamiento.R` | Anexo A | Agrupamiento de alelos separados por 1 pb |
| `A3_sensibilidad_loci.R` | Anexo A | Sensibilidad sin Ca25 ni Cs30 |
| `A4_sensibilidad_TO_Ca30.R` | Metodología | Efecto de los alelos de TO en Ca30 |
| `99_control.R` | (control) | Revisión de los valores reportados |

## Archivos intermedios

| Archivo | Lo produce | Lo usan |
|---|---|---|
| `datos_10loci_final.rds` | 01 | 02, 04, A1–A3 |
| `datos_10loci_preTO.rds` | 01 | A4 |
| `genind_10loci_final.rds` | 02 | 03, 99 |
| `fst_pareado_10loci_final.rds` | 03 | 04, 99 |
| `red_10loci.rds` | 04 | 05, 99 |
| `edges_df_10loci.rds` | 04 | 05, 06 |

## Semillas

Estructura (permutaciones, DAPC, sNMF y AMOVA): 123. divMigrate, estabilidad de la red,
agrupamiento aleatorio y sensibilidad: 2026.

## Datos

Los genotipos provienen de González-Rivillas et al. (2018) y Chacón-Vargas (2015).
