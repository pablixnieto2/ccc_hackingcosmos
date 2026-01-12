# ==============================================================================
#  The Geometry of the Echo: PMN-01 Model Source Code
#  ----------------------------------------------------------------------------
#  (c) 2025 Pablo Miguel Nieto Muñoz
#  License: MIT (See LICENSE file for details)
#  
#  Scientific Citation:
#  Nieto Muñoz, P. M. (2025). "The Geometry of the Echo: Observational 
#  Confirmation of the Chiral Dodecahedral Universe". 
#  Zenodo.
# ==============================================================================

library(ggplot2)
library(dplyr)
library(knitr)

# Configuración
INPUT_FILE <- "data/processed/fractal_metrics.csv"
OUTPUT_DIR <- "output/plots_survivors/"

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# 1. Cargar Datos
df <- read.csv(INPUT_FILE)

# 2. Calcular Coordenadas y Aplicar MÁSCARA GALÁCTICA
# Definimos "Zona Segura" como todo lo que esté a más de 20 grados del ecuador.
GALACTIC_CUT <- 20 

df <- df %>%
  mutate(
    theta_deg = theta * 180 / pi,
    phi_deg = phi * 180 / pi,
    lat = 90 - theta_deg,
    lon = phi_deg - 180
  )

# 3. Filtrar: Solo queremos lo que NO sea Galaxia y sea Anómalo
survivors <- df %>%
  filter(
    abs(lat) > GALACTIC_CUT,      # Fuera de la Vía Láctea
    hurst_I > 0.80,               # Fractal fuerte
    corr_IP > 0.20                # Correlación positiva (Temperatura y Polarización unidas)
  ) %>%
  arrange(desc(corr_IP)) # Ordenar por los más fuertes

# 4. Resultados en Texto
print(paste("--- REPORTE DE SUPERVIVIENTES ---"))
print(paste("Total de anillos analizados:", nrow(df)))
print(paste("Anillos eliminados por la Máscara Galáctica (+/-", GALACTIC_CUT, "deg):", nrow(df %>% filter(abs(lat) <= GALACTIC_CUT))))
print(paste("Candidatos restantes en Cielo Profundo:", nrow(survivors)))

if (nrow(survivors) > 0) {
  print("TOP 5 CANDIDATOS:")
  print(survivors %>% select(id_anillo, lat, lon, hurst_I, corr_IP) %>% head(5))
  
  # Guardar CSV con los supervivientes
  write.csv(survivors, file.path(OUTPUT_DIR, "final_candidates.csv"), index = FALSE)
} else {
  print("Resultado: 0 Supervivientes. La señal era puramente galáctica.")
}

# 5. Visualización Final (Viz E)
# Pintamos la máscara en gris para ver qué hemos borrado
p_final <- ggplot() +
  # Fondo
  annotate("rect", xmin = -180, xmax = 180, ymin = -90, ymax = 90, fill = "black") +
  
  # Zona de Exclusión Galáctica (Sombreado gris)
  annotate("rect", xmin = -180, xmax = 180, ymin = -GALACTIC_CUT, ymax = GALACTIC_CUT, 
           fill = "grey20", alpha = 0.5) +
  
  # Todos los puntos (muy tenues)
  geom_point(data = df, aes(x = lon, y = lat), color = "grey10", size = 0.1) +
  
  # LOS SUPERVIVIENTES (Si hay alguno, brillará en CYAN neón)
  geom_point(data = survivors, aes(x = lon, y = lat, size = corr_IP), 
             color = "cyan", shape = 16, alpha = 1) +
  
  # Círculo rojo alrededor para encontrarlos rápido
  geom_point(data = survivors, aes(x = lon, y = lat), 
             color = "red", shape = 1, size = 6, stroke = 1) +

  geom_hline(yintercept = c(-GALACTIC_CUT, GALACTIC_CUT), linetype = "dashed", color = "white") +
  
  coord_fixed(ratio = 1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "The Final Survivors",
    subtitle = paste("Deep Sky Anomalies (|Lat| >", GALACTIC_CUT, "°)"),
    x = "Longitude",
    y = "Latitude",
    caption = "Grey Zone = Galactic Mask (Excluded)"
  ) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "black"),
    plot.background = element_rect(fill = "white"),
    text = element_text(color = "black")
  )

ggsave(file.path(OUTPUT_DIR, "viz_E_survivors.png"), plot = p_final, width = 12, height = 7, dpi = 300)
print(paste("Gráfica guardada en", OUTPUT_DIR))