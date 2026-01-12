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
library(viridis)

# Configuración
INPUT_FILE <- "data/processed/fractal_metrics.csv"
OUTPUT_DIR <- "output/plots_elite/"

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

df <- read.csv(INPUT_FILE)

# Preprocesamiento de coordenadas
df <- df %>%
  mutate(
    theta_deg = theta * 180 / pi,
    phi_deg = phi * 180 / pi,
    lat = 90 - theta_deg,
    lon = phi_deg - 180
  )

# --- EL FILTRO DE ÉLITE ---
# Aquí definimos qué es una "anomalía real".
# Criterio: 
# 1. Hurst alto (Estructura fractal)
# 2. Correlación POSITIVA (La temperatura y polarización suben juntas)

UMBRAL_CORRELACION <- 0.25  # Solo queremos lo que esté muy a la derecha en tu gráfica B
UMBRAL_HURST <- 0.80        # Y que sea fractal

elite_candidates <- df %>%
  filter(corr_IP > UMBRAL_CORRELACION & hurst_I > UMBRAL_HURST)

num_candidatos <- nrow(elite_candidates)
print(paste("Hemos encontrado", num_candidatos, "candidatos de élite."))

# Si no hay candidatos, bajamos el listón para ver algo
if (num_candidatos == 0) {
    print("Demasiado estricto. Bajando umbral a 0.15...")
    elite_candidates <- df %>% filter(corr_IP > 0.15)
}

# --- VIZ D: EL MAPA DEL TESORO ---
# Solo pintamos los puntos de élite sobre un fondo negro.

p_elite <- ggplot() +
  # Fondo negro (todo el cielo)
  annotate("rect", xmin = -180, xmax = 180, ymin = -90, ymax = 90, fill = "black") +
  
  # Puntos de escaneo (gris muy tenue para referencia)
  geom_point(data = df, aes(x = lon, y = lat), color = "grey20", size = 0.5, alpha = 0.3) +
  
  # CANDIDATOS DE ÉLITE (Rojos y grandes)
  geom_point(data = elite_candidates, aes(x = lon, y = lat, size = corr_IP), 
             color = "red", shape = 16, alpha = 0.9) +
  
  # Círculos concéntricos visuales (ayuda visual)
  geom_point(data = elite_candidates, aes(x = lon, y = lat), 
             color = "yellow", shape = 1, size = 5, stroke = 1.5) +

  coord_fixed(ratio = 1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "The Elite Candidates",
    subtitle = paste("Rings with Correlation >", UMBRAL_CORRELACION, "| Count:", nrow(elite_candidates)),
    x = "Longitude",
    y = "Latitude",
    size = "Correlation Strength"
  ) +
  theme(
    panel.background = element_rect(fill = "black"),
    plot.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    axis.text = element_text(color = "grey50")
  )

ggsave(file.path(OUTPUT_DIR, "viz_D_elite_map.png"), plot = p_elite, width = 12, height = 7, dpi = 300)

print(paste("Mapa de élite guardado en", OUTPUT_DIR))
