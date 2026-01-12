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

import pandas as pd
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

# CONFIGURACIÓN
INPUT_FILE = 'data/processed/vertex_trace_647.csv'
OUTPUT_IMG = 'output/vertex_microscope_647.png'

# CENTRO EXACTO (Para poner el (0,0) del gráfico)
CENTER_LAT = -41.81
CENTER_LON = 354.38
NSIDE_TRACE = 256  # El que usamos para generar los datos

def adjust_longitude(lon, center_lon):
    """Corrige el salto de 360 a 0 grados para que el gráfico no se rompa."""
    delta = lon - center_lon
    if delta > 180:  lon -= 360
    if delta < -180: lon += 360
    return lon

def main():
    print(f"--- GENERANDO RADIOGRAFÍA DEL VÉRTICE (4.700+ Segmentos) ---")
    
    try:
        df = pd.read_csv(INPUT_FILE)
    except:
        print("No se encuentra el archivo CSV.")
        return

    print(f"Procesando {len(df)} vectores...")

    # Configuración del Lienzo
    plt.figure(figsize=(12, 12), facecolor='black')
    ax = plt.subplot(111, facecolor='black')
    
    lines = []
    colors = []
    linewidths = []
    
    # Preparamos los datos para LineCollection (mucho más rápido que plotear uno a uno)
    for i, row in df.iterrows():
        idx = int(row['center_idx'])
        angle_deg = row['angle']
        length_deg = row['length']
        corr = row['corr_IP']
        
        # 1. Obtener coordenadas del píxel
        theta, phi = hp.pix2ang(NSIDE_TRACE, idx)
        lat = 90 - np.degrees(theta)
        lon = np.degrees(phi)
        
        # Ajustar Longitud para evitar el corte del mapa (360->0)
        lon = adjust_longitude(lon, CENTER_LON)
        
        # 2. Calcular Deltas relativos al centro (Zoom Local)
        # Convertimos a coordenadas cartesianas locales (grados de diferencia)
        dx_center = (lon - CENTER_LON) * np.cos(np.radians(lat)) # Ajuste de latitud
        dy_center = lat - CENTER_LAT
        
        # 3. Calcular el vector (la línea pequeña)
        # El ángulo 0 es Norte (eje Y positivo). 90 es Este (eje X positivo).
        # Convertimos a radianes matemáticos estándar (0 es derecha)
        # Ángulo matemático = 90 - Ángulo Geográfico
        math_angle_rad = np.radians(90 - angle_deg)
        
        # Longitud visual de la línea (escalada para que no sea enorme en el zoom)
        # Usamos un factor visual, no el grado real, para que se vea limpio
        vis_len = 0.15 * (length_deg / 4.0) 
        
        vx = (vis_len / 2) * np.cos(math_angle_rad)
        vy = (vis_len / 2) * np.sin(math_angle_rad)
        
        # Puntos de inicio y fin de la línea
        p1 = (dx_center - vx, dy_center - vy)
        p2 = (dx_center + vx, dy_center + vy)
        
        lines.append([p1, p2])
        
        # COLOR: Rojo (Positivo) vs Cian (Negativo)
        if corr > 0:
            colors.append((1.0, 0.2, 0.2, min(abs(corr)*3, 1.0))) # Rojo Neón
        else:
            colors.append((0.2, 1.0, 1.0, min(abs(corr)*3, 1.0))) # Cian Neón
            
        # Grosor según fuerza
        linewidths.append(1 + abs(corr)*5)

    # 4. Dibujar Todo de golpe
    lc = LineCollection(lines, colors=colors, linewidths=linewidths, capstyle='round')
    ax.add_collection(lc)
    
    # Ajustar límites del zoom (Radio de 10 grados aprox)
    ax.set_xlim(-8, 8)
    ax.set_ylim(-8, 8)
    
    # Decoración
    ax.grid(True, color='#333333', linestyle='--')
    ax.set_xlabel("Offset Longitud (Grados)", color='gray')
    ax.set_ylabel("Offset Latitud (Grados)", color='gray')
    ax.tick_params(colors='gray')
    
    # Marcar el centro (La Cicatriz Original)
    ax.plot(0, 0, '+', color='yellow', markersize=20, markeredgewidth=2)
    ax.text(0, 0.5, "ORIGIN (647)", color='yellow', ha='center', fontweight='bold')

    plt.title(f"TOPOLOGICAL DEFECT MICROSCOPY\nCenter: Lat {CENTER_LAT}, Lon {CENTER_LON}", color='white', fontsize=16)
    
    print(f"Guardando imagen en {OUTPUT_IMG}...")
    plt.savefig(OUTPUT_IMG, dpi=300, bbox_inches='tight', facecolor='black')
    print("¡Radiografía lista!")

if __name__ == "__main__":
    main()