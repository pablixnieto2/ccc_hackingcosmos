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

# CONFIGURACIÓN
INPUT_FILE = 'data/processed/line_metrics.csv'
OUTPUT_IMG = 'output/fracture_map_168.png'
NSIDE = 8 
LINE_LENGTH_DEG = 8.0 

def main():
    print(f"--- GENERANDO MAPA TÁCTICO DE LAS 168 LÍNEAS ---")
    
    # 1. Cargar tus 168 líneas
    try:
        df = pd.read_csv(INPUT_FILE)
    except:
        print("No encuentro el archivo csv. Asegúrate de haber corrido el Line Hunter.")
        return

    print(f"Cargando {len(df)} fracturas detectadas...")

    # Configuración del lienzo (Fondo negro estilo espacio)
    plt.figure(figsize=(15, 10), facecolor='black')
    ax = plt.subplot(111, projection='mollweide', facecolor='black')
    
    # Rejilla tenue para referencia
    plt.grid(True, color='dimgray', alpha=0.4, linestyle=':')
    
    # --- DIBUJAR LAS LÍNEAS ---
    for i, row in df.iterrows():
        idx = int(row['center_idx'])
        angle_deg = row['angle']
        corr = row['corr_IP']
        
        # Coordenadas del centro del píxel
        theta, phi = hp.pix2ang(NSIDE, idx)
        
        # Conversión a proyección Mollweide (Radianes, Lon [-pi, pi])
        lon_rad = phi - np.pi 
        lat_rad = np.pi/2 - theta
        
        # Calcular los extremos de la línea según el ángulo detectado
        # Matemáticas de rotación local
        rot_rad = np.radians(angle_deg)
        length_rad = np.radians(LINE_LENGTH_DEG)
        
        # Delta X y Delta Y proyectados
        # El coseno de latitud corrige la distorsión del ancho en los polos
        d_lat = (length_rad / 2) * np.cos(rot_rad)
        d_lon = (length_rad / 2) * np.sin(rot_rad) / np.cos(lat_rad)
        
        x1 = lon_rad - d_lon
        y1 = lat_rad - d_lat
        x2 = lon_rad + d_lon
        y2 = lat_rad + d_lat
        
        # COLOR CODING:
        # ROJO NEÓN = Correlación Positiva (Calor y Polarización suben juntos)
        # CIAN NEÓN = Correlación Negativa (Uno sube, otro baja)
        color = '#ff3333' if corr > 0 else '#33ffff' # Red / Cyan
        
        # Grosor y opacidad según la fuerza de la señal
        alpha = min(abs(corr) * 2.5, 1.0) 
        linewidth = 1.5 + (abs(corr) * 3)
        
        ax.plot([x1, x2], [y1, y2], color=color, linewidth=linewidth, alpha=alpha, solid_capstyle='round')

        # Highlight para el TOP 3 (Los más fuertes)
        # (Asumimos que el CSV no está ordenado, así que marcamos si corr > 0.23 que vimos antes)
        if abs(corr) > 0.23:
            ax.plot(lon_rad, lat_rad, 'o', color='yellow', markersize=6, alpha=0.8)
            # Solo etiquetar el famoso 647
            if idx == 647:
                plt.text(lon_rad, lat_rad + 0.15, "THE SCAR (647)", color='yellow', 
                         ha='center', fontsize=9, fontweight='bold')

    # Cosmética final
    plt.title(f"FRACTURE LINES MAP: {len(df)} DETECTIONS\n(|Lat| > 30° | Corr > 0.15)", 
              color='white', fontsize=15, pad=20)
    
    # Ejes invisibles pero con etiquetas
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    
    print("Guardando imagen...")
    plt.savefig(OUTPUT_IMG, dpi=300, facecolor='black', bbox_inches='tight')
    print(f"¡Mapa listo! Abre: {OUTPUT_IMG}")

if __name__ == "__main__":
    main()