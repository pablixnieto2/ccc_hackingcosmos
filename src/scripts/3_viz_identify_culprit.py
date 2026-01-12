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
OUTPUT_IMG = 'output/fracture_identification.png'
NSIDE = 8 
LINE_LENGTH_DEG = 8.0 

# COORDENADAS DE LOS SOSPECHOSOS HABITUALES (Galactic Coordinates l, b)
# Large Magellanic Cloud (LMC)
LMC_L = 280.46
LMC_B = -32.89  # OJO: En coordenadas galácticas es -32, no -69 (eso es celestial)

# Small Magellanic Cloud (SMC)
SMC_L = 302.80
SMC_B = -44.30

# Centaurus A (Galaxia Radioactiva potente)
CEN_A_L = 309.5
CEN_A_B = 19.4

def main():
    print("--- IDENTIFICACIÓN DE OBJETIVOS ---")
    df = pd.read_csv(INPUT_FILE)
    
    plt.figure(figsize=(15, 10), facecolor='black')
    ax = plt.subplot(111, projection='mollweide', facecolor='black')
    plt.grid(True, color='dimgray', alpha=0.3)
    
    # 1. DIBUJAR TUS LÍNEAS (Igual que antes)
    for i, row in df.iterrows():
        idx = int(row['center_idx'])
        angle_deg = row['angle']
        corr = row['corr_IP']
        theta, phi = hp.pix2ang(NSIDE, idx)
        lon_rad = phi - np.pi 
        lat_rad = np.pi/2 - theta
        
        rot_rad = np.radians(angle_deg)
        length_rad = np.radians(LINE_LENGTH_DEG)
        d_lat = (length_rad / 2) * np.cos(rot_rad)
        d_lon = (length_rad / 2) * np.sin(rot_rad) / np.cos(lat_rad)
        
        color = '#ff3333' if corr > 0 else '#33ffff'
        alpha = min(abs(corr) * 2.5, 0.8) 
        ax.plot([lon_rad - d_lon, lon_rad + d_lon], [lat_rad - d_lat, lat_rad + d_lat], 
                color=color, linewidth=2, alpha=alpha)

    # 2. DIBUJAR LOS SOSPECHOSOS (ETIQUETAS VERDES)
    suspects = [
        ("LMC (Magallanes)", LMC_L, LMC_B),
        ("SMC (Pequeña Magallanes)", SMC_L, SMC_B),
        ("Centaurus A", CEN_A_L, CEN_A_B)
    ]
    
    for name, l, b in suspects:
        # Convertir L, B (Grados) a Radianes para Mollweide
        # Healpy phi es 0..2pi. Mollweide espera -pi..pi.
        # L va de 0 a 360. Si L > 180, restamos 360 para que quede negativo (lado derecho)
        
        # Ajuste de coordenadas:
        phi_rad = np.radians(l)
        if phi_rad > np.pi: phi_rad -= 2*np.pi
        
        # Mollweide en Matplotlib usa Lon [-pi, pi] donde 0 es el centro.
        # Pero OJO: En Astronomía Galáctica, el centro (0,0) es el centro.
        # En Healpy/Matplotlib a veces la convención de signo de Longitud invierte Este/Oeste.
        # Vamos a plotear tal cual.
        
        lon_proj = phi_rad - np.pi # Ajuste usual para centrar el mapa en 0
        # Espera, la conversión anterior era: lon_rad = phi - np.pi. 
        # Healpy phi=0 es L=0? Sí.
        
        theta_rad = np.radians(90 - b) # Colatitude
        lat_proj = np.pi/2 - theta_rad
        
        # Dibujar punto
        ax.plot(lon_proj, lat_proj, 'o', color='#00FF00', markersize=15, mfc='none', markeredgewidth=2)
        ax.text(lon_proj, lat_proj - 0.15, name, color='#00FF00', ha='center', fontsize=12, fontweight='bold')

    plt.title(f"THE VERDICT: Anomalies vs Known Objects", color='white', fontsize=16)
    plt.savefig(OUTPUT_IMG, dpi=300, facecolor='black')
    print(f"Mapa de evidencia guardado en {OUTPUT_IMG}")

if __name__ == "__main__":
    main()