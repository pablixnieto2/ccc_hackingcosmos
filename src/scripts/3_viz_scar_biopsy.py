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

import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

# CONFIGURACIÓN
INPUT_FILE = 'data/raw/COM_CMB_IQU-sevem_2048_R4.00.fits'
OUTPUT_IMG = 'output/scar_biopsy_647.png'

# OBJETIVO: EL NÚCLEO DE LA CICATRIZ
# Usamos el índice 647 (NSIDE 8) que detectó el Line Hunter
TARGET_PIX_LOWRES = 647
NSIDE_LOWRES = 8

def main():
    print(f"--- INICIANDO BIOPSIA DE LA CICATRIZ (Pixel {TARGET_PIX_LOWRES}) ---")
    
    # 1. Calcular coordenadas exactas del centro de la cicatriz
    theta, phi = hp.pix2ang(NSIDE_LOWRES, TARGET_PIX_LOWRES)
    lat_deg = 90 - np.degrees(theta)
    lon_deg = np.degrees(phi)
    
    # Ajuste para rotación (Healpy usa Lon 0..360, Gnomview prefiere -180..180 a veces)
    print(f"Coordenadas Objetivo: Lat {lat_deg:.2f}°, Lon {lon_deg:.2f}°")

    # 2. Cargar el mapa FULL RESOLUTION (2048)
    print("Cargando mapa de alta resolución (puede tardar un poco)...")
    try:
        maps = hp.read_map(INPUT_FILE, field=None, hdu=1, verbose=False, memmap=True)
        if len(maps) == 1 or maps.ndim == 1:
            map_I = maps
            maps_pol = hp.read_map(INPUT_FILE, field=None, hdu=2, verbose=False, memmap=True)
            map_Q, map_U = maps_pol[0], maps_pol[1]
        else:
            map_I, map_Q, map_U = maps[0], maps[1], maps[2]
            
        # Calcular Magnitud de Polarización
        map_P = np.sqrt(map_Q**2 + map_U**2)
        
    except Exception as e:
        print(f"Error cargando mapas: {e}")
        return

    # 3. VISUALIZACIÓN DUAL (Intensidad vs Polarización)
    plt.figure(figsize=(16, 8), facecolor='black')
    
    # Definir rotación para centrar la cámara en la cicatriz
    # Rot = [Lon, Lat, Psi] -> Centramos en (lon_deg, lat_deg)
    rot_coords = [lon_deg, lat_deg, 0]
    
    # --- PANEL 1: TEMPERATURA (INTENSIDAD) ---
    # Buscamos estructuras físicas (nubes, líneas)
    hp.gnomview(map_I, rot=rot_coords, xsize=400, ysize=400, reso=3.0, 
                sub=(1, 2, 1), title=f'Intensity (Heat) @ {lon_deg:.1f}, {lat_deg:.1f}',
                unit='K_cmb', cmap='inferno', hold=True, notext=False)
    
    # --- PANEL 2: POLARIZACIÓN (VIBRACIÓN) ---
    # Buscamos "grano" y dirección. Si es CCC, debería tener una textura diferente.
    hp.gnomview(map_P, rot=rot_coords, xsize=400, ysize=400, reso=3.0, 
                sub=(1, 2, 2), title='Polarization Magnitude (Energy)',
                unit='K_cmb', cmap='viridis', norm='hist', hold=True, notext=False)
    
    # Supertítulo
    plt.suptitle(f"THE SCAR BIOPSY: High-Res Analysis of Anomaly 647", color='white', fontsize=18, y=0.95)
    
    print(f"Guardando imagen de alta resolución en {OUTPUT_IMG}...")
    plt.savefig(OUTPUT_IMG, dpi=200, facecolor='black')
    print("¡Biopsia completada!")

if __name__ == "__main__":
    main()