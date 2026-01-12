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

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import os

# --- CONFIGURACIÓN ---
INPUT_FILE = 'data/raw/COM_CMB_IQU-sevem_2048_R4.00.fits'
# Coordenadas del Vértice 647 (Refinado por el Sabueso/Hydra)
START_LAT, START_LON = -41.81, 354.38
# Dirección de la Rama Sólida (Detectada por Hydra)
BEARING = 204.5 
STEP_SIZE = 0.1
MAX_DISTANCE = 45 # Sabemos que anda por 35, miramos hasta 45 por si acaso

def get_structural_signal(lat, lon, map_I, map_P, nside):
    """
    Detecta si estamos en un NODO.
    Un nodo se caracteriza por tener alta complejidad (intersección).
    Usamos la varianza local como proxy de 'cruce de caminos'.
    """
    pix = hp.ang2pix(nside, np.pi/2 - np.radians(lat), np.radians(lon))
    # Vecinos del pixel
    neighbors = hp.get_all_neighbours(nside, pix)
    vals = map_I[neighbors]
    # Si hay mucha varianza, es que hay 'lío' (cruce de paredes de dominio)
    return np.std(vals) * map_P[pix]

def move(lat, lon, angle, step_deg):
    lat_r, lon_r = np.radians(lat), np.radians(lon)
    ang_r, dist_r = np.radians(angle), np.radians(step_deg)
    new_lat_r = np.arcsin(np.sin(lat_r)*np.cos(dist_r) + np.cos(lat_r)*np.sin(dist_r)*np.cos(ang_r))
    new_lon_r = lon_r + np.arctan2(np.sin(ang_r)*np.sin(dist_r)*np.cos(lat_r), np.cos(dist_r)-np.sin(lat_r)*np.sin(new_lat_r))
    return np.degrees(new_lat_r), np.degrees(new_lon_r)

def cosmic_ruler():
    print("📏 LA REGLA CÓSMICA: Midiendo la celda del universo...")
    
    if not os.path.exists(INPUT_FILE):
        print(f"❌ Falta el archivo {INPUT_FILE}")
        return

    maps = hp.read_map(INPUT_FILE, field=[0,1,2], verbose=False)
    map_I, map_P = maps[0], np.sqrt(maps[1]**2 + maps[2]**2)
    nside = hp.get_nside(map_I)

    print(f"   📍 Saliendo de Vértice 647 con Rumbo {BEARING}º...")
    
    curr_lat, curr_lon = START_LAT, START_LON
    dist_travelled = 0
    signal_log = []
    
    # Caminamos paso a paso
    steps = int(MAX_DISTANCE / STEP_SIZE)
    
    for _ in range(steps):
        # 1. Medir "Nodalidad" (¿Parece esto un vértice?)
        sig = get_structural_signal(curr_lat, curr_lon, map_I, map_P, nside)
        signal_log.append(sig)
        
        # 2. Avanzar
        curr_lat, curr_lon = move(curr_lat, curr_lon, BEARING, STEP_SIZE)
        dist_travelled += STEP_SIZE
    
    # --- ANÁLISIS DE LA REGLA ---
    # Buscamos el pico de señal DESPUÉS de haber salido del origen (digamos > 10 grados)
    # Suavizar señal
    signal_smooth = np.convolve(signal_log, np.ones(5)/5, mode='valid')
    x_axis = np.linspace(0, dist_travelled, len(signal_smooth))
    
    # Ignorar los primeros 10 grados (es el vértice de salida)
    mask = x_axis > 10
    search_zone_sig = signal_smooth[mask]
    search_zone_x = x_axis[mask]
    
    if len(search_zone_x) == 0:
        print("❌ Algo falló. No hay datos suficientes.")
        return

    # Encontrar el pico máximo en la zona de búsqueda
    peak_idx = np.argmax(search_zone_sig)
    vertex_distance = search_zone_x[peak_idx]
    
    print("\n=== MEDICIÓN FINAL ===")
    print(f"   📐 Distancia al Siguiente Vértice: {vertex_distance:.4f} GRADOS")
    
    # Gráfico
    plt.style.use('dark_background')
    plt.figure(figsize=(10, 5))
    plt.plot(x_axis, signal_smooth, color='cyan', label='Señal de Estructura (Nodalidad)')
    plt.axvline(vertex_distance, color='yellow', linestyle='--', linewidth=2, label=f'Vértice Destino ({vertex_distance:.1f}º)')
    plt.axvspan(34, 36, color='white', alpha=0.1, label='Zona Teórica (35º)')
    
    plt.title(f"LA REGLA CÓSMICA: ARISTA MEDIDA = {vertex_distance:.2f}º", fontsize=14)
    plt.xlabel("Distancia desde Vértice 647 (Grados)")
    plt.ylabel("Intensidad de Cruce")
    plt.legend()
    plt.savefig('cosmic_ruler_result.png')
    print("💾 Gráfica guardada: cosmic_ruler_result.png")
    plt.show()

if __name__ == "__main__":
    cosmic_ruler()