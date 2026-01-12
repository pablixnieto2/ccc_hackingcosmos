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
import pandas as pd
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import os

# --- CONFIGURACIÓN TÁCTICA ---
INPUT_FILE = 'data/raw/COM_CMB_IQU-sevem_2048_R4.00.fits'
BRANCH_3_FILE = 'data/processed/branch_3.csv' # Cargamos la rama sur
OUTPUT_EXTENDED = 'data/processed/branch_3_extended.csv'
STEP_SIZE = 0.1
SCAN_RADIUS = 30.0  # Cuánto vamos a profundizar (en grados)
MAX_STEPS = int(SCAN_RADIUS / STEP_SIZE)

def get_local_corr(lat, lon, angle, map_I, map_P, nside):
    dist = 0.5
    lat_r, lon_r = np.radians(lat), np.radians(lon)
    ang_r, dist_r = np.radians(angle), np.radians(dist)
    nl_r = np.arcsin(np.sin(lat_r)*np.cos(dist_r) + np.cos(lat_r)*np.sin(dist_r)*np.cos(ang_r))
    nlo_r = lon_r + np.arctan2(np.sin(ang_r)*np.sin(dist_r)*np.cos(lat_r), np.cos(dist_r)-np.sin(lat_r)*np.sin(nl_r))
    pix = hp.ang2pix(nside, np.pi/2 - nl_r, nlo_r)
    return abs(map_I[pix] * map_P[pix])

def main():
    print("🐶 --- SABUESO HYDRA V6: DEEP SCAN (RAMA 3) ---")
    
    # 1. Cargar Datos
    print("🛰️ Cargando mapas y rastro previo...")
    maps = hp.read_map(INPUT_FILE, field=[0,1,2], verbose=False)
    map_I, map_P = maps[0], np.sqrt(maps[1]**2 + maps[2]**2)
    nside = hp.get_nside(map_I)
    
    df_prev = pd.read_csv(BRANCH_3_FILE)
    curr_lat = df_prev.iloc[-1]['lat']
    curr_lon = df_prev.iloc[-1]['lon']
    # Estimar rumbo actual basado en los últimos puntos
    prev_lat = df_prev.iloc[-5]['lat']
    prev_lon = df_prev.iloc[-5]['lon']
    curr_ang = np.degrees(np.arctan2(curr_lon - prev_lon, curr_lat - prev_lat))

    extended_path = df_prev.to_dict('records')
    vertex_found = False

    # 2. RASTREO DE LARGA DISTANCIA
    print(f"🚀 Iniciando exploración desde Lat {curr_lat:.2f}, Lon {curr_lon:.2f}")
    for s in range(MAX_STEPS):
        # Escaneo de dirección para corregir rumbo
        best_a, max_v = curr_ang, 0
        for a in np.arange(curr_ang-20, curr_ang+20, 1):
            v = get_local_corr(curr_lat, curr_lon, a, map_I, map_P, nside)
            if v > max_v:
                max_v, best_a = v, a
        
        curr_ang = best_a
        
        # Moverse 0.1º
        lat_r, lon_r = np.radians(curr_lat), np.radians(curr_lon)
        d_r, a_r = np.radians(STEP_SIZE), np.radians(curr_ang)
        curr_lat = np.degrees(np.arcsin(np.sin(lat_r)*np.cos(d_r) + np.cos(lat_r)*np.sin(d_r)*np.cos(a_r)))
        curr_lon = np.degrees(lon_r + np.arctan2(np.sin(a_r)*np.sin(d_r)*np.cos(lat_r), np.cos(d_r)-np.sin(lat_r)*np.sin(np.radians(curr_lat))))
        
        extended_path.append({'lat': curr_lat, 'lon': curr_lon, 'corr': max_v})

        # DETECTOR DE VÉRTICES: Si hay un pico lateral fuerte, sospechamos vértice
        side_v = get_local_corr(curr_lat, curr_lon, curr_ang + 120, map_I, map_P, nside)
        if side_v > max_v * 0.8 and s > 50:
            print(f"\n⚠️ ¡VÉRTICE POTENCIAL DETECTADO! en Lat {curr_lat:.2f}, Lon {curr_lon:.2f}")
            vertex_found = True
            # No paramos, seguimos para mapear la salida

        if s % 50 == 0:
            print(f"👣 Paso {s}/{MAX_STEPS} | Lat: {curr_lat:.2f} | Rumbo: {curr_ang:.1f}º")

    # 3. Guardar y Dibujar
    df_final = pd.DataFrame(extended_path)
    df_final.to_csv(OUTPUT_EXTENDED, index=False)
    print(f"\n💾 Datos extendidos guardados en: {OUTPUT_EXTENDED}")

    plt.figure(figsize=(10, 8))
    plt.plot(df_final['lon'], df_final['lat'], 'g-', label='Rama 3 Extendida')
    plt.scatter([df_final.iloc[0]['lon']], [df_final.iloc[0]['lat']], c='red', label='Vértice 647')
    if vertex_found:
        plt.scatter([curr_lon], [curr_lat], c='yellow', marker='*', s=200, label='Candidato Vértice 648')
    plt.title("Mapeo de Larga Distancia - Rama 3")
    plt.xlabel("Longitud")
    plt.ylabel("Latitud")
    plt.legend()
    plt.grid(True)
    plt.savefig('data/processed/long_range_scan.png')
    print("🖼️ Mapa generado: data/processed/long_range_scan.png")

if __name__ == "__main__":
    main()