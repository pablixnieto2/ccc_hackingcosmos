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
import matplotlib.pyplot as plt
import os

# --- COORDENADAS DEL VECINO 1 (EL GANADOR) ---
TARGET_LAT = -70.8927
TARGET_LON = 136.2065
INPUT_FILE = 'data/raw/COM_CMB_IQU-sevem_2048_R4.00.fits'

# Configuración
SONAR_RADIUS = 15.0   # Radio de búsqueda de pared
STEP_SIZE = 0.2
CORR_THRESHOLD_FACTOR = 0.4

def get_signal_strength(lat, lon, map_I, map_P, nside):
    # Evitar errores de rango
    if lat > 90: lat = 180 - lat
    if lat < -90: lat = -180 - lat
    
    pix = hp.ang2pix(nside, np.radians(90 - lat), np.radians(lon))
    return abs(map_I[pix] * map_P[pix])

def move_geodesic(lat, lon, angle, dist):
    lat_r, lon_r = np.radians(lat), np.radians(lon)
    ang_r, dist_r = np.radians(angle), np.radians(dist)
    new_lat_r = np.arcsin(np.sin(lat_r)*np.cos(dist_r) + np.cos(lat_r)*np.sin(dist_r)*np.cos(ang_r))
    new_lon_r = lon_r + np.arctan2(np.sin(ang_r)*np.sin(dist_r)*np.cos(lat_r), np.cos(dist_r)-np.sin(lat_r)*np.sin(new_lat_r))
    return np.degrees(new_lat_r), np.degrees(new_lon_r)

def main():
    print("🛸 SABUESO V10: NEIGHBOR TRACER - OBJETIVO SUR")
    print(f"📍 Desplegando en Vecino 1: Lat {TARGET_LAT:.4f}, Lon {TARGET_LON:.4f}")
    
    maps = hp.read_map(INPUT_FILE, field=[0,1,2], verbose=False)
    map_I, map_P = maps[0], np.sqrt(maps[1]**2 + maps[2]**2)
    nside = hp.get_nside(map_I)

    # 1. FASE DE SONAR (Buscar la pared)
    print("📡 Escaneando perímetro para encontrar muro de carga...")
    best_wall_lat, best_wall_lon = None, None
    max_sig = 0
    wall_bearing = 0
    
    # Escaneo radial
    for ang in range(0, 360, 10):
        for r in [8, 10, 12, 15]: # Distancias típicas al borde
            l, lo = move_geodesic(TARGET_LAT, TARGET_LON, ang, r)
            sig = get_signal_strength(l, lo, map_I, map_P, nside)
            if sig > max_sig:
                max_sig = sig
                best_wall_lat, best_wall_lon = l, lo
                wall_bearing = ang

    if best_wall_lat is None:
        print("❌ No se detectó muro claro. La señal es difusa.")
        return

    print(f"✅ ¡MURO DETECTADO! Lat {best_wall_lat:.2f}, Lon {best_wall_lon:.2f} (Señal: {max_sig:.2e})")
    print("🕷️ Iniciando rastreo perimetral...")

    # 2. FASE DE RASTREO (Modo Araña)
    current_lat, current_lon = best_wall_lat, best_wall_lon
    current_bearing = (wall_bearing + 90) % 360 
    
    path = [{'lat': current_lat, 'lon': current_lon, 'type': 'START_WALL'}]
    vertices = 0
    steps_since_vertex = 0
    threshold = max_sig * CORR_THRESHOLD_FACTOR
    
    for s in range(1500): 
        # Mirar adelante
        next_lat, next_lon = move_geodesic(current_lat, current_lon, current_bearing, 0.5)
        sig_ahead = get_signal_strength(next_lat, next_lon, map_I, map_P, nside)
        
        # Lógica de Vértice
        if sig_ahead < threshold and steps_since_vertex > 25:
            print(f"   ⚠️ Vértice detectado en Paso {s}. Reorientando...")
            
            # Escaneo de giro
            best_new_ang = current_bearing
            local_max = 0
            reverse = (current_bearing + 180) % 360
            
            for scan_ang in range(0, 360, 10):
                if abs(scan_ang - reverse) < 45: continue
                l, lo = move_geodesic(current_lat, current_lon, scan_ang, 0.5)
                s_val = get_signal_strength(l, lo, map_I, map_P, nside)
                if s_val > local_max:
                    local_max = s_val
                    best_new_ang = scan_ang
            
            current_bearing = best_new_ang
            path.append({'lat': current_lat, 'lon': current_lon, 'type': 'VERTEX'})
            vertices += 1
            steps_since_vertex = 0
            print(f"   ↪️ Nuevo Rumbo: {current_bearing:.1f}º")
            
            if vertices >= 5:
                print("🏆 ¡PENTÁGONO VECINO CERRADO!")
                break
        else:
            # Autocorrección
            best_adj = 0
            m_val = -1
            for adj in [-15, -5, 0, 5, 15]:
                l, lo = move_geodesic(current_lat, current_lon, current_bearing + adj, STEP_SIZE)
                v = get_signal_strength(l, lo, map_I, map_P, nside)
                if v > m_val:
                    m_val = v
                    best_adj = adj
            current_bearing += best_adj
            steps_since_vertex += 1

        # Mover
        current_lat, current_lon = move_geodesic(current_lat, current_lon, current_bearing, STEP_SIZE)
        path.append({'lat': current_lat, 'lon': current_lon, 'type': 'PATH'})

    # 3. RESULTADOS
    df = pd.DataFrame(path)
    df.to_csv('data/processed/neighbor1_track.csv', index=False)
    
    plt.figure(figsize=(10, 8))
    plt.plot(df['lon'], df['lat'], 'g-', linewidth=2, label='Rastro Vecino 1')
    v_df = df[df['type'] == 'VERTEX']
    plt.scatter(v_df['lon'], v_df['lat'], c='orange', s=100, zorder=5, label='Vértices')
    plt.scatter([TARGET_LON], [TARGET_LAT], c='blue', marker='x', s=150, label='Centro Calculado (Moran)')

    plt.title(f"Mapeo del Vecino 1\nLat {TARGET_LAT}, Lon {TARGET_LON}")
    plt.xlabel("Longitud")
    plt.ylabel("Latitud")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig('data/processed/neighbor1_result.png')
    print("🖼️ Mapa guardado: data/processed/neighbor1_result.png")

if __name__ == "__main__":
    main()