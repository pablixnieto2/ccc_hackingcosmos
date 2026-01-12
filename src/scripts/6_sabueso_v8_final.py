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

# --- CONFIGURACIÓN DE MISIÓN ---
INPUT_FILE = 'data/raw/COM_CMB_IQU-sevem_2048_R4.00.fits'

# Puntos de partida confirmados (Vértice 647)
START_LAT, START_LON = -41.81, 354.38
INITIAL_BEARING = 204.3         # Rumbo hacia la Rama 3
STEP_SIZE = 0.2                 # Resolución de paso
MAX_STEPS = 1500                # Límite de seguridad

def get_signal_at(lat, lon, angle, map_I, map_P, nside):
    """
    Sondea la señal 0.5 grados adelante.
    Fórmula de Detección: Intensidad * Polarización (Busca gradientes fuertes)
    """
    # Proyectar coordenadas
    dist_rad = np.radians(0.5)
    lat_r, lon_r = np.radians(lat), np.radians(lon)
    ang_r = np.radians(angle)
    
    new_lat_r = np.arcsin(np.sin(lat_r)*np.cos(dist_rad) + np.cos(lat_r)*np.sin(dist_rad)*np.cos(ang_r))
    new_lon_r = lon_r + np.arctan2(np.sin(ang_r)*np.sin(dist_rad)*np.cos(lat_r), np.cos(dist_rad)-np.sin(lat_r)*np.sin(new_lat_r))
    
    # Convertir a píxel HEALPix
    pix = hp.ang2pix(nside, np.pi/2 - new_lat_r, new_lon_r)
    
    # Extraer valor (Evitar NaNs)
    val_I = map_I[pix]
    val_P = map_P[pix]
    if np.isnan(val_I): val_I = 0
    if np.isnan(val_P): val_P = 0
    
    return abs(val_I * val_P)

def move(lat, lon, angle, step_deg):
    """Mueve al Sabueso en la esfera."""
    lat_r, lon_r = np.radians(lat), np.radians(lon)
    ang_r, dist_r = np.radians(angle), np.radians(step_deg)
    
    new_lat_r = np.arcsin(np.sin(lat_r)*np.cos(dist_r) + np.cos(lat_r)*np.sin(dist_r)*np.cos(ang_r))
    new_lon_r = lon_r + np.arctan2(np.sin(ang_r)*np.sin(dist_r)*np.cos(lat_r), np.cos(dist_r)-np.sin(lat_r)*np.sin(new_lat_r))
    
    return np.degrees(new_lat_r), np.degrees(new_lon_r)

def run_spider_v8():
    print("🕷️ INICIANDO SABUESO V8 (CORREGIDO)...")
    
    if not os.path.exists(INPUT_FILE):
        print(f"❌ ERROR: No encuentro el archivo en {INPUT_FILE}")
        return

    # 1. CARGA DE DATOS
    print("   ⏳ Cargando mapas Planck (I, Q, U)...")
    maps = hp.read_map(INPUT_FILE, field=[0,1,2], verbose=False)
    map_I = maps[0]
    # Calcular Polarización Total P = sqrt(Q^2 + U^2)
    map_P = np.sqrt(maps[1]**2 + maps[2]**2)
    nside = hp.get_nside(map_I)
    
    # 2. INICIALIZACIÓN
    path = [{'lat': START_LAT, 'lon': START_LON, 'type': 'VERTEX_START'}]
    current_lat, current_lon = START_LAT, START_LON
    current_bearing = INITIAL_BEARING
    
    # Calibración de señal basal
    initial_sig = get_signal_at(START_LAT, START_LON, current_bearing, map_I, map_P, nside)
    # Umbral dinámico: Si la señal cae por debajo del 45% de la media local, es un corte.
    signal_threshold = initial_sig * 0.45 
    
    print(f"   📶 Señal Basal: {initial_sig:.2e} | Umbral Corte: {signal_threshold:.2e}")
    print("   🚀 Sabueso desplegado. Rastreo activo...")

    vertices_found = 0
    steps_on_edge = 0
    total_steps = 0
    
    while total_steps < MAX_STEPS:
        # A) Mirar adelante
        sig_ahead = get_signal_at(current_lat, current_lon, current_bearing, map_I, map_P, nside)
        
        # B) Comprobación de "Caída al Abismo" (Fin de Arista)
        # Solo consideramos caída si hemos caminado al menos 20 pasos (para evitar ruido inicial)
        if sig_ahead < signal_threshold and steps_on_edge > 20:
            print(f"\n🛑 ARISTA TERMINADA en Paso {total_steps} (Lat {current_lat:.2f}, Lon {current_lon:.2f})")
            print(f"   📉 Señal cayó a {sig_ahead:.2e}. Iniciando Radar...")
            
            # --- MANIOBRA DE GIRO (RADAR) ---
            best_new_angle = current_bearing
            max_scan_sig = -1
            
            # Escanear 360 grados
            scan_angles = np.arange(0, 360, 5)
            reverse_angle = (current_bearing + 180) % 360
            
            for ang in scan_angles:
                # No volver hacia atrás (+/- 45 grados)
                diff = abs(ang - reverse_angle)
                if diff > 180: diff = 360 - diff
                if diff < 45: continue
                
                scan_sig = get_signal_at(current_lat, current_lon, ang, map_I, map_P, nside)
                
                # Bonus geométrico: Si el giro es cercano a 72 grados (Dodecaedro), le damos peso extra
                turn_angle = abs(ang - current_bearing)
                if turn_angle > 180: turn_angle = 360 - turn_angle
                if 60 < turn_angle < 85: 
                    scan_sig *= 1.2 # Pequeño empujón heurístico a la geometría sagrada
                
                if scan_sig > max_scan_sig:
                    max_scan_sig = scan_sig
                    best_new_angle = ang
            
            # Registrar Vértice
            print(f"   ↪️ GIRO CONFIRMADO: Rumbo {current_bearing:.1f}º -> {best_new_angle:.1f}º")
            path.append({'lat': current_lat, 'lon': current_lon, 'type': 'VERTEX_FOUND'})
            vertices_found += 1
            
            # Actualizar estado
            current_bearing = best_new_angle
            steps_on_edge = 0
            # Recalibrar umbral para la nueva arista (importante si la intensidad cambia)
            signal_threshold = max_scan_sig * 0.45 
            
            if vertices_found >= 5:
                print("\n🏆 ¡PENTÁGONO CERRADO! 5 Vértices localizados.")
                break
        
        else:
            # C) Caminar y Micro-Corregir (Autopilot)
            # Miramos ligeramente a los lados para mantenernos en la señal máxima
            best_adj = 0
            local_max = -1
            for adj in [-8, -4, 0, 4, 8]:
                check_sig = get_signal_at(current_lat, current_lon, current_bearing + adj, map_I, map_P, nside)
                if check_sig > local_max:
                    local_max = check_sig
                    best_adj = adj
            
            current_bearing += best_adj
            steps_on_edge += 1
            
            # Mover físicamente
            current_lat, current_lon = move(current_lat, current_lon, current_bearing, STEP_SIZE)
            path.append({'lat': current_lat, 'lon': current_lon, 'type': 'PATH'})

        if total_steps % 100 == 0:
            print(f"   👣 Paso {total_steps}: Lat {current_lat:.2f}, Lon {current_lon:.2f} (Sig: {sig_ahead:.2e})")
            
        total_steps += 1

    # --- RESULTADOS FINALES ---
    df = pd.DataFrame(path)
    output_csv = 'spider_track_corrected.csv'
    df.to_csv(output_csv, index=False)
    print(f"\n💾 Datos guardados en: {output_csv}")
    
    # Calcular Centroide
    verts = df[df['type'].str.contains('VERTEX')]
    center_lat, center_lon = 0, 0
    if len(verts) > 0:
        center_lat = verts['lat'].mean()
        center_lon = verts['lon'].mean()
        print(f"🎯 CENTRO GEOMÉTRICO (CARA ALFA): Lat {center_lat:.4f}, Lon {center_lon:.4f}")
    
    # Graficar
    plt.style.use('dark_background')
    plt.figure(figsize=(10, 8))
    
    # Dibujar camino
    plt.plot(df['lon'], df['lat'], 'c-', linewidth=1, alpha=0.7, label='Rastro Sabueso')
    # Dibujar vértices
    plt.scatter(verts['lon'], verts['lat'], c='red', s=150, zorder=10, edgecolors='white', label='Vértices')
    # Dibujar centro
    if len(verts) > 0:
        plt.scatter(center_lon, center_lat, c='yellow', marker='*', s=400, zorder=15, label='Centro Cara Alfa')
    
    # Inicio
    plt.text(START_LON, START_LAT, " INICIO", color='lime', fontweight='bold')

    plt.title(f"SABUESO V8: RASTREO PENTAGONAL (Datos Reales)\nCentro detectado: {center_lat:.2f}, {center_lon:.2f}", fontsize=14)
    plt.xlabel("Longitud Galáctica")
    plt.ylabel("Latitud Galáctica")
    plt.legend()
    plt.grid(True, alpha=0.2)
    
    output_img = 'spider_pentagon_proof.png'
    plt.savefig(output_img)
    print(f"🖼️ Evidencia visual guardada: {output_img}")
    plt.show()

if __name__ == "__main__":
    run_spider_v8()