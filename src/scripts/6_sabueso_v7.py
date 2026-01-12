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

# --- CONFIGURACIÓN DE NAVEGACIÓN ---
INPUT_FILE = 'data/raw/COM_CMB_IQU-sevem_2048_R4.00.fits'
# Punto de Partida: Vértice 647
START_LAT = -41.81
START_LON = 354.38
# Rumbo Inicial de la Rama 3 (detectado por Hydra V6)
INITIAL_BEARING = 204.3 
# HIPÓTESIS FÍSICA: Longitud estimada de la arista del poliedro
# 35-40 grados es un valor típico para modelos de topología cósmica.
TARGET_DISTANCE_DEG = 35.0
# Radio de búsqueda alrededor del punto predicho
SEARCH_RADIUS_DEG = 3.0  

def geodesic_destination(lat1, lon1, bearing, distance_deg):
    """Calcula el punto final tras viajar una distancia con un rumbo fijo."""
    lat1_r, lon1_r = np.radians(lat1), np.radians(lon1)
    bearing_r = np.radians(bearing)
    dist_r = np.radians(distance_deg) # Distancia angular en la esfera unitaria

    lat2_r = np.arcsin(np.sin(lat1_r) * np.cos(dist_r) + 
                       np.cos(lat1_r) * np.sin(dist_r) * np.cos(bearing_r))
    
    lon2_r = lon1_r + np.arctan2(np.sin(bearing_r) * np.sin(dist_r) * np.cos(lat1_r), 
                                 np.cos(dist_r) - np.sin(lat1_r) * np.sin(lat2_r))
    
    return np.degrees(lat2_r), np.degrees(lon2_r)

def scan_area_for_vertex(center_lat, center_lon, radius_deg, map_I, map_P, nside):
    """Escanea un área circular buscando el punto de máxima señal combinada."""
    print(f"🔎 Escaneando área de {radius_deg}º radio alrededor de Lat {center_lat:.2f}, Lon {center_lon:.2f}...")
    
    best_lat, best_lon, max_signal = None, None, 0
    
    # Crear una cuadrícula de búsqueda local
    # Resolución de la cuadrícula (0.1 grados)
    grid_step = 0.1
    lat_range = np.arange(center_lat - radius_deg, center_lat + radius_deg, grid_step)
    lon_range = np.arange(center_lon - radius_deg, center_lon + radius_deg, grid_step)
    
    results_lat, results_lon, results_signal = [], [], []

    for lat in lat_range:
        for lon in lon_range:
            # Verificar si está dentro del radio circular (aprox)
            if (lat - center_lat)**2 + (lon - center_lon)**2 > radius_deg**2:
                continue
                
            pix = hp.ang2pix(nside, np.radians(90 - lat), np.radians(lon))
            # Proxy de señal de vértice: Intensidad * Polarización (busca picos de energía)
            # Un vértice real debería ser un cruce de señales fuertes.
            signal = abs(map_I[pix] * map_P[pix])
            
            results_lat.append(lat)
            results_lon.append(lon)
            results_signal.append(signal)

            if signal > max_signal:
                max_signal = signal
                best_lat, best_lon = lat, lon
                
    return best_lat, best_lon, max_signal, (results_lat, results_lon, results_signal)

def main():
    print("🐶 --- SABUESO V7: INTERCEPTOR DE VÉRTICES ---")
    if not os.path.exists('data/processed'): os.makedirs('data/processed')

    # 1. Cargar Mapas
    print("🛰️ Cargando datos del CMB...")
    try:
        maps = hp.read_map(INPUT_FILE, field=[0,1,2], verbose=False)
        map_I, map_P = maps[0], np.sqrt(maps[1]**2 + maps[2]**2)
        nside = hp.get_nside(map_I)
    except FileNotFoundError:
        print(f"❌ ERROR: No se encuentra el archivo {INPUT_FILE}")
        return

    # 2. Cálculo de Navegación Inercial
    print(f"📐 Calculando destino geodésico...")
    print(f"   Inicio: Vértice 647 (Lat {START_LAT}, Lon {START_LON})")
    print(f"   Rumbo: {INITIAL_BEARING}º | Distancia Objetivo: {TARGET_DISTANCE_DEG}º")
    
    pred_lat, pred_lon = geodesic_destination(START_LAT, START_LON, INITIAL_BEARING, TARGET_DISTANCE_DEG)
    print(f"🎯 Zona de Intercepción Predicha: Lat {pred_lat:.2f}, Lon {pred_lon:.2f}")

    # 3. Escaneo de Área (Radar)
    cand_lat, cand_lon, cand_sig, scan_data = scan_area_for_vertex(pred_lat, pred_lon, SEARCH_RADIUS_DEG, map_I, map_P, nside)
    
    print("\n" + "="*40)
    print("      RESULTADO DE LA INTERCEPCIÓN")
    print("="*40)
    if cand_lat is not None:
        print(f"✅ ¡CANDIDATO A VÉRTICE 648 LOCALIZADO!")
        print(f"📍 Coordenadas: Lat {cand_lat:.4f}º, Lon {cand_lon:.4f}º")
        print(f"📶 Fuerza de Señal: {cand_sig:.2e}")
        print(f"📉 Desviación de la predicción: {np.sqrt((cand_lat-pred_lat)**2 + (cand_lon-pred_lon)**2):.2f}º")
        
        # Guardar el hallazgo
        with open('data/processed/vertex_648_candidate.csv', 'w') as f:
            f.write("lat,lon,signal_strength,notes\n")
            f.write(f"{cand_lat},{cand_lon},{cand_sig},Predicted from V647 bearing {INITIAL_BEARING} dist {TARGET_DISTANCE_DEG}\n")
        print("💾 Coordenadas guardadas en: data/processed/vertex_648_candidate.csv")
        
        # Visualización del Escaneo
        plt.figure(figsize=(10, 8))
        # Puntos del escaneo (mapa de calor)
        plt.scatter(scan_data[1], scan_data[0], c=scan_data[2], cmap='viridis', s=10, alpha=0.5, label='Señal Escaneada')
        plt.colorbar(label='Intensidad de Señal (|T*P|)')
        # Puntos clave
        plt.scatter([pred_lon], [pred_lat], c='red', marker='x', s=150, linewidth=2, label='Predicción Teórica')
        plt.scatter([cand_lon], [cand_lat], c='yellow', marker='*', s=250, edgecolor='black', label='Candidato V648 (Máximo)')
        
        plt.title(f"Escaneo de Intercepción - Zona Vértice 648\n(Distancia {TARGET_DISTANCE_DEG}º desde V647)")
        plt.xlabel("Longitud")
        plt.ylabel("Latitud")
        plt.legend()
        plt.grid(True, alpha=0.3)
        output_img = 'data/processed/vertex_648_scan.png'
        plt.savefig(output_img, dpi=150)
        print(f"🖼️ Mapa del hallazgo generado: {output_img}")
        
    else:
        print("❌ No se encontró una señal clara en el área de búsqueda.")

if __name__ == "__main__":
    main()