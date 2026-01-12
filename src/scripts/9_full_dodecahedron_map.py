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
from scipy.spatial.transform import Rotation as R
import pandas as pd

# --- CONFIGURACIÓN ---
INPUT_FILE = 'data/raw/COM_CMB_IQU-sevem_2048_R4.00.fits'
ALPHA_LAT = -43.3116
ALPHA_LON = 348.6708
PATCH_SIZE_DEG = 20.0  # Tamaño de la ventana de análisis

def get_dodecahedron_centers():
    # Coordenadas base de los 12 centros de un dodecaedro (radio unitario)
    phi = (1 + np.sqrt(5)) / 2
    verts = []
    # (0, ±1, ±phi)
    for i in [-1, 1]:
        for j in [-1, 1]:
            verts.append([0, i, j*phi])
            verts.append([j*phi, 0, i])
            verts.append([i, j*phi, 0])
    verts = np.array(verts)
    norm = np.linalg.norm(verts, axis=1, keepdims=True)
    return verts / norm

def get_rotation_to_target(source_vec, target_lat, target_lon):
    t_lat_r, t_lon_r = np.radians(target_lat), np.radians(target_lon)
    target_vec = np.array([
        np.cos(t_lat_r) * np.cos(t_lon_r),
        np.cos(t_lat_r) * np.sin(t_lon_r),
        np.sin(t_lat_r)
    ])
    axis = np.cross(source_vec, target_vec)
    axis_norm = np.linalg.norm(axis)
    if axis_norm < 1e-6: return R.identity()
    axis = axis / axis_norm
    angle = np.arccos(np.dot(source_vec, target_vec))
    return R.from_rotvec(axis * angle)

def calculate_texture_score(patch):
    # Índice de Moran simplificado (Proxy de estructura)
    # Buscamos varianza local alta organizada
    vals = patch.flatten()
    vals = vals[vals != 0] # Ignorar máscara
    if len(vals) < 100: return 0
    
    # Z-score
    z = (vals - np.mean(vals)) / np.std(vals)
    
    # Autocorrelación simple (shift 1 pixel)
    # Si es ruido, z * z_shift tiende a 0. Si es estructura, es alto.
    # Usamos un shift arbitrario para detectar "manchas"
    score = np.mean(z * np.roll(z, 1)) 
    return score * 100 # Escalar para legibilidad

def main():
    print("🌐 INICIANDO PROTOCOLO DODECA-SCANNER (MAPA COMPLETO) 🌐")
    
    # 1. Cargar Datos
    # Suppress verbose output if possible or ignore warning
    maps = hp.read_map(INPUT_FILE, field=[0,1,2])
    map_comb = maps[0] * np.sqrt(maps[1]**2 + maps[2]**2)
    nside = hp.get_nside(maps[0])
    
    # Máscara Galáctica
    npix = hp.nside2npix(nside)
    theta, _ = hp.pix2ang(nside, np.arange(npix))
    lat_gal = 90 - np.degrees(theta)
    mask = (np.abs(lat_gal) > 20.0)
    map_comb[~mask] = 0.0

    # 2. Alinear Dodecaedro con Cara Alfa
    base_centers = get_dodecahedron_centers()
    # Asumimos que el vertice 0 es nuestra Alfa
    rot_matrix = get_rotation_to_target(base_centers[0], ALPHA_LAT, ALPHA_LON)
    real_centers_vec = rot_matrix.apply(base_centers)
    
    # Convertir a Lat/Lon
    lats = np.degrees(np.arcsin(real_centers_vec[:, 2]))
    lons = np.degrees(np.arctan2(real_centers_vec[:, 1], real_centers_vec[:, 0]))

    # 3. Escanear las 12 Caras
    results = []
    print(f"   ... Escaneando las 12 caras teóricas ...")
    
    for i in range(12):
        lat, lon = lats[i], lons[i]
        
        # Extraer Patch
        patch = hp.gnomview(map_comb, rot=[lon, lat, 0], xsize=100, ysize=100, 
                            reso=12, return_projected_map=True, no_plot=True)
        patch = patch.filled(0)
        
        # Calcular Score
        score = calculate_texture_score(patch)
        
        # Clasificación
        status = "🟢 SÓLIDO" if score > 50 else "🟡 DÉBIL" if score > 20 else "🔴 RUIDO"
        
        results.append({
            'id': i+1, 'lat': lat, 'lon': lon, 'score': score, 'status': status
        })
        print(f"   -> Cara {i+1:02d}: Lat {lat:6.1f}, Lon {lon:6.1f} | Score: {score:6.2f} | {status}")

    # 4. Generar Mapa Visual
    print("   ... Dibujando el Mapa del Universo ...")
    plt.figure(figsize=(15, 10))
    
    # Fondo: Mapa Mollweide del cielo completo (baja res para fondo)
    hp.mollview(map_comb, title="EL MAPA DEL UNIVERSO DODECAÉDRICO\n(12 Centros de Cara Identificados)", 
                cmap='Greys', min=0, max=np.max(map_comb)/2, cbar=False, hold=True)
    
    # Superponer los 12 centros
    for res in results:
        # Convertir lat/lon a coordenadas proj Mollweide
        # healpy projplot usa (lon, lat)
        hp.projtext(res['lon'], res['lat'], f"{res['id']}", lonlat=True, 
                    fontsize=12, fontweight='bold', color='white', ha='center')
        
        # Círculo indicando el score (Color = Intensidad)
        color = 'lime' if res['score'] > 50 else 'gold' if res['score'] > 20 else 'red'
        size = max(50, res['score'] * 2)
        hp.projscatter(res['lon'], res['lat'], lonlat=True, s=size, 
                       edgecolor='white', facecolor=color, alpha=0.7, zorder=5)

    plt.savefig('data/processed/full_dodecahedron_map.png')
    print("✨ ¡MAPA COMPLETADO! Abre: data/processed/full_dodecahedron_map.png")
    
    # Guardar CSV con las coordenadas para el futuro
    pd.DataFrame(results).to_csv('data/processed/dodecahedron_faces_coordinates.csv', index=False)

if __name__ == "__main__":
    main()
