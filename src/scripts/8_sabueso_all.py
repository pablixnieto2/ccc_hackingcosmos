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
from scipy.spatial import cKDTree

# --- CONFIGURACIÓN ---
INPUT_FILE = 'data/raw/COM_CMB_IQU-sevem_2048_R4.00.fits'
ALPHA_LAT = -43.3116
ALPHA_LON = 348.6708
PATCH_RADIUS = 15.0 # Grados

def get_dodecahedron_neighbors():
    # En un dodecaedro, una cara tiene 5 vecinos adyacentes.
    # Si la Cara 0 está en el Polo Norte (0,0,1), los vecinos están en un anillo
    # a latitud approx 26.57 grados (o 63.43 desde el polo).
    phi = (1 + np.sqrt(5)) / 2
    # Ángulo diedro del dodecaedro es ~116.56 grados.
    # Distancia entre centros de caras adyacentes es ~63.43 grados.
    
    # Generamos los 5 centros vecinos relativos a un centro en el Polo Z
    # Latitud de los vecinos desde el ecuador:
    lat_neighbor = np.degrees(np.arctan(0.5)) # ~26.565 grados
    
    neighbors = []
    # Los 5 vecinos están rotados 72 grados entre sí en longitud
    for i in range(5):
        lon = i * 72.0
        # Convertir a vector
        theta = np.radians(90 - lat_neighbor)
        phi_ang = np.radians(lon)
        x = np.sin(theta) * np.cos(phi_ang)
        y = np.sin(theta) * np.sin(phi_ang)
        z = np.cos(theta)
        neighbors.append([x, y, z])
    
    return np.array(neighbors)

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

def calculate_local_moran(patch_data):
    # Calcula el Moran's I rápido para una imagen/patch 2D
    # Aplanamos
    values = patch_data.flatten()
    # Eliminamos nulos (máscara)
    values = values[values != 0]
    
    if len(values) < 10: return 0.0
    
    # Normalización Z-score
    mean = np.mean(values)
    std = np.std(values)
    if std == 0: return 0.0
    z = (values - mean) / std
    
    # Moran simplificado: Autocorrelación de valores altos con altos
    # Asumimos que los vecinos en el array plano son representativos (heurística rápida)
    # Para hacerlo bien, necesitaríamos matriz de pesos, pero para comparar parches
    # la varianza espacial (std) combinada con la magnitud suele ser un buen proxy.
    
    # MEJOR: Usar la desviación estándar espacial (Textura)
    # Si hay estructura, hay contraste (bordes).
    # Moran I ~ Autocorrelación.
    
    # Vamos a usar una métrica híbrida: "Energía Estructural"
    # Suma de (valor * valor_vecino)
    # Usamos roll para simular vecinos
    z_right = np.roll(z, 1)
    z_down = np.roll(z, int(np.sqrt(len(z))))
    
    moran_proxy = np.sum(z * z_right) + np.sum(z * z_down)
    return moran_proxy / len(z)

def main():
    print("🕵️ SABUESO MORAN SCANNER: BUSCANDO LA CARA VECINA MÁS SÓLIDA")
    
    # 1. Cargar Mapa
    maps = hp.read_map(INPUT_FILE, field=[0,1,2], verbose=False)
    map_comb = maps[0] * np.sqrt(maps[1]**2 + maps[2]**2)
    nside = hp.get_nside(maps[0])
    
    # Máscara Galáctica Rápida
    npix = hp.nside2npix(nside)
    theta, _ = hp.pix2ang(nside, np.arange(npix))
    lat_gal = 90 - np.degrees(theta)
    mask = (np.abs(lat_gal) > 20.0)
    map_comb[~mask] = 0.0

    # 2. Generar Candidatos (Los 5 Vecinos de Alfa)
    print("   ... Calculando geometría local ...")
    # Vector base (Polo Norte)
    base_pole = np.array([0,0,1])
    # Rotación para llevar el Polo a nuestra Cara Alfa
    rot_alpha = get_rotation_to_target(base_pole, ALPHA_LAT, ALPHA_LON)
    
    # Obtener vecinos base y rotarlos
    neighbors_base = get_dodecahedron_neighbors()
    real_neighbors = rot_alpha.apply(neighbors_base)
    
    # 3. Escanear cada vecino
    print("   ... Escaneando texturas topológicas ...")
    
    results = []
    
    plt.figure(figsize=(15, 5))
    
    for i in range(5):
        vec = real_neighbors[i]
        lat = np.degrees(np.arcsin(vec[2]))
        lon = np.degrees(np.arctan2(vec[1], vec[0]))
        
        # Extraer Patch (Gnomview)
        patch = hp.gnomview(map_comb, rot=[lon, lat, 0], xsize=100, ysize=100, 
                            reso=15, return_projected_map=True, no_plot=True)
        patch = patch.filled(0)
        
        # Calcular Moran Score (Estructura)
        score = calculate_local_moran(patch)
        
        results.append({'id': i+1, 'lat': lat, 'lon': lon, 'score': score, 'patch': patch})
        print(f"      -> Vecino {i+1}: Lat {lat:.1f}, Lon {lon:.1f} | Moran Score: {score:.4f}")
        
        # Visualizar
        plt.subplot(1, 5, i+1)
        plt.imshow(patch, cmap='inferno', origin='lower')
        plt.title(f"Vecino {i+1}\nScore: {score:.2f}")
        plt.axis('off')

    plt.tight_layout()
    plt.savefig('data/processed/moran_neighbors_scan.png')
    
    # 4. Veredicto
    best = max(results, key=lambda x: x['score'])
    print("\n" + "="*50)
    print(f"🏆 GANADOR: VECINO {best['id']}")
    print(f"📍 Coordenadas Objetivo: Lat {best['lat']:.4f}, Lon {best['lon']:.4f}")
    print(f"📈 Estructura (Moran): {best['score']:.4f}")
    print("="*50)
    print("👉 Este es tu próximo objetivo. Aquí es donde la 'costura' del universo es más visible.")

if __name__ == "__main__":
    main()