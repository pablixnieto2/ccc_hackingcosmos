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
from scipy.ndimage import rotate

# --- CONFIGURACIÓN ---
INPUT_FILE = 'data/raw/COM_CMB_IQU-sevem_2048_R4.00.fits'
ALPHA_LAT = -43.3116
ALPHA_LON = 348.6708
STACK_RADIUS_DEG = 15.0
IMG_SIZE = 200

def get_icosahedron_vertices():
    phi = (1 + np.sqrt(5)) / 2
    verts = []
    for i in [-1, 1]:
        for j in [-1, 1]:
            verts.append([0, i, j*phi])
            verts.append([j*phi, 0, i])
            verts.append([i, j*phi, 0])
    verts = np.array(verts)
    return verts / np.linalg.norm(verts, axis=1, keepdims=True)

def cartesian_to_spherical(xyz):
    lat = np.degrees(np.arcsin(xyz[:, 2]))
    lon = np.degrees(np.arctan2(xyz[:, 1], xyz[:, 0]))
    return lat, lon

def get_rotation_to_target(source_vec, target_lat, target_lon):
    t_lat_r, t_lon_r = np.radians(target_lat), np.radians(target_lon)
    target_vec = np.array([np.cos(t_lat_r) * np.cos(t_lon_r), np.cos(t_lat_r) * np.sin(t_lon_r), np.sin(t_lat_r)])
    axis = np.cross(source_vec, target_vec)
    axis_norm = np.linalg.norm(axis)
    if axis_norm < 1e-6: return R.identity()
    return R.from_rotvec((axis / axis_norm) * np.arccos(np.dot(source_vec, target_vec)))

def main():
    print("🧹 STACKING RITUAL V2: LIMPIEZA GALÁCTICA 🧹")
    
    # 1. Cargar Mapas
    print("   ... Cargando mapa ...")
    maps = hp.read_map(INPUT_FILE, field=[0,1,2], verbose=False)
    # Combinar I, Q, U
    map_comb = maps[0] * np.sqrt(maps[1]**2 + maps[2]**2)
    nside = hp.get_nside(maps[0])
    
    # --- MÁSCARA GALÁCTICA ---
    print("   ... Aplicando escudo anti-galáctico ...")
    # Creamos una máscara simple: Cortamos +/- 20 grados del ecuador galáctico
    # En coordenadas galácticas, el ecuador es lat=0.
    # El mapa ya viene en coordenadas galácticas usualmente.
    
    # Convertir índices de píxel a coordenadas para ver dónde está cada uno
    theta, phi = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
    lat_gal = 90 - np.degrees(theta)
    
    # Definir zona de exclusión (La Vía Láctea es ruidosa)
    # Enmascaramos todo lo que esté entre -20 y +20 grados de latitud
    galactic_mask = (np.abs(lat_gal) > 20.0) 
    
    # Aplicar máscara (ponemos a 0 o NaN lo sucio)
    map_clean = map_comb.copy()
    map_clean[~galactic_mask] = 0 # Ponemos a 0 la galaxia para que no sume ruido
    
    # 2. Geometría
    base_centers = get_icosahedron_vertices()
    alpha_rot = get_rotation_to_target(base_centers[0], ALPHA_LAT, ALPHA_LON)
    real_centers_vec = alpha_rot.apply(base_centers)
    real_lats, real_lons = cartesian_to_spherical(real_centers_vec)

    print(f"   ... Procesando 12 caras con máscara activa ...")
    
    stack_accumulator = np.zeros((IMG_SIZE, IMG_SIZE))
    valid_counts = np.zeros((IMG_SIZE, IMG_SIZE)) # Para promediar solo píxeles válidos
    
    reso_arcmin = (STACK_RADIUS_DEG * 2 * 60) / IMG_SIZE

    for i in range(12):
        lat = real_lats[i]
        lon = real_lons[i]
        
        # Proyectar
        try:
            patch = hp.gnomview(map_clean, rot=[lon, lat, 0], xsize=IMG_SIZE, ysize=IMG_SIZE, 
                                reso=reso_arcmin, return_projected_map=True, no_plot=True)
        except:
            continue # Skip si falla
            
        if hasattr(patch, 'filled'): patch = patch.filled(0)
        else: patch = np.nan_to_num(patch)

        # Quiralidad
        dot_prod = np.dot(real_centers_vec[i], real_centers_vec[0])
        if dot_prod < 0: 
             patch = np.fliplr(patch)
             patch = rotate(patch, 36, reshape=False, mode='nearest')
             tag = "[ANTIPODA]"
        else:
            tag = "[FRONTAL]"
            
        # Acumular (Solo si el pixel no es 0 absoluto de la máscara)
        # Creamos una máscara local para este patch
        local_mask = (np.abs(patch) > 1e-9).astype(float)
        
        stack_accumulator += patch
        valid_counts += local_mask # Contamos cuántas caras contribuyen a cada pixel
        
        print(f"      -> Cara {i+1} asimilada {tag}")

    # Normalizar (Promedio inteligente)
    # Evitar división por cero
    valid_counts[valid_counts == 0] = 1
    final_image = stack_accumulator / valid_counts
    
    # VISUALIZACIÓN
    plt.figure(figsize=(10, 10))
    extent = [-STACK_RADIUS_DEG, STACK_RADIUS_DEG, -STACK_RADIUS_DEG, STACK_RADIUS_DEG]
    
    # Usamos vmin/vmax para recortar picos extremos y ver el detalle tenue
    limit = np.percentile(np.abs(final_image), 98) 
    
    plt.imshow(final_image, cmap='RdBu_r', origin='lower', extent=extent, vmin=-limit, vmax=limit)
    plt.colorbar(label='Señal Promedio Filtrada')
    plt.title(f"THE COSMIC DODECAHEDRON (CLEAN)\nStacking con Máscara Galáctica")

    # Guía Pentagonal
    angles = np.linspace(0, 2*np.pi, 6) + np.radians(90)
    r_pent = 6.0 
    x_p = r_pent * np.cos(angles)
    y_p = r_pent * np.sin(angles)
    plt.plot(x_p, y_p, 'y--', linewidth=2, label='Pentágono Teórico')
    
    plt.legend()
    plt.savefig('data/processed/final_dodecahedron_stack_clean.png', dpi=150)
    print("✨ ¡IMAGEN LIMPIA GENERADA! Abre: data/processed/final_dodecahedron_stack_clean.png")

if __name__ == "__main__":
    main()