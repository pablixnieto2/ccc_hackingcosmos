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
import sys
from multiprocessing import Pool, cpu_count

# --- CONFIGURACIÓN ---
INPUT_FILE = 'data/raw/COM_CMB_IQU-sevem_2048_R4.00.fits'
OUTPUT_FILE = 'data/processed/line_metrics.csv'

NSIDE_SCAN = 8       # Resolución de puntos centrales
LINE_LENGTH = 8.0    # Largo de la línea
LINE_WIDTH = 0.5     # Grosor de la línea

# MEJORA 1: PRECISIÓN DE FRANCOTIRADOR
# Bajamos de 15º a 2º. Ahora no se nos escapa nada entre medias.
ANGLE_STEP = 2       

# MEJORA 2: FILTRO ANTI-GALAXIA
# Si la latitud es menor a esto, ni nos molestamos en taladrar. Ahorra tiempo y falsas ilusiones.
MIN_LAT_FILTER = 30.0 

def get_rotated_rect_pixels(nside, center_vec, length_deg, width_deg, angle_deg):
    """Crea un rectángulo fino rotado sobre la esfera."""
    z_axis = center_vec / np.linalg.norm(center_vec)
    aux = np.array([0, 0, 1]) if abs(z_axis[2]) < 0.9 else np.array([0, 1, 0])
    x_axis = np.cross(aux, z_axis)
    x_axis /= np.linalg.norm(x_axis)
    y_axis = np.cross(z_axis, x_axis)
    
    L = np.radians(length_deg) / 2.0 
    W = np.radians(width_deg) / 2.0  
    rad_angle = np.radians(angle_deg)
    
    c, s = np.cos(rad_angle), np.sin(rad_angle)
    rot_matrix = np.array([[c, -s], [s, c]])
    
    corners_2d = np.array([[L, W], [-L, W], [-L, -W], [L, -W]])
    
    sphere_vertices = []
    for corner in corners_2d:
        rotated = np.dot(rot_matrix, corner)
        vertex_vec = z_axis + (x_axis * rotated[0]) + (y_axis * rotated[1])
        vertex_vec /= np.linalg.norm(vertex_vec) 
        sphere_vertices.append(vertex_vec)
        
    return hp.query_polygon(nside, np.array(sphere_vertices))

def process_point(args):
    idx, center_vec, map_I, map_Q, map_U, nside = args
    results = []
    
    # --- FILTRO DE LATITUD (OPTIMIZACIÓN) ---
    # Convertimos vector a latitud para ver si estamos en zona segura
    # healpy usa colatitude (theta), lat = 90 - theta
    cx, cy, cz = center_vec
    radius = np.sqrt(cx*cx + cy*cy + cz*cz)
    theta = np.arccos(cz/radius)
    lat_deg = 90.0 - np.degrees(theta)
    
    # Si estamos en la Galaxia (entre -30 y 30), abortamos misión en este punto.
    if abs(lat_deg) < MIN_LAT_FILTER:
        return []

    # --- BARRIDO DE ALTA PRECISIÓN ---
    # Probamos cada 2 grados. Si la línea está en 1.5º, el escaneo de 2.0º la tocará.
    for angle in range(0, 180, ANGLE_STEP):
        
        pixels = get_rotated_rect_pixels(nside, center_vec, LINE_LENGTH, LINE_WIDTH, angle)
        
        if len(pixels) < 20: continue 
        
        vals_I = map_I[pixels]
        vals_Q = map_Q[pixels]
        vals_U = map_U[pixels]
        vals_P = np.sqrt(vals_Q**2 + vals_U**2)
        
        # Filtro básico de varianza para evitar errores numéricos
        if len(vals_I) > 2 and np.std(vals_I) > 0 and np.std(vals_P) > 0:
            corr, _ = pearsonr(vals_I, vals_P)
            
            # Guardamos candidatos fuertes
            if abs(corr) > 0.15: 
                results.append({
                    'center_idx': idx,
                    'lat': lat_deg,      # Guardamos Lat para verificar luego
                    'angle': angle,
                    'corr_IP': corr,
                    'pixel_count': len(pixels)
                })
                
    return results

def main():
    print("--- LINE HUNTER V3 (PRECISION + FILTERS) ---")
    print(f"Angle Step: {ANGLE_STEP} deg | Galactic Filter: |Lat| > {MIN_LAT_FILTER}")
    
    print("Loading Data...")
    try:
        maps = hp.read_map(INPUT_FILE, field=None, hdu=1, verbose=False, memmap=True)
        if len(maps) == 1 or maps.ndim == 1:
            map_I = maps
            maps_pol = hp.read_map(INPUT_FILE, field=None, hdu=2, verbose=False, memmap=True)
            map_Q, map_U = maps_pol[0], maps_pol[1]
        else:
            map_I, map_Q, map_U = maps[0], maps[1], maps[2]
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
        
    nside = hp.get_nside(map_I)
    
    npix_scan = hp.nside2npix(NSIDE_SCAN)
    scan_vectors = [hp.pix2vec(NSIDE_SCAN, i) for i in range(npix_scan)]
    print(f"Scanning {len(scan_vectors)} points...")
    
    tasks = [(i, vec, map_I, map_Q, map_U, nside) for i, vec in enumerate(scan_vectors)]
    
    all_hits = []
    print("Drilling Deep Sky...")
    
    with Pool(processes=cpu_count()) as pool:
        for res in pool.imap_unordered(process_point, tasks, chunksize=20):
            if res:
                all_hits.extend(res)
                
    if all_hits:
        df = pd.DataFrame(all_hits)
        df.to_csv(OUTPUT_FILE, index=False)
        print(f"SUCCESS: Found {len(df)} High-Confidence Lines.")
        print("\nTOP 5 DEEP SKY LINES:")
        # Ordenamos por fuerza absoluta de correlación
        df['abs_corr'] = df['corr_IP'].abs()
        print(df.sort_values(by='abs_corr', ascending=False).head(5))
    else:
        print("No significant lines found in Deep Sky.")

if __name__ == "__main__":
    main()