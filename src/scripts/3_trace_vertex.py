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

# --- CONFIGURACIÓN "MODO MICROSCOPIO" ---
INPUT_FILE = 'data/raw/COM_CMB_IQU-sevem_2048_R4.00.fits'
OUTPUT_FILE = 'data/processed/vertex_trace_647.csv'

# OBJETIVO: LA CICATRIZ 647
TARGET_LAT = -41.81
TARGET_LON = 354.38  

# ZONA DE OPERACIONES
ROI_RADIUS = 10.0    # Radio de búsqueda en grados

# ¡CAMBIO 1! DENSIDAD EXTREMA
# NSIDE 256 = Puntos separados por ~0.2 grados.
NSIDE_TRACE = 256     

# ¡CAMBIO 2! SIN HUECOS ANGULARES
# Probamos CADA grado. Si está en 33º, lo veremos.
ANGLE_STEP = 1       

# HERRAMIENTAS
LINE_WIDTH = 0.5
LINE_LENGTHS = [4.0, 8.0] 

# UMBRAL DE SENSIBILIDAD
# Lo mantenemos bajo para ver hasta los "ecos" más débiles
CORR_THRESHOLD = 0.08  

def get_rotated_rect_pixels(nside, center_vec, length_deg, width_deg, angle_deg):
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
    best_corr = 0
    best_line = None
    
    # BARRIDO COMPLETO DE ÁNGULOS (GRADO A GRADO)
    for angle in range(0, 180, ANGLE_STEP):
        for length in LINE_LENGTHS:
            pixels = get_rotated_rect_pixels(nside, center_vec, length, LINE_WIDTH, angle)
            if len(pixels) < 20: continue 
            
            vals_I = map_I[pixels]
            vals_Q = map_Q[pixels]
            vals_U = map_U[pixels]
            vals_P = np.sqrt(vals_Q**2 + vals_U**2)
            
            if len(vals_I) > 2 and np.std(vals_I) > 1e-6 and np.std(vals_P) > 1e-6:
                corr, _ = pearsonr(vals_I, vals_P)
                
                # Guardamos el mejor candidato en este píxel
                if abs(corr) > abs(best_corr):
                    best_corr = corr
                    best_line = {
                        'center_idx': idx,
                        'angle': angle,
                        'length': length,
                        'corr_IP': corr
                    }
    
    if best_line and abs(best_corr) > CORR_THRESHOLD:
        return best_line
    return None

def main():
    print(f"--- RASTREO DE ALTA RESOLUCIÓN (MICRO-GRID) ---")
    print(f"Objetivo: Lat {TARGET_LAT}, Lon {TARGET_LON} | NSIDE: {NSIDE_TRACE} | Angle Step: {ANGLE_STEP}º")
    
    # 1. Cargar Datos
    print("Cargando mapas...")
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
    nside_map = hp.get_nside(map_I)

    # 2. Definir Zona de Operaciones
    theta_rad = np.radians(90.0 - TARGET_LAT)
    phi_rad = np.radians(TARGET_LON)
    target_vec = hp.ang2vec(theta_rad, phi_rad)
    
    roi_pixels = hp.query_disc(NSIDE_TRACE, target_vec, np.radians(ROI_RADIUS))
    roi_vectors = [hp.pix2vec(NSIDE_TRACE, p) for p in roi_pixels]
    
    print(f"Analizando {len(roi_vectors)} puntos de la red local...")

    # 3. Ejecutar Rastreo Paralelo
    tasks = [(roi_pixels[i], vec, map_I, map_Q, map_U, nside_map) for i, vec in enumerate(roi_vectors)]
    
    traces = []
    print("Iniciando barrido angular grado a grado...")
    
    with Pool(processes=cpu_count()) as pool:
        # Usamos chunksize pequeño para actualizar más a menudo
        for i, res in enumerate(pool.imap_unordered(process_point, tasks, chunksize=5)):
            if res:
                traces.append(res)
            if i % 100 == 0:
                sys.stdout.write(f"\rProgreso: {i}/{len(roi_vectors)} puntos analizados. Detectados: {len(traces)}")
                sys.stdout.flush()
                
    # 4. Guardar Resultados
    print("\n")
    if traces:
        df = pd.DataFrame(traces)
        df.to_csv(OUTPUT_FILE, index=False)
        print(f"¡HECHO! Se han cartografiado {len(df)} segmentos de fractura.")
        print(f"Datos en: {OUTPUT_FILE}")
    else:
        print("Resultado vacío. La zona parece limpia de correlaciones lineales > 0.08.")

if __name__ == "__main__":
    main()