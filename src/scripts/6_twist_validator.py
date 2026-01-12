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

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# --- CARGAR LOS DOS HALLAZGOS ---
FILE_ALPHA = 'data/processed/spider_track.csv'
FILE_GHOST = 'data/processed/ghost_face_track.csv'

def get_centroid(df):
    return df['lat'].mean(), df['lon'].mean()

def rotate_points(lats, lons, degrees):
    theta = np.radians(degrees)
    c, s = np.cos(theta), np.sin(theta)
    # Rotación 2D estándar
    new_lons = lons * c - lats * s
    new_lats = lons * s + lats * c
    return new_lats, new_lons

def main():
    print("🧬 TWIST VALIDATOR 2.0: BUSCANDO LA QUIRALIDAD DEL UNIVERSO")
    
    # 1. Cargar y Centrar
    df_alpha = pd.read_csv(FILE_ALPHA)
    df_ghost = pd.read_csv(FILE_GHOST)
    
    # Filtrar vértices/saltos
    pts_alpha = df_alpha[df_alpha['type'].isin(['PATH', 'VERTEX'])]
    pts_ghost = df_ghost[df_ghost['type'].isin(['PATH', 'VERTEX'])]
    
    # Centrar en (0,0)
    lat_a_c, lon_a_c = get_centroid(pts_alpha)
    lat_g_c, lon_g_c = get_centroid(pts_ghost)
    
    alpha_y = pts_alpha['lat'] - lat_a_c
    alpha_x = pts_alpha['lon'] - lon_a_c
    
    ghost_y_raw = pts_ghost['lat'] - lat_g_c
    ghost_x_raw = pts_ghost['lon'] - lon_g_c

    # 2. GENERAR ESCENARIOS
    fig, axes = plt.subplots(2, 2, figsize=(15, 15))
    scenarios = [
        ("+36º (Dextrógiro)", 36, False, axes[0,0]),
        ("-36º (Levógiro)", -36, False, axes[0,1]),
        ("Espejo + 36º", 36, True, axes[1,0]),
        ("Espejo - 36º", -36, True, axes[1,1])
    ]
    
    print("\nEvaluando escenarios...")

    for title, angle, mirror, ax in scenarios:
        # Preparar datos fantasma
        g_x, g_y = ghost_x_raw.copy(), ghost_y_raw.copy()
        
        # Aplicar Espejo si corresponde (Invertir Eje X / Longitud)
        if mirror:
            g_x = -g_x
            
        # Aplicar Rotación
        rot_y, rot_x = rotate_points(g_y, g_x, angle)
        
        # Graficar
        # Cara Alfa (Fija - Cian)
        ax.plot(alpha_x, alpha_y, 'c-', linewidth=3, label='Cara Alfa', alpha=0.7)
        # Cara Fantasma Transformada (Magenta)
        ax.plot(rot_x, rot_y, 'm-', linewidth=3, label='Fantasma')
        
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.axis('equal')
        ax.legend()
        
        # Métrica de similitud visual (Distancia promedio al punto más cercano)
        # Esto es una heurística simple
        from scipy.spatial import cKDTree
        tree = cKDTree(np.column_stack((alpha_x, alpha_y)))
        dists, _ = tree.query(np.column_stack((rot_x, rot_y)))
        score = np.mean(dists)
        ax.set_xlabel(f"Error de Ajuste: {score:.4f}")
        
        print(f"   👉 {title}: Error = {score:.4f}")

    plt.tight_layout()
    output_file = 'data/processed/chirality_check.png'
    plt.savefig(output_file)
    print(f"\n✅ Análisis guardado en: {output_file}")
    print("Busca el gráfico con el 'Error de Ajuste' más bajo (y visualmente más parecido).")

if __name__ == "__main__":
    main()