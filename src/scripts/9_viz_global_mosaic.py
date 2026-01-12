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
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

# --- CONFIGURACIÓN ---
FILE_WIRE = 'data/processed/dodecahedron_wireframe.csv'
FILE_ALPHA = 'data/processed/spider_track_corrected.csv'
FILE_NEIGHBOR = 'data/processed/neighbor1_track.csv'
# Si tienes el archivo fantasma:
FILE_GHOST = 'data/processed/ghost_face_result.csv' 

def latlon2xyz(lat, lon):
    lat_rad = np.radians(lat)
    lon_rad = np.radians(lon)
    x = np.cos(lat_rad) * np.cos(lon_rad)
    y = np.cos(lat_rad) * np.sin(lon_rad)
    z = np.sin(lat_rad)
    return x, y, z

def main():
    print("🌍 GENERANDO MOSAICO TECTÓNICO GLOBAL (viz_global_mosaic.py) ...")
    
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_facecolor('black')
    
    # 1. DIBUJAR ESFERA BASE
    u, v = np.mgrid[0:2*np.pi:50j, 0:np.pi:25j]
    x = 0.98 * np.cos(u) * np.sin(v)
    y = 0.98 * np.sin(u) * np.sin(v)
    z = 0.98 * np.cos(v)
    ax.plot_surface(x, y, z, color='darkblue', alpha=0.1, rstride=1, cstride=1)

    # 2. DIBUJAR JAULA TEÓRICA (Gris o Roja según preferencia)
    # El usuario preguntó por "aristas rojas". Aquí definimos el color de la teoría.
    try:
        df_wire = pd.read_csv(FILE_WIRE)
        print(f"   -> Cargando Jaula Teórica: {len(df_wire)} puntos")
        for edge_id in df_wire['edge_id'].unique():
            segment = df_wire[df_wire['edge_id'] == edge_id]
            xw, yw, zw = latlon2xyz(segment['lat'], segment['lon'])
            # Color GRIS por defecto en notas, pero si quieres resaltar la teoría:
            ax.plot(xw, yw, zw, color='gray', linewidth=0.8, alpha=0.5, linestyle='--')
    except Exception as e:
        print(f"⚠️ No se pudo cargar la jaula teórica: {e}")

    # 3. DIBUJAR CARA ALFA (ROJO) - "Las Aristas Rojas"
    try:
        df_alpha = pd.read_csv(FILE_ALPHA)
        print(f"   -> Cargando Cara Alfa (Rojo): {len(df_alpha)} puntos")
        xa, ya, za = latlon2xyz(df_alpha['lat'], df_alpha['lon'])
        # ESTAS SON LAS ARISTAS ROJAS DE LA REALIDAD
        ax.plot(xa, ya, za, color='red', linewidth=3, label='Cara Alfa (Realidad)')
        
        # Marcar inicio
        ax.scatter(xa.iloc[0], ya.iloc[0], za.iloc[0], color='red', s=50, marker='o')
    except Exception as e:
        print(f"⚠️ No se pudo cargar Cara Alfa: {e}")

    # 4. DIBUJAR VECINO 1 (VERDE)
    try:
        df_neigh = pd.read_csv(FILE_NEIGHBOR)
        print(f"   -> Cargando Vecino 1 (Verde): {len(df_neigh)} puntos")
        xn, yn, zn = latlon2xyz(df_neigh['lat'], df_neigh['lon'])
        ax.plot(xn, yn, zn, color='#00FF00', linewidth=3, label='Vecino 1 (Deformado)')
    except Exception as e:
        print(f"⚠️ No se pudo cargar Vecino 1: {e}")

    # 5. DIBUJAR CARA FANTASMA (MAGENTA)
    try:
        df_ghost = pd.read_csv(FILE_GHOST)
        print(f"   -> Cargando Cara Fantasma (Magenta): {len(df_ghost)} puntos")
        xg, yg, zg = latlon2xyz(df_ghost['lat'], df_ghost['lon'])
        ax.plot(xg, yg, zg, color='magenta', linewidth=3, label='Cara Fantasma (Eco)')
    except:
        # Si no existe, simulamos la antípoda de Alfa para visualización
        if 'df_alpha' in locals():
            print("   -> Generando Fantasma (Simulado Antípoda)")
            xg, yg, zg = latlon2xyz(-df_alpha['lat'], df_alpha['lon'] + 180)
            ax.plot(xg, yg, zg, color='magenta', linewidth=2, linestyle=':', label='Fantasma (Teórico)')

    # Configuración Final
    ax.set_axis_off()
    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.set_zlim([-1, 1])
    plt.legend(loc='lower left', frameon=False, labelcolor='white')
    plt.title("MOSAICO TECTÓNICO GLOBAL\nAristas Rojas (Alfa) vs Teoría (Gris)", color='white', fontsize=14)
    
    # Guardar Vistas
    views = [(30, -10), (90, 0), (30, 170)]
    names = ['front', 'top', 'back']
    
    for (elev, azim), name in zip(views, names):
        ax.view_init(elev=elev, azim=azim)
        outfile = f'data/processed/global_mosaic_{name}.png'
        plt.savefig(outfile, facecolor='black', dpi=150)
        print(f"📸 Guardado: {outfile}")

if __name__ == "__main__":
    main()
