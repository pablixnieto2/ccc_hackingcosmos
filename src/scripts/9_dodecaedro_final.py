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
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial.transform import Rotation as R
import pandas as pd
import os

# --- CONFIGURACIÓN ---
ALPHA_LAT = -43.3116
ALPHA_LON = 348.6708
CSV_FILE = 'data/processed/dodecahedron_faces_coordinates.csv'

def get_face_status():
    # Leer el estado de las caras desde el CSV generado por el escáner
    if not os.path.exists(CSV_FILE):
        print("⚠️ No se encontró el CSV. Usando valores por defecto.")
        return [1] * 12 # Todos verdes por defecto
    
    df = pd.read_csv(CSV_FILE)
    # Asumimos que el CSV está ordenado por ID 1..12
    # Convertir 'status' a 1 (Sólido/Débil) o 0 (Ruido)
    status_list = []
    for _, row in df.iterrows():
        s = row['status']
        if "RUIDO" in s:
            status_list.append(0)
        else:
            status_list.append(1)
    return status_list

def get_dodecahedron_geometry():
    # El número áureo
    phi = (1 + np.sqrt(5)) / 2
    
    # Vértices de un dodecaedro (20 vértices)
    # Grupo 1: (±1, ±1, ±1)
    vertices = []
    for i in [-1, 1]:
        for j in [-1, 1]:
            for k in [-1, 1]:
                vertices.append([i, j, k])
    
    # Grupo 2: (0, ±phi, ±1/phi) y permutaciones cíclicas
    inv_phi = 1 / phi
    for i in [-1, 1]:
        for j in [-1, 1]:
            vertices.append([0, i*phi, j*inv_phi])
            vertices.append([j*inv_phi, 0, i*phi])
            vertices.append([i*phi, j*inv_phi, 0])
            
    vertices = np.array(vertices)
    # Normalizar al radio unitario
    vertices = vertices / np.linalg.norm(vertices[0])
    
    return vertices

def get_faces(verts):
    # Algoritmo simple para agrupar vértices en caras pentagonales
    faces = []
    centers = []
    
    # Buscamos los centros de las caras (que son vértices de un icosaedro dual)
    phi = (1 + np.sqrt(5)) / 2
    ico_verts = []
    for i in [-1, 1]:
        for j in [-1, 1]:
            ico_verts.append([0, i, j*phi])
            ico_verts.append([j*phi, 0, i])
            ico_verts.append([i, j*phi, 0])
    ico_verts = np.array(ico_verts)
    ico_verts /= np.linalg.norm(ico_verts, axis=1, keepdims=True)
    
    # Para cada centro teórico, buscamos los 5 vértices del dodecaedro más cercanos
    for center in ico_verts:
        dists = np.linalg.norm(verts - center, axis=1)
        # Los 5 más cercanos forman la cara
        idx = np.argsort(dists)[:5]
        face_verts = verts[idx]
        
        # Ordenar vértices angularmente
        v_centered = face_verts - center
        z_axis = center
        x_axis = np.cross(np.array([0,0,1]), z_axis)
        if np.linalg.norm(x_axis) < 0.1: x_axis = np.array([1,0,0])
        x_axis /= np.linalg.norm(x_axis)
        y_axis = np.cross(z_axis, x_axis)
        
        angles = np.arctan2(np.dot(v_centered, y_axis), np.dot(v_centered, x_axis))
        sort_order = np.argsort(angles)
        
        faces.append(face_verts[sort_order])
        centers.append(center)
        
    return faces, np.array(centers)

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

def main():
    print("💎 GENERANDO HOLOGRAMA DEL UNIVERSO DODECAÉDRICO...")
    
    # 0. Obtener Estado de Caras
    face_status = get_face_status()
    print(f"   -> Estado de caras cargado: {face_status}")
    
    # 1. Geometría Base
    verts = get_dodecahedron_geometry()
    faces_raw, centers_raw = get_faces(verts)
    
    # 2. Alinear con Cara Alfa
    # Importante: Asegurar que el orden de faces_raw coincida con el orden del escáner (1..12)
    # El escáner genera coordenadas basadas en rotar los centros base.
    # Aquí rotamos todo el objeto para que la Cara 0 (la primera generada) apunte a Alfa.
    rot = get_rotation_to_target(centers_raw[0], ALPHA_LAT, ALPHA_LON)
    
    rotated_faces = []
    for face in faces_raw:
        rotated_faces.append(rot.apply(face))
        
    # 3. Visualización 3D
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_facecolor('black')
    
    print("   ... Construyendo estructura de cristal ...")
    
    for i, poly in enumerate(rotated_faces):
        # Determinar color
        status = face_status[i] if i < len(face_status) else 1
        
        if status == 1:
            # Sólido (Verde)
            color = '#00FF00' # Verde Matrix
            edge = 'lime'
            alpha = 0.4
            lw = 1
        else:
            # Ruido (Rojo - Oculto por Galaxia)
            # AUMENTAMOS VISIBILIDAD PARA QUE SE VEA BIEN
            color = '#550000' # Rojo Oscuro
            edge = 'red'      # Arista Roja Brillante
            alpha = 0.5       # Más opaco para que se note la "tapa"
            lw = 2            # Arista más gruesa
        
        # Crear polígono
        poly3d = Poly3DCollection([poly], alpha=alpha, linewidths=lw)
        poly3d.set_facecolor(color)
        poly3d.set_edgecolor(edge)
        ax.add_collection3d(poly3d)
        
        # Etiqueta de número de cara
        center = np.mean(poly, axis=0) * 1.15
        text_color = 'white' if status == 1 else 'red'
        ax.text(center[0], center[1], center[2], f"{i+1}", color=text_color, fontsize=12, fontweight='bold', ha='center')

    # Configuración de cámara
    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.set_zlim([-1, 1])
    ax.set_axis_off()
    
    solid_count = sum(face_status)
    title = f"EL UNIVERSO CONFIRMADO\n{solid_count}/12 Caras Detectadas (Verde)\n{12-solid_count} Ocultas por Galaxia (Rojo)"
    plt.title(title, color='white', fontsize=14)
    
    # Guardar vistas
    views = [(30, 0), (30, 90), (30, 180), (90, 0)]
    names = ['front', 'side', 'back', 'top']
    
    for (elev, azim), name in zip(views, names):
        ax.view_init(elev=elev, azim=azim)
        filename = f'data/processed/hologram_{name}.png'
        plt.savefig(filename, facecolor='black', dpi=150)
        print(f"📸 Holograma capturado: {filename}")

if __name__ == "__main__":
    main()
