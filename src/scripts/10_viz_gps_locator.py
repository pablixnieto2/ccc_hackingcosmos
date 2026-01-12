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
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial.transform import Rotation as R
import random

# --- DATOS DEL DESCUBRIMIENTO ---
# 10 Caras detectadas (Verde), 2 Ocultas por Galaxia (Rojo)
FACE_STATUS = [1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1] 
ALPHA_LAT = -43.3116
ALPHA_LON = 348.6708
OBSERVER_POS = [0.12, -0.04, 0.28] 

def neon_glow_line(ax, p1, p2, color, core_width=1):
    """Simula efecto neón dibujando múltiples líneas con decreciente opacidad"""
    # Brillo exterior (Halo)
    ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], 
            color=color, linewidth=core_width*6, alpha=0.05)
    ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], 
            color=color, linewidth=core_width*3, alpha=0.15)
    # Núcleo brillante
    ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], 
            color='white', linewidth=core_width*0.5, alpha=0.9)
    # Color principal
    ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], 
            color=color, linewidth=core_width, alpha=1.0)

def get_dodecahedron_geometry():
    phi = (1 + np.sqrt(5)) / 2
    vertices = []
    for i in [-1, 1]:
        for j in [-1, 1]:
            for k in [-1, 1]:
                vertices.append([i, j, k])
    inv_phi = 1 / phi
    for i in [-1, 1]:
        for j in [-1, 1]:
            vertices.append([0, i*phi, j*inv_phi])
            vertices.append([j*inv_phi, 0, i*phi])
            vertices.append([i*phi, j*inv_phi, 0])
    vertices = np.array(vertices)
    vertices = vertices / np.linalg.norm(vertices[0])
    return vertices

def get_faces(verts):
    faces = []
    centers = []
    phi = (1 + np.sqrt(5)) / 2
    ico_verts = []
    for i in [-1, 1]:
        for j in [-1, 1]:
            ico_verts.append([0, i, j*phi])
            ico_verts.append([j*phi, 0, i])
            ico_verts.append([i, j*phi, 0])
    ico_verts = np.array(ico_verts)
    ico_verts /= np.linalg.norm(ico_verts, axis=1, keepdims=True)
    
    for center in ico_verts:
        dists = np.linalg.norm(verts - center, axis=1)
        idx = np.argsort(dists)[:5]
        face_verts = verts[idx]
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
    print("🎬 INICIANDO MOTOR GRÁFICO 'CINEMA'...")
    
    verts = get_dodecahedron_geometry()
    faces_raw, centers_raw = get_faces(verts)
    rot = get_rotation_to_target(centers_raw[0], ALPHA_LAT, ALPHA_LON)
    
    rotated_faces = []
    for face in faces_raw:
        rotated_faces.append(rot.apply(face))
        
    # --- CONFIGURACIÓN VISUAL ---
    plt.style.use('dark_background')
    fig = plt.figure(figsize=(12, 12)) 
    ax = fig.add_subplot(111, projection='3d')
    ax.set_facecolor('black')
    

    # 1. DIBUJAR CARAS CON NEÓN
    print("   ⚡ Activando celdas de energía...")
    for i, poly in enumerate(rotated_faces):
        status = FACE_STATUS[i] if i < len(FACE_STATUS) else 1
        
        if status == 1:
            edge_color = '#00FFFF' # Cyan eléctrico
            face_color = '#00FFaa' 
            glow_strength = 0.8
        else:
            edge_color = '#FF0033' # Rojo Neón
            face_color = '#330000'
            glow_strength = 0.3

        poly3d = Poly3DCollection([poly], alpha=0.1 * glow_strength)
        poly3d.set_facecolor(face_color)
        ax.add_collection3d(poly3d)
        
        for j in range(len(poly)):
            p1 = poly[j]
            p2 = poly[(j+1)%len(poly)]
            neon_glow_line(ax, p1, p2, edge_color, core_width=1.5)

    # 2. DIBUJAR AL OBSERVADOR Y CENTRO
    print(f"   📍 Posicionando baliza humana...")
    
    # Centro del Universo (Referencia)
    ax.scatter(0, 0, 0, color='gray', s=50, marker='+', alpha=0.5, label='CENTRO ABSOLUTO')
    
    # Nosotros (Observador)
    ax.scatter(OBSERVER_POS[0], OBSERVER_POS[1], OBSERVER_POS[2], 
               color='white', s=600, alpha=0.3)
    ax.scatter(OBSERVER_POS[0], OBSERVER_POS[1], OBSERVER_POS[2], 
               color='#FFD700', s=300, marker='*', label='NOSOTROS')
    
    # Vector de Desplazamiento (Línea desde el centro a nosotros)
    neon_glow_line(ax, [0,0,0], OBSERVER_POS, '#FFD700', core_width=1.0)
    
    # Etiquetas de Coordenadas en 3D
    label_text = f"  NOSOTROS\n  X: {OBSERVER_POS[0]}\n  Y: {OBSERVER_POS[1]}\n  Z: {OBSERVER_POS[2]}"
    ax.text(OBSERVER_POS[0], OBSERVER_POS[1], OBSERVER_POS[2], label_text, color='white', fontsize=9)

    # Líneas punteadas de proyección a los ejes (para dar sentido de profundidad)
    ax.plot([OBSERVER_POS[0], OBSERVER_POS[0]], [OBSERVER_POS[1], OBSERVER_POS[1]], [0, OBSERVER_POS[2]], color='gray', linestyle='--', alpha=0.5)
    ax.plot([0, OBSERVER_POS[0]], [OBSERVER_POS[1], OBSERVER_POS[1]], [0, 0], color='gray', linestyle='--', alpha=0.5)
    ax.plot([OBSERVER_POS[0], OBSERVER_POS[0]], [0, OBSERVER_POS[1]], [0, 0], color='gray', linestyle='--', alpha=0.5)

    ax.set_xlim([-1.2, 1.2])
    ax.set_ylim([-1.2, 1.2])
    ax.set_zlim([-1.2, 1.2])
    
    ax.set_axis_off()
    ax.grid(False)
    
    # Ajustar vista para resaltar la altura Z
    ax.view_init(elev=20, azim=135)
    
    title = "EL GPS CÓSMICO: NUESTRA POSICIÓN REAL"
    plt.title(title, color='white', fontsize=18, fontname='Arial', weight='bold', pad=-20)
    
    # --- CORRECCIÓN AQUÍ: Eliminado letter_spacing ---
    coord_str = f"COORDENADAS: X={OBSERVER_POS[0]} | Y={OBSERVER_POS[1]} | Z={OBSERVER_POS[2]}"
    plt.figtext(0.5, 0.05, coord_str, ha="center", color="#FFD700", fontsize=12, weight='bold')
    plt.figtext(0.5, 0.02, "Desplazamiento vertical (Z) explica la distorsión del Vecino 1", ha="center", color="gray", fontsize=8)


    print(f"\n💾 Guardando gráfica de la fórmula...")
    plt.savefig('universal_formula_plot.png', dpi=300)

    print("🚀 RENDER COMPLETADO.")

if __name__ == "__main__":
    main()
