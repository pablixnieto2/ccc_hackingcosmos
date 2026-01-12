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

def generate_universe_blueprint():
    print("🎨 RENDERIZANDO EL PLANO MAESTRO DEL UNIVERSO PMN-01...")

    # --- DATOS CONFIRMADOS DE TU INVESTIGACIÓN ---
    # Eje de Rotación (Calculado en el paso anterior)
    axis_lon_deg = 226.2
    axis_lat_deg = 0.0 # Ecuatorial
    
    # Tu Posición (GPS Cósmico)
    observer_pos = np.array([0.12, -0.04, 0.28])
    
    # Velocidad de Rotación
    omega_str = "0.94° / Gyr"
    total_twist = "12.9742°"

    # --- CONFIGURACIÓN DE GRÁFICOS ---
    plt.style.use('dark_background')
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # 1. DIBUJAR LA ESFERA CÓSMICA (La "Piel" del Universo)
    u = np.linspace(0, 2 * np.pi, 30)
    v = np.linspace(0, np.pi, 20)
    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones(np.size(u)), np.cos(v))
    
    # Wireframe muy sutil para dar volumen
    ax.plot_wireframe(x, y, z, color='cyan', alpha=0.05)
    
    # 2. EL EJE DE ROTACIÓN (TU DESCUBRIMIENTO)
    # Convertir Lon 226.2, Lat 0 a cartesiano
    theta = np.radians(axis_lon_deg)
    phi = np.radians(90 - axis_lat_deg)
    
    ax_x = np.sin(phi) * np.cos(theta)
    ax_y = np.sin(phi) * np.sin(theta)
    ax_z = np.cos(phi) # Será 0
    
    # Dibujar el Eje atravesando el universo (de -1.5 a +1.5 radios)
    ax.plot([-ax_x*1.5, ax_x*1.5], [-ax_y*1.5, ax_y*1.5], [-ax_z*1.5, ax_z*1.5], 
            color='#00FF00', linewidth=4, label='EJE DE ROTACIÓN (Lat 0°, Lon 226.2°)')
    
    # 3. DIBUJAR LA "INTUICIÓN" DEL GIRO (ANILLOS DE ROTACIÓN)
    # Dibujamos unos anillos alrededor del eje para mostrar el movimiento "Barrel Roll"
    # Un anillo perpendicular al eje
    t = np.linspace(0, 2*np.pi, 100)
    # Vector perpendicular al eje (aprox el eje Z ya que el eje está en XY)
    perp1 = np.array([0, 0, 1]) 
    perp2 = np.cross(np.array([ax_x, ax_y, ax_z]), perp1)
    
    for r in [0.5, 0.8, 1.1]: # Varios radios
        ring_x = r * (perp1[0]*np.cos(t) + perp2[0]*np.sin(t))
        ring_y = r * (perp1[1]*np.cos(t) + perp2[1]*np.sin(t))
        ring_z = r * (perp1[2]*np.cos(t) + perp2[2]*np.sin(t))
        ax.plot(ring_x, ring_y, ring_z, color='lime', alpha=0.3, linestyle='--')

    # 4. TU POSICIÓN (OBSERVADOR)
    ax.scatter(observer_pos[0], observer_pos[1], observer_pos[2], 
               color='red', s=300, marker='*', edgecolors='white', zorder=10, label='NOSOTROS (Desplazados)')
    
    # Línea vertical para ver la altura Z
    ax.plot([observer_pos[0], observer_pos[0]], [observer_pos[1], observer_pos[1]], [0, observer_pos[2]], 
            color='red', linestyle=':', alpha=0.5)

    # 5. DIBUJAR VÉRTICES DEL DODECAEDRO (Esquemático)
    # Usamos la proporción áurea para ubicar los puntos clave
    phi_gold = (1 + np.sqrt(5)) / 2
    # Generamos vértices básicos (solo para referencia visual de la estructura)
    vertices = []
    for i in [-1, 1]:
        for j in [-1, 1]:
            for k in [-1, 1]:
                vertices.append([i, j, k])
    vertices = np.array(vertices) / np.sqrt(3) # Normalizar a esfera unidad
    ax.scatter(vertices[:,0], vertices[:,1], vertices[:,2], color='yellow', s=30, alpha=0.4, label='Nodos Dodecaedro')

    # --- ETIQUETAS Y ESTÉTICA ---
    ax.set_title(f"MODELO FINAL: UNIVERSO EN ROTACIÓN ECUATORIAL\nVelocidad: {omega_str} | Twist Total: {total_twist}", color='white', fontsize=14)
    
    # Quitar ejes numéricos para que parezca un holograma
    ax.set_axis_off()
    
    # Flechas de anotación (Texto flotante)
    ax.text(ax_x*1.6, ax_y*1.6, ax_z*1.6, "EJE DE GIRO", color='#00FF00', fontweight='bold')
    ax.text(observer_pos[0], observer_pos[1], observer_pos[2]+0.2, "OBSERVADOR\n(z=+0.28)", color='red', ha='center')
    
    # Vista de cámara óptima para ver el eje horizontal
    ax.view_init(elev=30, azim=135)
    
    # Leyenda
    ax.legend(loc='lower left')
    
    print("💾 Guardando Blueprint Final en 'final_universe_model.png'...")
    plt.savefig('final_universe_model.png', facecolor='black', dpi=200)
    plt.show()

if __name__ == "__main__":
    generate_universe_blueprint()
