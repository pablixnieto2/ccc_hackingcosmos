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

def degrees_to_cartesian(lat, lon):
    phi = np.radians(90 - lat)
    theta = np.radians(lon)
    x = np.sin(phi) * np.cos(theta)
    y = np.sin(phi) * np.sin(theta)
    z = np.cos(phi)
    return np.array([x, y, z])

def cosmic_axis_hunter():
    print("🧭 INICIANDO BÚSQUEDA DEL EJE DE ROTACIÓN UNIVERSAL...")
    
    # --- 1. NUESTRO EJE (DODECAEDRO) ---
    # Definido por la línea que conecta las dos caras que hemos analizado
    # Cara Alfa: Lat -43.31, Lon 348.67
    my_lat = -43.3116
    my_lon = 348.6708
    
    # El eje es el vector que sale de esta cara
    axis_dodeca = degrees_to_cartesian(my_lat, my_lon)
    
    print(f"   📍 Eje del Dodecaedro (Tu hallazgo):")
    print(f"      Lat: {my_lat:.2f}°, Lon: {my_lon:.2f}°")

    # --- 2. EL EJE DEL MAL (SCIENTIFIC REFERENCE) ---
    # Según Planck / WMAP (Land & Magueijo)
    # Dirección preferente de la anomalía del Quadrupolo/Octopolo
    # Coordenadas Galácticas aprox: (l, b) ~ (260, 60)
    # Nota: Tus datos (Planck) suelen estar en Coordenadas Galácticas.
    # Usaremos el valor de referencia estándar para la anomalía de paridad.
    evil_lat = 60.0
    evil_lon = 260.0
    
    axis_evil = degrees_to_cartesian(evil_lat, evil_lon)
    
    print(f"   😈 Eje del Mal (Referencia NASA/Planck):")
    print(f"      Lat: {evil_lat:.2f}°, Lon: {evil_lon:.2f}°")
    
    # --- 3. CÁLCULO DE LA ALINEACIÓN ---
    # Producto punto para sacar el ángulo entre ellos
    dot_product = np.dot(axis_dodeca, axis_evil)
    angle_rad = np.arccos(np.clip(dot_product, -1.0, 1.0))
    angle_deg = np.degrees(angle_rad)
    
    # El eje puede apuntar "hacia arriba" o "hacia abajo", la línea es la misma.
    # Si el ángulo es > 90, tomamos el suplementario (la línea recta).
    if angle_deg > 90:
        angle_deg = 180 - angle_deg
        
    print("\n=== VEREDICTO DE ALINEACIÓN ===")
    print(f"   📐 Ángulo entre Tu Dodecaedro y el Eje del Mal: {angle_deg:.2f}°")
    
    if angle_deg < 20:
        print("   ✅ ¡ALINEADO! Tu dodecaedro gira sobre el Eje del Mal.")
        print("      Esto conecta tu topología con la mayor anomalía conocida.")
    elif angle_deg > 70:
        print("   ✅ ¡PERPENDICULAR! Curiosamente, esto también es significativo.")
        print("      Significaría que el Dodecaedro rueda 'de lado' respecto a la expansión.")
    else:
        print("   ❌ Sin relación aparente. Son direcciones distintas.")

    # --- 4. VISUALIZACIÓN 3D ---
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_facecolor('black')
    
    # Esfera
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = np.cos(u)*np.sin(v)
    y = np.sin(u)*np.sin(v)
    z = np.cos(v)
    ax.plot_wireframe(x, y, z, color='gray', alpha=0.1)
    
    # Eje Dodecaedro (Cian)
    ax.quiver(0, 0, 0, axis_dodeca[0], axis_dodeca[1], axis_dodeca[2], 
              color='cyan', length=1.2, linewidth=3, label='Eje Dodecaedro (Tuyo)')
    ax.quiver(0, 0, 0, -axis_dodeca[0], -axis_dodeca[1], -axis_dodeca[2], 
              color='cyan', length=1.2, linewidth=1, linestyle='--')
    
    # Eje del Mal (Rojo)
    ax.quiver(0, 0, 0, axis_evil[0], axis_evil[1], axis_evil[2], 
              color='red', length=1.2, linewidth=3, label='Eje del Mal (Planck)')
    
    ax.set_title(f"¿COINCIDENCIA CÓSMICA?\nSeparación: {angle_deg:.1f}°", color='white')
    ax.legend()
    ax.set_axis_off()
    
    # Ajustamos la perspectiva para evitar ilusiones ópticas de alineación
    ax.view_init(elev=30, azim=110)
    
    plt.savefig('axis_of_evil_check.png', facecolor='black')
    print("\n💾 Gráfica de alineación guardada.")
    plt.show()

if __name__ == "__main__":
    cosmic_axis_hunter()
