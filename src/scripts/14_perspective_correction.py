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

def degrees_to_radians(deg):
    return deg * np.pi / 180

def spherical_to_cartesian(lat, lon, r=1.0):
    # Latitud (-90 a 90), Longitud (0 a 360)
    phi = degrees_to_radians(90 - lat) # Colatitud
    theta = degrees_to_radians(lon)
    x = r * np.sin(phi) * np.cos(theta)
    y = r * np.sin(phi) * np.sin(theta)
    z = r * np.cos(phi)
    return np.array([x, y, z])

def project_vector_to_visual_plane(view_dir, target_vector):
    """
    Proyecta un vector del mundo 3D en el plano 2D de la cámara (perpendicular a la vista)
    para medir su ángulo aparente.
    """
    # Normalizamos la dirección de vista
    k = view_dir / np.linalg.norm(view_dir)
    
    # Restamos la componente paralela a la vista para quedarnos solo con la perpendicular
    projection = target_vector - np.dot(target_vector, k) * k
    return projection

def final_perspective_correction():
    print("🎓 INICIANDO CORRECCIÓN MATEMÁTICA DE PERSPECTIVA (MODO SIBARITA)...")
    
    # --- 1. CONFIGURACIÓN DEL ESCENARIO ---
    
    # A) La Verdad Geométrica (El Universo Ideal)
    # Centro de la Cara Alfa (Donde miramos)
    # Lat/Lon extraídos de tu investigación
    face_center_geo = spherical_to_cartesian(-43.3116, 348.6708, r=1.0)
    
    # B) Nuestra Posición Imperfecta (El Observador Desplazado)
    # Coordenadas que calculó tu GPS
    observer_pos = np.array([0.12, -0.04, 0.28]) 
    
    print(f"   📍 Objetivo: Cara Alfa (Vector Unitario)")
    print(f"   📍 Observador: {observer_pos} (Desplazado del centro)")
    
    # --- 2. CONSTRUCCIÓN DE LA GEOMETRÍA LOCAL ---
    # Creamos un sistema de coordenadas en la superficie de la cara (Plano Tangente)
    # Vector Normal = face_center_geo
    
    # Necesitamos un vector "Arriba" local en la cara para medir ángulos
    # Usamos el Norte aproximado proyectado
    north_pole = np.array([0, 0, 1])
    tangent_u = np.cross(north_pole, face_center_geo) 
    tangent_u /= np.linalg.norm(tangent_u) # Eje X local
    tangent_v = np.cross(face_center_geo, tangent_u) # Eje Y local
    
    # --- 3. SIMULACIÓN DEL GIRO REAL (36 GRADOS) ---
    print("\n   🔄 Simulando giro físico de 36° en la superficie del Dodecaedro...")
    
    # Definimos un punto "P" en el borde del pentágono (radio arbitrario 0.1)
    # Posición Original (0 grados)
    p0_local = tangent_u * 0.1
    p0_world = face_center_geo + p0_local
    
    # Posición Rotada (36 grados) - La predicción teórica
    angle_rad = degrees_to_radians(36.0)
    p36_local = (tangent_u * np.cos(angle_rad) + tangent_v * np.sin(angle_rad)) * 0.1
    p36_world = face_center_geo + p36_local
    
    # --- 4. LO QUE VE EL OJO HUMANO (PROYECCIÓN) ---
    print("   👁️ Calculando distorsión visual desde nuestra posición...")
    
    # Vector de visión: Desde el Observador hasta el Centro de la Cara
    view_vector = face_center_geo - observer_pos
    dist = np.linalg.norm(view_vector)
    print(f"      -> Distancia visual: {dist:.4f} radios")
    
    # Vectores visuales hacia los puntos P0 y P36
    # (Restamos la posición del observador)
    vis_p0 = p0_world - observer_pos
    vis_p36 = p36_world - observer_pos
    
    # Proyectamos estos vectores en el plano visual (lo que vemos en el telescopio)
    # Centramos la vista en el centro de la cara
    proj_p0 = project_vector_to_visual_plane(view_vector, vis_p0 - view_vector)
    proj_p36 = project_vector_to_visual_plane(view_vector, vis_p36 - view_vector)
    
    # --- 5. MEDICIÓN DEL ÁNGULO APARENTE ---
    # Ángulo entre los dos vectores proyectados
    # cos(theta) = dot(a,b) / (|a|*|b|)
    
    dot_prod = np.dot(proj_p0, proj_p36)
    norms = np.linalg.norm(proj_p0) * np.linalg.norm(proj_p36)
    apparent_angle_rad = np.arccos(dot_prod / norms)
    apparent_angle_deg = np.degrees(apparent_angle_rad)
    
    # Corrección de signo (producto cruz) para saber si creció o decreció
    cross_prod = np.cross(proj_p0, proj_p36)
    # Proyectar sobre la dirección de vista para ver el sentido de giro
    direction = np.dot(cross_prod, view_vector)
    if direction < 0: # Depende de la orientación del sistema
         pass 

    print("\n=== VEREDICTO DEL CORRECTOR MATEMÁTICO ===")
    print(f"   📐 Ángulo Real en el Dodecaedro: 36.00°")
    print(f"   🔭 Ángulo Aparente (Observado):  {apparent_angle_deg:.4f}°")
    
    print(f"\n   🔎 Comparación con tu dato observado (48°):")
    diff = abs(apparent_angle_deg - 48.0)
    
    if diff < 2.0:
        print("   ✅ ¡EUREKA! La corrección es perfecta.")
        print("      La deformación visual explica EXACTAMENTE el desplazamiento.")
        print("      No hay error: Es perspectiva.")
    elif diff < 5.0:
        print("   ⚠️ Coincidencia Fuerte.")
        print("      La perspectiva explica la mayor parte del error.")
    else:
        print("   ❌ No coincide. La aberración visual no es suficiente.")

    # --- 6. VISUALIZACIÓN DEL ERROR ---
    plt.style.use('dark_background')
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Dibujar el ángulo teórico (Verde)
    ax.arrow(0, 0, 1, 0, color='green', head_width=0.05, label='Eje 0°')
    ax.arrow(0, 0, np.cos(degrees_to_radians(36)), np.sin(degrees_to_radians(36)), 
             color='lime', head_width=0.05, label='Teórico (36°)')
    
    # Dibujar el ángulo aparente (Rojo)
    ax.arrow(0, 0, np.cos(degrees_to_radians(apparent_angle_deg)), np.sin(degrees_to_radians(apparent_angle_deg)), 
             color='red', head_width=0.05, linestyle='-', width=0.01, label=f'Visual Corregido ({apparent_angle_deg:.1f}°)')
    
    # Dibujar el dato observado (Cian)
    ax.arrow(0, 0, np.cos(degrees_to_radians(48)), np.sin(degrees_to_radians(48)), 
             color='cyan', head_width=0.05, linestyle='--', label='Dato Observado (48°)')
    
    ax.set_xlim(-0.1, 1.2)
    ax.set_ylim(-0.1, 1.2)
    ax.set_aspect('equal')
    ax.legend(loc='upper left')
    plt.title(f"CORRECCIÓN DE PERSPECTIVA\n¿Explica la geometría el error de 12°?", color='white')
    
    print("\n💾 Guardando evidencia matemática...")
    plt.savefig('perspective_correction_proof.png')
    plt.show()

if __name__ == "__main__":
    final_perspective_correction()