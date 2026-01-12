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

# CONFIGURACIÓN
INPUT_FILE = 'data/raw/COM_CMB_IQU-sevem_2048_R4.00.fits'
OUTPUT_IMG = 'output/kaiser_stebbins_profile.png'

# UBICACIÓN DE LA CICATRIZ (El centro del vórtice que encontramos)
SCAR_LAT = -41.81
SCAR_LON = 354.38

# GEOMETRÍA DEL CORTE
# Vamos a cortar perpendicularmente a la fractura.
# En la imagen visual, la frontera (rojo/cian) parecía ir en diagonal (aprox 45 grados).
# Así que cortaremos en el ángulo opuesto (135 grados) para cruzarla de lleno.
PROFILE_LENGTH = 10.0 # Grados (5 a cada lado)
PROFILE_ANGLE = 135   # Ángulo del corte (Noreste a Suroeste aprox)
SAMPLES = 200         # Puntos de medición a lo largo del corte

def get_profile_coords(center_lat, center_lon, length_deg, angle_deg, n_samples):
    """Genera las coordenadas (lat, lon) a lo largo de una línea de corte."""
    # Convertir a vectores unitarios
    theta = np.radians(90 - center_lat)
    phi = np.radians(center_lon)
    center_vec = hp.ang2vec(theta, phi)
    
    # Vectores base locales
    z_axis = center_vec
    aux = np.array([0,0,1]) if abs(z_axis[2]) < 0.9 else np.array([0,1,0])
    x_axis = np.cross(aux, z_axis); x_axis /= np.linalg.norm(x_axis)
    y_axis = np.cross(z_axis, x_axis)
    
    # Definir dirección del corte
    cut_angle_rad = np.radians(angle_deg)
    cut_dir = x_axis * np.cos(cut_angle_rad) + y_axis * np.sin(cut_angle_rad)
    
    # Generar puntos
    coords = []
    offsets = np.linspace(-length_deg/2, length_deg/2, n_samples)
    
    for off in offsets:
        off_rad = np.radians(off)
        # Moverse desde el centro en la dirección del corte
        # Aproximación de plano tangente (válida para 10 grados)
        p = center_vec + cut_dir * off_rad
        p /= np.linalg.norm(p) # Proyectar de nuevo a la esfera
        coords.append(p)
        
    return np.array(coords), offsets

def main():
    print(f"--- PERFIL KAISER-STEBBINS (CORTE TRANSVERSAL) ---")
    print(f"Cruzando la cicatriz en Lat {SCAR_LAT}, Lon {SCAR_LON}")
    
    # 1. Cargar Mapas
    print("Cargando datos...")
    maps = hp.read_map(INPUT_FILE, field=[0, 1, 2], verbose=False, memmap=True)
    map_I, map_Q, map_U = maps[0], maps[1], maps[2]
    nside = hp.get_nside(map_I)
    
    # 2. Obtener Coordenadas del Corte
    vecs, x_axis_deg = get_profile_coords(SCAR_LAT, SCAR_LON, PROFILE_LENGTH, PROFILE_ANGLE, SAMPLES)
    
    # 3. Extraer Datos a lo largo de la línea
    # Usamos interpolación bilineal para que sea suave
    temps = hp.get_interp_val(map_I, [v[0] for v in vecs], [v[1] for v in vecs], [v[2] for v in vecs])
    q_vals = hp.get_interp_val(map_Q, [v[0] for v in vecs], [v[1] for v in vecs], [v[2] for v in vecs])
    u_vals = hp.get_interp_val(map_U, [v[0] for v in vecs], [v[1] for v in vecs], [v[2] for v in vecs])
    
    # Calcular Polarización Total y Ángulo
    p_intensity = np.sqrt(q_vals**2 + u_vals**2)
    
    # 4. Graficar
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)
    
    # Gráfica 1: Temperatura (Buscamos el Escalón)
    # Suavizamos un poco para quitar ruido de píxel individual
    window = 5
    temps_smooth = np.convolve(temps, np.ones(window)/window, mode='valid')
    x_smooth = x_axis_deg[window//2 : -window//2 + 1]
    
    ax1.plot(x_axis_deg, temps * 1e6, color='gray', alpha=0.3, label='Raw Signal') # MicroKelvins
    ax1.plot(x_smooth, temps_smooth * 1e6, color='#ffcc00', linewidth=2, label='Smoothed Temp')
    
    ax1.set_ylabel('Temperatura (µK)')
    ax1.set_title(f'¿Efecto Kaiser-Stebbins? Perfil de Temperatura a través de la Cicatriz', fontsize=14)
    ax1.axvline(0, color='red', linestyle='--', alpha=0.5, label='Centro de la Cicatriz')
    ax1.legend()
    ax1.grid(True, alpha=0.2)
    
    # Gráfica 2: Polarización (El Muro de Dominio)
    # Aquí deberíamos ver un cambio de comportamiento
    ax2.plot(x_axis_deg, p_intensity * 1e6, color='cyan', linewidth=2, label='Intensidad Polarización')
    
    ax2.set_ylabel('Polarización (µK)')
    ax2.set_xlabel('Distancia desde el centro (Grados)')
    ax2.set_title('Perfil de Energía de Polarización', fontsize=14)
    ax2.axvline(0, color='red', linestyle='--', alpha=0.5)
    ax2.legend()
    ax2.grid(True, alpha=0.2)
    
    plt.tight_layout()
    plt.savefig(OUTPUT_IMG, dpi=150)
    print(f"Gráfica guardada en {OUTPUT_IMG}")
    
    # Análisis numérico rápido del salto
    left_mean = np.mean(temps[:SAMPLES//2-20])
    right_mean = np.mean(temps[SAMPLES//2+20:])
    step = abs(left_mean - right_mean) * 1e6
    print(f"\nANÁLISIS DE SALTO TÉRMICO:")
    print(f"Promedio Izquierda: {left_mean*1e6:.2f} µK")
    print(f"Promedio Derecha:   {right_mean*1e6:.2f} µK")
    print(f"MAGNITUD DEL ESCALÓN: {step:.2f} µK")

if __name__ == "__main__":
    main()