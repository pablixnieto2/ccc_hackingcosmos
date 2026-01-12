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

def dark_energy_buster():
    print("💥 INICIANDO PROTOCOLO: DARK ENERGY BUSTER...")
    print("   Objetivo: Comprobar si la rotación del Dodecaedro explica la expansión.")

    # --- 1. DATOS DE TU UNIVERSO (PMN-01) ---
    # Velocidad angular calculada basada en la torsión precisa
    twist_degrees = 12.9742
    age_universe_gyr = 13.8
    omega_deg_per_gyr = twist_degrees / age_universe_gyr
    
    # Conversión a unidades físicas (Sistema Internacional - SI)
    # 1 Gyr (Giga-año) en segundos
    gyr_in_seconds = 1e9 * 365.25 * 24 * 3600
    # Grados a Radianes
    omega_rad_s = np.radians(omega_deg_per_gyr) / gyr_in_seconds
    
    print(f"   🌪️ Velocidad de Rotación (omega): {omega_rad_s:.2e} rad/s")

    # --- 2. FÍSICA TEÓRICA ---
    # La aceleración centrífuga es a = omega^2 * r
    # Necesitamos un radio de referencia. Usaremos el Radio de Hubble (el horizonte visible).
    # Radio de Hubble ~ 14.4 Giga-parsecs ~ 4.4e26 metros
    r_hubble = 4.4e26 
    
    # Aceleración centrífuga generada por TU rotación
    a_centrifugal = (omega_rad_s ** 2) * r_hubble
    
    print(f"   🏎️ Aceleración Centrífuga generada: {a_centrifugal:.2e} m/s^2")

    # --- 3. COMPARACIÓN CON LA "ENERGÍA OSCURA" (Lambda) ---
    # La aceleración cósmica medida por la NASA (H0^2 * Omega_Lambda * r)
    # H0 ~ 70 km/s/Mpc ~ 2.2e-18 / s
    H0 = 2.2e-18
    Omega_Lambda = 0.69 # El 69% del universo es Energía Oscura según ellos
    
    # Aceleración debida a la Energía Oscura (según modelo estándar)
    a_dark_energy = (H0 ** 2) * Omega_Lambda * r_hubble
    
    print(f"   👻 Aceleración por Energía Oscura (NASA): {a_dark_energy:.2e} m/s^2")
    
    # --- 4. VEREDICTO ---
    ratio = a_centrifugal / a_dark_energy
    percent_explained = ratio * 100
    
    print("\n=== VEREDICTO DE UNIFICACIÓN ===")
    print(f"   Tu rotación explica el {percent_explained:.4f}% de la Energía Oscura.")
    
    # Gráfico de Fuerzas
    plt.style.use('dark_background')
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Barras
    labels = ['Energía Oscura (NASA)', 'Fuerza Centrífuga (TUYA)']
    values = [a_dark_energy, a_centrifugal]
    colors = ['purple', 'lime']
    
    bars = ax.bar(labels, values, color=colors, alpha=0.7)
    
    # Escala Logarítmica si son muy diferentes
    # ax.set_yscale('log') 
    
    # Etiquetas de valor
    for bar in bars:
        yval = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, yval, f'{yval:.2e}', ha='center', va='bottom', color='white')

    plt.title(f"¿ES LA ROTACIÓN LA CAUSA DE LA EXPANSIÓN?\nExplicación: {percent_explained:.2f}%", fontsize=14)
    plt.ylabel("Aceleración Cósmica (m/s^2)")
    
    output_path = 'FINAL_PMN/src/18_DARK_ENERGY/images/dark_energy_vs_rotation.png'
    print(f"💾 Guardando resultado en '{output_path}'...")
    plt.savefig(output_path)
    # plt.show() # Commented out for headless execution
    
    if percent_explained > 1.0:
        print("\n✅ RESULTADO POSITIVO: Tu rotación contribuye físicamente a la expansión.")
        print("   No necesitas tanta 'energía mágica' como dicen.")
    else:
        print("\n⚠️ RESULTADO BAJO: La rotación es real, pero demasiado lenta para explicar toda la expansión.")
        print("   Quizás la Energía Oscura es la tensión de la pared del Dodecaedro (Efecto Casimir).")

if __name__ == "__main__":
    dark_energy_buster()
