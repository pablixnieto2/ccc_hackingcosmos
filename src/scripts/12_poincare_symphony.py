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
import healpy as hp
from scipy.io import wavfile
import os

# --- RUTA A TU ARCHIVO REAL ---
# Asegúrate de que esta ruta sea correcta en tu Mac
FITS_FILE = "data/raw/COM_CMB_IQU-sevem_2048_R4.00.fits"

def real_data_symphony():
    print(f"📡 CONECTANDO CON EL ARCHIVO: {FITS_FILE}")
    
    if not os.path.exists(FITS_FILE):
        print("❌ ERROR: No encuentro el archivo. Revisa la ruta.")
        return

    # 1. CARGAR EL MAPA REAL (PLANCK SEVEM)
    print("   ⏳ Leyendo mapa FITS (esto puede tardar unos segundos)...")
    # Campo 0 = Temperatura (I_STOKES)
    map_I = hp.read_map(FITS_FILE, field=0, verbose=False)
    
    # 2. EXTRAER EL ESPECTRO DE POTENCIA (ANAFAST)
    print("   🧮 Calculando Espectro de Potencia Real (C_l) con AnaFast...")
    # Calculamos hasta l=2000 (suficiente para audio)
    cl_real = hp.anafast(map_I, lmax=2000)
    l_axis = np.arange(len(cl_real))
    
    # Normalizamos para visualización y audio (evitar valores minúsculos tipo 1e-12)
    # Multiplicamos por l(l+1)/2pi que es el estándar cosmológico D_l
    dl_real = l_axis * (l_axis + 1) * cl_real / (2 * np.pi)
    
    # Limpieza básica: quitar monopolos/dipolos (l=0, 1) que suelen ser ruido enorme
    dl_real[0:2] = 0 

    # 3. ANÁLISIS FORENSE DE BAJOS (TU PREGUNTA SOBRE QUADRUPOLO)
    print("\n🔍 --- DATOS REALES EXTRAÍDOS ---")
    print(f"   Quadrupolo (l=2) Real: {dl_real[2]:.6f}")
    print(f"   Octopolo   (l=3) Real: {dl_real[3]:.6f}")
    
    # Comparativa rápida con un modelo estándar simple (aprox) para ver la anomalía
    # (Valor promedio esperado aprox ~1000-1200 en estas unidades)
    print(f"   ¿Es bajo? (Referencia aprox ~1000):")
    if dl_real[2] < 500:
        print("   -> SÍ. El Quadrupolo es ANÓMALAMENTE BAJO (Confirmado).")
    else:
        print("   -> No parece tan bajo en este mapa.")

    # 4. GENERAR AUDIO CON DATOS REALES
    print("\n🎹 Generando 'La Verdadera Música de las Esferas'...")
    duration = 10 # segundos
    sample_rate = 44100
    t = np.linspace(0, duration, int(sample_rate * duration), endpoint=False)
    
    audio_signal = np.zeros_like(t)
    
    # Síntesis aditiva usando el espectro REAL
    # Usamos l de 2 a 1500
    for ell in range(2, 1500, 5): # Saltos de 5 para velocidad
        amp = dl_real[ell]
        if amp > 0:
            # Mapeo: l=2 -> 50Hz (Graves profundos)
            freq = 50 + (ell * 0.5) 
            # Normalizar amplitud localmente para que se oiga
            volume = np.sqrt(amp) 
            audio_signal += volume * np.sin(2 * np.pi * freq * t)

    # Normalizar master
    audio_signal = audio_signal / np.max(np.abs(audio_signal))
    wavfile.write("REAL_PLANCK_SYMPHONY.wav", sample_rate, audio_signal.astype(np.float32))
    print("   ✅ Audio guardado: REAL_PLANCK_SYMPHONY.wav")

    # 5. GRAFICAR LA REALIDAD
    plt.style.use('dark_background')
    plt.figure(figsize=(12, 6))
    
    plt.plot(l_axis, dl_real, color='cyan', alpha=0.8, label='Datos Reales (SEVEM)')
    plt.fill_between(l_axis, 0, dl_real, color='cyan', alpha=0.2)
    
    # Marcar los puntos críticos
    plt.scatter(2, dl_real[2], color='red', s=100, zorder=10, label='Quadrupolo Real')
    plt.scatter(3, dl_real[3], color='orange', s=100, zorder=10, label='Octopolo Real')
    
    plt.title(f"ESPECTRO DE POTENCIA REAL (Extraído de {os.path.basename(FITS_FILE)})", fontsize=14)
    plt.xlabel("Multipolo (l)")
    plt.ylabel("Potencia D_l")
    plt.legend()
    plt.xlim(0, 2000)
    plt.grid(True, alpha=0.2)
    
    plt.savefig('REAL_POWER_SPECTRUM.png', dpi=150)
    print("   ✅ Gráfica guardada: REAL_POWER_SPECTRUM.png")
    plt.show()

if __name__ == "__main__":
    real_data_symphony()