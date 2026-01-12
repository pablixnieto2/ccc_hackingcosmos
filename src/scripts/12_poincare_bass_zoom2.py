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
import os

# --- CONFIGURACIÓN ---
# Ruta a TU archivo real (ajusta si es necesario)
FITS_FILE = "data/raw/COM_CMB_IQU-sevem_2048_R4.00.fits"

def real_forensic_analysis():
    print(f"🕵️‍♂️ INICIANDO FORENSE CON DATOS REALES: {FITS_FILE}")
    
    if not os.path.exists(FITS_FILE):
        print(f"❌ ERROR CRÍTICO: No encuentro el archivo en {FITS_FILE}")
        return

    # 1. CARGA Y CÁLCULO REAL
    print("   ⏳ Cargando mapa y calculando espectro (esto tarda unos segundos)...")
    # Leemos el campo de Temperatura (I)
    map_I = hp.read_map(FITS_FILE, field=0, verbose=False)
    
    # Calculamos el espectro C_l hasta l=30 (solo bajos)
    cl_real = hp.anafast(map_I, lmax=35)
    l_axis = np.arange(len(cl_real))
    
    # Convertimos a D_l (Unidades estándar muK^2)
    # D_l = l * (l + 1) * C_l / (2 * pi)
    # * Factor de escala del Planck (habitualmente 1e12 para muK^2 si viene en K)
    # Asumimos unidades estándar del FITS (K_cmb). Ajustamos escala visual.
    dl_real = l_axis * (l_axis + 1) * cl_real / (2 * np.pi) * 1e12 
    
    # Limpiamos l=0 y l=1 (no son física cosmológica)
    dl_real[0] = 0
    dl_real[1] = 0

    # 2. MODELO ESTÁNDAR (REFERENCIA INFINITA)
    # Curva teórica aproximada Lambda-CDM para comparar
    standard_model = 1000 * (l_axis / (l_axis + 1.5)) + 850
    # Ajuste de altura para que coincida con los picos altos de tu data (autocalibración)
    scale_factor = np.mean(dl_real[10:30]) / np.mean(standard_model[10:30])
    standard_model = standard_model * scale_factor

    # 3. VISUALIZACIÓN FORENSE
    print("   📊 Generando la escena del crimen...")
    plt.style.use('dark_background')
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # --- ZONA ROJA (PROHIBIDA) ---
    ax.axvspan(1.5, 3.5, color='red', alpha=0.15, label='Zona de Supresión (Dodecaedro)')

    # Modelo Teórico (Verde)
    ax.plot(l_axis[2:31], standard_model[2:31], color='#00FF00', linestyle='--', alpha=0.5, label='Teoría Universo Infinito')
    ax.fill_between(l_axis[2:31], 0, standard_model[2:31], color='#00FF00', alpha=0.05)
    
    # TUS DATOS REALES (Azul)
    ax.plot(l_axis[2:31], dl_real[2:31], color='cyan', alpha=0.4, linewidth=2, label='TUS DATOS (Planck Real)')
    ax.scatter(l_axis[2:31], dl_real[2:31], color='cyan', s=40, alpha=0.8)

    # --- ANÁLISIS AUTOMÁTICO DE SOSPECHOSOS ---
    
    # CUADRUPOLO (l=2)
    q_val = dl_real[2]
    q_theo = standard_model[2]
    # Calculamos la pérdida real
    loss_q = int(100 - (q_val/q_theo)*100)
    
    color_q = 'red' if loss_q > 50 else 'yellow' # Rojo si falta mucho, amarillo si es normal
    ax.scatter(2, q_val, color=color_q, s=300, edgecolors='white', zorder=10, label='Quadrupolo Real (l=2)')
    
    ax.annotate(f'QUADRUPOLO REAL\nFalta {loss_q}% de Energía', 
                xy=(2, q_val), xytext=(2.5, q_val + (q_theo*0.2)),
                arrowprops=dict(facecolor=color_q, shrink=0.05),
                color=color_q, fontweight='bold', bbox=dict(boxstyle="round", fc="black", ec=color_q))

    # OCTOPOLO (l=3)
    o_val = dl_real[3]
    o_theo = standard_model[3]
    loss_o = int(100 - (o_val/o_theo)*100)
    
    color_o = 'orange' if loss_o > 40 else 'yellow'
    ax.scatter(3, o_val, color=color_o, s=250, edgecolors='white', zorder=10, label='Octopolo Real (l=3)')
    
    ax.annotate(f'OCTOPOLO REAL\nFalta {loss_o}% de Energía', 
                xy=(3, o_val), xytext=(4.5, o_val + (o_theo*0.1)),
                arrowprops=dict(facecolor=color_o, shrink=0.05),
                color=color_o, fontweight='bold', bbox=dict(boxstyle="round", fc="black", ec=color_o))

    # Líneas de referencia (La caída)
    ax.vlines(2, q_val, q_theo, colors=color_q, linestyles='solid', linewidth=2, alpha=0.5)
    ax.vlines(3, o_val, o_theo, colors=color_o, linestyles='solid', linewidth=2, alpha=0.5)

    # Decoración Final
    plt.title(f"FORENSE DE DATOS REALES: {os.path.basename(FITS_FILE)}", fontsize=16, color='white', pad=20)
    plt.xlabel("Multipolo (l)", fontsize=12)
    plt.ylabel("Potencia Real ($\mu K^2$)", fontsize=12)
    plt.xlim(1.5, 15)
    plt.grid(True, alpha=0.15)
    plt.legend(loc='upper right')
    
    output_img = 'real_data_bass_forensics.png'
    print(f"💾 Guardando la verdad en {output_img}...")
    plt.savefig(output_img, dpi=150, facecolor='black')
    print("✅ GRÁFICA GENERADA CON ÉXITO.")
    
    # Imprimir veredicto en texto
    print("\n--- VEREDICTO DE TUS DATOS ---")
    print(f"Quadrupolo (l=2): {q_val:.2f} (Teórico: {q_theo:.2f}) -> PÉRDIDA: {loss_q}%")
    print(f"Octopolo   (l=3): {o_val:.2f} (Teórico: {o_theo:.2f}) -> PÉRDIDA: {loss_o}%")
    
    if loss_q > 70 and loss_o > 50:
        print("\n>>> CONCLUSIÓN: ¡LA ANOMALÍA ES REAL EN TU ARCHIVO!")
        print(">>> Los bajos están suprimidos. El Dodecaedro es compatible.")
    else:
        print("\n>>> CONCLUSIÓN: Datos ambiguos. La supresión no es tan fuerte como se esperaba.")

    plt.show()

if __name__ == "__main__":
    real_forensic_analysis()