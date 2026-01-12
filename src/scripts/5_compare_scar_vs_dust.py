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

import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

# CONFIGURACIÓN
INPUT_FILE = 'data/raw/COM_CMB_IQU-sevem_2048_R4.00.fits'
OUTPUT_IMG = 'output/scar_vs_dust_comparison.png'

# COORDENADAS
SCAR_LAT = -41.81
SCAR_LON = 354.38

# COORDENADAS DEL "CONTROL SUCIO" (LA GALAXIA)
# Cogemos el punto más cercano en el ecuador galáctico
DUST_LAT = 0.00
DUST_LON = 354.38

def get_patch_pixels(map_data, lat, lon, radius_deg=5.0):
    """Extrae los valores de los píxeles en un disco de forma robusta."""
    # CORRECCIÓN: Conversión manual a Radianes (Theta, Phi)
    # Theta (Colatitude) va de 0 (Norte) a 180 (Sur).
    # Phi (Longitud) va de 0 a 360 (convertido a radianes).
    
    theta = np.radians(90.0 - lat)
    phi = np.radians(lon)
    
    # Usamos ang2vec estándar (espera radianes)
    vec = hp.ang2vec(theta, phi)
    
    nside = hp.get_nside(map_data)
    # Extraemos píxeles del disco
    pix_indices = hp.query_disc(nside, vec, np.radians(radius_deg))
    return map_data[pix_indices]

def main():
    print(f"--- COMPARATIVA FORENSE: CICATRIZ vs POLVO GALÁCTICO ---")
    
    # 1. Cargar Mapa (Solo Intensidad esta vez)
    print("Cargando mapa...")
    try:
        # Intentamos cargar campo 0 (Intensidad)
        maps = hp.read_map(INPUT_FILE, field=0, verbose=False, memmap=True)
    except:
        # Fallback por si acaso
        maps = hp.read_map(INPUT_FILE, field=None, hdu=1, verbose=False, memmap=True)
    
    # Asegurar que es 1D array
    if maps.ndim > 1: map_I = maps[0]
    else: map_I = maps

    # 2. Extraer Muestras
    print(f"Extrayendo biopsia de Cicatriz (Lat {SCAR_LAT})...")
    scar_vals = get_patch_pixels(map_I, SCAR_LAT, SCAR_LON)
    
    print(f"Extrayendo biopsia de Polvo (Lat {DUST_LAT})...")
    dust_vals = get_patch_pixels(map_I, DUST_LAT, DUST_LON)
    
    # Convertir a MicroKelvins
    scar_vals *= 1e6
    dust_vals *= 1e6

    # 3. Visualización
    print("Generando gráficos comparativos...")
    fig = plt.figure(figsize=(15, 10), facecolor='black')
    
    # --- VISUAL: LA CICATRIZ ---
    # Para gnomview, rot = (Lon, Lat) en grados funciona bien
    rot_scar = [SCAR_LON, SCAR_LAT, 0]
    hp.gnomview(map_I, rot=rot_scar, xsize=300, ysize=300, reso=3.0, 
                sub=(2, 2, 1), title=f"MUESTRA A: LA CICATRIZ\n(Lat {SCAR_LAT})",
                unit='K_cmb', cmap='inferno', hold=True, notext=False)
    
    # --- VISUAL: EL POLVO ---
    rot_dust = [DUST_LON, DUST_LAT, 0]
    hp.gnomview(map_I, rot=rot_dust, xsize=300, ysize=300, reso=3.0, 
                sub=(2, 2, 2), title=f"MUESTRA B: POLVO GALÁCTICO\n(Lat {DUST_LAT})",
                unit='K_cmb', cmap='inferno', hold=True, notext=False)

    # --- ESTADÍSTICA: HISTOGRAMAS ---
    ax_hist = fig.add_subplot(2, 1, 2, facecolor='#111111')
    
    # Plot Histograma Cicatriz
    ax_hist.hist(scar_vals, bins=100, color='#ffcc00', alpha=0.6, density=True, label='Cicatriz (Muestra A)')
    
    # Plot Histograma Polvo
    # El polvo tiene picos muy altos, limitamos el rango visual para comparar la forma de la base
    ax_hist.hist(dust_vals, bins=100, color='cyan', alpha=0.4, density=True, label='Polvo Galáctico (Muestra B)')
    
    ax_hist.set_title("HUELLA DACTILAR ESPECTRAL (Histograma de Intensidad)", color='white', fontsize=14)
    ax_hist.set_xlabel("Temperatura (µK)", color='gray')
    ax_hist.set_ylabel("Frecuencia Relativa", color='gray')
    ax_hist.legend()
    ax_hist.grid(True, alpha=0.2)
    ax_hist.tick_params(colors='gray')
    
    # Estadísticas
    scar_skew = np.mean((scar_vals - np.mean(scar_vals))**3) / np.std(scar_vals)**3
    dust_skew = np.mean((dust_vals - np.mean(dust_vals))**3) / np.std(dust_vals)**3
    
    stats_text = (
        f"--- ANÁLISIS ESTADÍSTICO ---\n"
        f"SKEWNESS (Asimetría):\n"
        f"  Cicatriz: {scar_skew:.2f}\n"
        f"  Polvo:    {dust_skew:.2f}\n"
        f"\n"
        f"STD DEV (Energía):\n"
        f"  Cicatriz: {np.std(scar_vals):.2f} uK\n"
        f"  Polvo:    {np.std(dust_vals):.2f} uK"
    )
    plt.figtext(0.02, 0.02, stats_text, color='white', fontsize=12, bbox=dict(facecolor='black', alpha=0.7))

    print(f"Guardando comparativa en {OUTPUT_IMG}...")
    plt.savefig(OUTPUT_IMG, dpi=150, facecolor='black')
    print("¡Análisis completado!")

if __name__ == "__main__":
    main()