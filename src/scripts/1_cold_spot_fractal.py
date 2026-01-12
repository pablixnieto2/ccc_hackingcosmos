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
import pandas as pd
from scipy.stats import entropy, pearsonr
import sys
import os

# --- CONFIGURACIÓN DEL FRANCOTIRADOR ---
# Coordenadas de la Mancha Fría (Cold Spot)
TARGET_L = 209.0
TARGET_B = -57.0
INPUT_FILE = 'data/raw/COM_CMB_IQU-sevem_2048_R4.00.fits' # Asegúrate que este nombre es correcto
OUTPUT_FILE = 'data/processed/cold_spot_profile.csv'

def calculate_hurst(ts):
    ts = np.array(ts)
    if len(ts) < 20: return np.nan
    try:
        lags = range(2, min(len(ts)//2, 50))
        tau = []
        for lag in lags:
            rs_values = []
            for start in range(0, len(ts) - lag, lag):
                chunk = ts[start:start+lag]
                if len(chunk) < 2: continue
                R = np.max(np.cumsum(chunk - np.mean(chunk))) - np.min(np.cumsum(chunk - np.mean(chunk)))
                S = np.std(chunk, ddof=1)
                if S > 0: rs_values.append(R/S)
            if rs_values: tau.append((lag, np.mean(rs_values)))
        if len(tau) < 3: return 0.5
        tau_df = pd.DataFrame(tau, columns=['lag', 'rs'])
        return np.polyfit(np.log10(tau_df['lag']), np.log10(tau_df['rs']), 1)[0]
    except: return np.nan

def main():
    print(f"--- COLD SPOT SNIPER MISSION ---")
    print(f"Target: l={TARGET_L}, b={TARGET_B}")

    # 1. Cargar Mapa (Modo Inteligente)
    print("Cargando mapa SEVEM...")
    try:
        # Intentamos leer HDU 1 (Intensidad) y HDU 2 (Polarización) por separado si es necesario
        # Para SEVEM PR4, a veces están juntos o separados.
        # Probamos lectura estándar primero
        maps = hp.read_map(INPUT_FILE, field=None, hdu=1, verbose=False, memmap=True)
        
        # Si devuelve solo 1 campo, buscamos Q/U en el siguiente HDU
        if len(maps) == 1 or maps.ndim == 1:
            print("Formato Split detectado. Leyendo Q/U de HDU 2...")
            map_I = maps
            maps_pol = hp.read_map(INPUT_FILE, field=None, hdu=2, verbose=False, memmap=True)
            map_Q = maps_pol[0]
            map_U = maps_pol[1]
        else:
            map_I = maps[0]
            map_Q = maps[1]
            map_U = maps[2]
            
        nside = hp.get_nside(map_I)
        print(f"Mapa cargado. NSIDE: {nside}")
        
    except Exception as e:
        print(f"ERROR: {e}")
        sys.exit(1)

    # 2. Definir el vector objetivo
    # Convertir Galácticas (l, b) a Theta, Phi
    rot = hp.Rotator(coord=['G', 'C']) # Galáctico a Celestial (Equatorial) no es necesario para pixels, healpy usa coord interna
    # Pero ojo: Healpy por defecto suele estar en coordenadas Galácticas si el FITS lo dice.
    # Asumimos que el mapa está en coordenadas GALÁCTICAS (estándar Planck).
    theta = np.radians(90.0 - TARGET_B)
    phi = np.radians(TARGET_L)
    
    # Calculamos el vector con los radianes ya convertidos
    center_vec = hp.ang2vec(theta, phi)

    # 3. Biopsia Radial
    # Vamos a expandirnos desde el centro grado a grado
    print("Iniciando escaneo radial...")
    results = []
    
    # Escaneamos desde 0.5 grados hasta 25 grados de radio
    radii = np.arange(0.5, 25.0, 0.5) 
    
    for r_deg in radii:
        r_rad = np.radians(r_deg)
        width_rad = np.radians(0.5) # Anillo fino
        
        # Obtener píxeles del anillo
        inner = r_rad - width_rad/2
        outer = r_rad + width_rad/2
        
        pix_outer = hp.query_disc(nside, center_vec, outer)
        pix_inner = hp.query_disc(nside, center_vec, inner)
        ring_pixels = np.setdiff1d(pix_outer, pix_inner)
        
        if len(ring_pixels) < 100: continue
        
        # Extraer valores
        vals_I = map_I[ring_pixels]
        vals_Q = map_Q[ring_pixels]
        vals_U = map_U[ring_pixels]
        vals_P = np.sqrt(vals_Q**2 + vals_U**2)
        
        # Métricas
        h_I = calculate_hurst(vals_I)
        h_P = calculate_hurst(vals_P)
        corr, _ = pearsonr(vals_I, vals_P)
        
        print(f"Radio {r_deg:5.1f}° | Hurst I: {h_I:.3f} | Corr I-P: {corr:.3f}")
        
        results.append({
            'radio_deg': r_deg,
            'hurst_I': h_I,
            'hurst_P': h_P,
            'corr_IP': corr
        })

    # 4. Guardar
    pd.DataFrame(results).to_csv(OUTPUT_FILE, index=False)
    print(f"Biopsia completada. Datos en {OUTPUT_FILE}")

if __name__ == "__main__":
    main()