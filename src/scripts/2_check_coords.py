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

# Configuración basada en las Notas del Capítulo 2
NSIDE = 8  # Resolución baja usada para la búsqueda inicial de polígonos
CANDIDATOS = [368, 336, 399]

print("--- INTERROGATORIO DE SOSPECHOSOS ---")

for idx in CANDIDATOS:
    # Convertir índice de píxel a coordenadas angulares (theta, phi)
    theta, phi = hp.pix2ang(NSIDE, idx)
    
    # Convertir a Latitud/Longitud (Grados)
    # Theta va de 0 (Norte) a 180 (Sur). Latitud = 90 - Theta
    # Phi va de 0 a 360. Longitud = Phi
    
    lat = 90 - np.degrees(theta)
    lon = np.degrees(phi)
    
    # Ajuste visual para lon > 180 (opcional, para estilo -180 a 180)
    if lon > 180:
        lon -= 360
        
    print(f"Candidato {idx}: Lon {lon:6.2f} | Lat {lat:6.2f}")

print("-------------------------------------")
print("CRITERIO:")
print("Si Lat está cerca de 0 (+/- 20), es la Galaxia.")
print("Si Lat es mayor de 30, tenemos algo serio.")
