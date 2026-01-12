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

import pandas as pd
import numpy as np
import healpy as hp
from scipy.spatial import cKDTree
import time

# CONFIGURACIÓN
INPUT_FILE = 'data/processed/vertex_trace_647.csv'
N_SIMULATIONS = 500000 

def calculate_moran_i(values, weights_sparse, W_sum):
    """
    Calcula el Índice de Moran (Autocorrelación Espacial) de forma optimizada.
    Mide cuán 'agrupados' están los valores altos con altos y bajos con bajos.
    """
    # Estandarizar valores (Mean 0, Std 1)
    y = values - np.mean(values)
    y_sq_sum = np.sum(y**2)
    
    # Numerador: Suma de (xi * xj) para vecinos
    # Usamos la matriz de pesos pre-calculada (índices de vecinos)
    # weights_sparse es una lista de pares (i, j)
    
    # Para velocidad en Python puro con bucles, vectorizamos:
    # Pero aquí weights_sparse lo haremos como lista de vecinos
    pass 

def fast_moran(values, neighbor_indices, W):
    """Versión vectorizada ultra-rápida para Monte Carlo"""
    y = values - values.mean()
    z = y / values.std()
    
    # Producto cruzado de vecinos
    # Creamos un vector largo con los valores de los vecinos
    # Reordenamos z según neighbor_indices (flattened)
    
    # Nota: neighbor_indices es una matriz (N, k) con los índices de los k vecinos
    # z[neighbor_indices] crea una matriz (N, k) con los valores de los vecinos
    
    z_neighbors = z[neighbor_indices] 
    sum_neighbors = np.sum(z_neighbors, axis=1) # Suma de valores de vecinos para cada punto
    
    numerator = np.dot(z, sum_neighbors)
    return numerator / W

def main():
    print(f"--- ANÁLISIS ESTADÍSTICO DE SIGNIFICANCIA (SIGMA) ---")
    print(f"Objetivo: {N_SIMULATIONS} simulaciones Monte Carlo (Bootstrap Espacial)")
    
    # 1. Cargar Datos
    try:
        df = pd.read_csv(INPUT_FILE)
    except:
        print("Error: No se encuentra el archivo CSV.")
        return
        
    print(f"Datos cargados: {len(df)} puntos de fractura.")
    
    # Extraer coordenadas y valores
    # Necesitamos convertir HEALPix indices a vectores 3D para buscar vecinos
    nside = 256 # El que usamos antes
    indices = df['center_idx'].values
    correlations = df['corr_IP'].values
    
    # Convertir a x,y,z
    vecs = np.array(hp.pix2vec(nside, indices)).T
    
    # 2. Construir Red de Vecinos (KDTree)
    print("Construyendo red de vecindad espacial...")
    tree = cKDTree(vecs)
    
    # Buscamos los 8 vecinos más cercanos para cada punto
    k = 8 
    dists, neighbor_indices = tree.query(vecs, k=k+1) # k+1 porque el primero es él mismo
    
    # Eliminamos la primera columna (distancia 0, él mismo)
    neighbor_indices = neighbor_indices[:, 1:]
    
    # W es el número total de conexiones (N * k)
    W = len(correlations) * k
    
    # 3. Calcular Moran's I Real (El de tu imagen)
    real_moran = fast_moran(correlations, neighbor_indices, W)
    print(f"\n>>> ÍNDICE DE MORAN REAL (TU ESTRUCTURA): {real_moran:.5f}")
    print("    (Valor > 0 indica agrupación. Valor ~0 es aleatorio)")
    
    # 4. Monte Carlo Loop
    print(f"\nIniciando {N_SIMULATIONS} simulaciones aleatorias...")
    print("Barajando valores para romper la estructura pero mantener la estadística...")
    
    start_time = time.time()
    random_morans = []
    
    # Copia para barajar
    shuffled_corrs = correlations.copy()
    
    for i in range(N_SIMULATIONS):
        np.random.shuffle(shuffled_corrs) # Barajar in-place (muy rápido)
        
        r_moran = fast_moran(shuffled_corrs, neighbor_indices, W)
        random_morans.append(r_moran)
        
        if i % 1000 == 0 and i > 0:
            elapsed = time.time() - start_time
            rate = i / elapsed
            remain = (N_SIMULATIONS - i) / rate
            print(f"   Simulación {i}/{N_SIMULATIONS} | T-Restante: {remain:.1f}s", end='\r')
            
    random_morans = np.array(random_morans)
    
    # 5. Resultados y Cálculo de Sigma
    print(f"\n\n--- RESULTADOS FINALES ---")
    
    # P-Valor: Proporción de simulaciones que superaron al real
    # Si Moran Real es positivo (agrupado), buscamos cuántos aleatorios fueron MAYORES
    # Si Moran Real fuera negativo (tablero ajedrez), buscaríamos menores. 
    # Asumimos positivo por la imagen (bloques).
    
    n_higher = np.sum(random_morans >= real_moran)
    p_value = (n_higher + 1) / (N_SIMULATIONS + 1) # +1 para evitar p=0 estricto (Corrección)
    
    # Cálculo aproximado de Sigma (Z-Score)
    mean_rnd = np.mean(random_morans)
    std_rnd = np.std(random_morans)
    sigma = (real_moran - mean_rnd) / std_rnd
    
    print(f"Mean Random Moran: {mean_rnd:.5f} +/- {std_rnd:.5f}")
    print(f"Tu Moran Real:     {real_moran:.5f}")
    print(f"Diferencia:        {real_moran - mean_rnd:.5f}")
    print("-" * 30)
    print(f"P-VALUE: {p_value:.6f}")
    print(f"SIGMA (Z-SCORE): {sigma:.2f} σ")
    print("-" * 30)
    
    if sigma > 5:
        print("VEREDICTO: DESCUBRIMIENTO CONFIRMADO (GOLD STANDARD > 5σ)")
    elif sigma > 3:
        print("VEREDICTO: EVIDENCIA FUERTE (> 3σ). Publicable.")
    elif sigma > 2:
        print("VEREDICTO: Indicio interesante, pero no concluyente (~2σ).")
    else:
        print("VEREDICTO: Ruido compatible con el azar.")

if __name__ == "__main__":
    main()