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

import os
import sys
import gc
import glob
import numpy as np
import healpy as hp
import pandas as pd
from scipy.stats import entropy, pearsonr
from tqdm import tqdm
import multiprocessing as mp
from functools import partial
from astropy.io import fits

# Configuration
INPUT_PATTERN = 'data/raw/*.fits' 
OUTPUT_FILE = 'data/processed/fractal_metrics.csv'

def calculate_hurst(ts):
    """Calculates the Hurst Exponent using Rescaled Range (R/S) analysis."""
    ts = np.array(ts)
    if len(ts) < 20: return np.nan
    try:
        lags = range(2, min(len(ts)//2, 100))
        tau = []
        for lag in lags:
            rs_values = []
            for start in range(0, len(ts) - lag, lag):
                chunk = ts[start:start+lag]
                if len(chunk) < 2: continue
                mean = np.mean(chunk)
                R = np.max(np.cumsum(chunk - mean)) - np.min(np.cumsum(chunk - mean))
                S = np.std(chunk, ddof=1)
                if S > 0: rs_values.append(R/S)
            if rs_values: tau.append((lag, np.mean(rs_values)))
        
        if len(tau) < 3: return 0.5
        tau_df = pd.DataFrame(tau, columns=['lag', 'rs'])
        return np.polyfit(np.log10(tau_df['lag']), np.log10(tau_df['rs']), 1)[0]
    except:
        return np.nan

def calculate_entropy(ts):
    """Calculates Shannon Entropy."""
    if len(ts) == 0: return np.nan
    hist, _ = np.histogram(ts, bins='fd', density=True)
    return entropy(hist)

def get_ring_pixels(nside, center_vec, radius_rad, width_rad=0.01):
    """Returns pixel indices for a ring."""
    inner_rad = radius_rad - width_rad/2
    outer_rad = radius_rad + width_rad/2
    pixels = hp.query_disc(nside, center_vec, outer_rad)
    if inner_rad > 0:
        inner_pixels = hp.query_disc(nside, center_vec, inner_rad)
        pixels = np.setdiff1d(pixels, inner_pixels)
    return pixels

# Global storage for workers
shared_I = None
shared_Q = None
shared_U = None

def init_worker(I, Q, U):
    global shared_I, shared_Q, shared_U
    shared_I = I
    shared_Q = Q
    shared_U = U

def process_ring_memory_optimized(args):
    """
    args: (ring_id, theta, phi, radius_deg, reduced_indices)
    Uses pre-computed reduced indices to access shared memory arrays.
    """
    ring_id, theta, phi, radius_deg, idxs = args
    global shared_I, shared_Q, shared_U
    
    try:
        if len(idxs) < 50: return None
        
        vals_I = shared_I[idxs]
        vals_Q = shared_Q[idxs]
        vals_U = shared_U[idxs]
        
        vals_P = np.sqrt(vals_Q**2 + vals_U**2)
        
        # Check for zero variance
        if np.std(vals_I) == 0 or np.std(vals_P) == 0:
            return None

        hurst_I = calculate_hurst(vals_I)
        entropy_I = calculate_entropy(vals_I)
        hurst_P = calculate_hurst(vals_P)
        corr_IP, _ = pearsonr(vals_I, vals_P)
            
        return {
            'id_anillo': ring_id,
            'theta': theta,
            'phi': phi,
            'radio': radius_deg,
            'hurst_I': hurst_I,
            'entropy_I': entropy_I,
            'hurst_P': hurst_P,
            'corr_IP': corr_IP
        }
    except Exception:
        return None

def get_nside_from_header(filepath):
    try:
        with fits.open(filepath) as hdul:
            # Check HDU 1 first
            if len(hdul) > 1:
                header = hdul[1].header
                if 'NSIDE' in header:
                    return header['NSIDE']
                elif 'NAXIS2' in header:
                    return hp.npix2nside(header['NAXIS2'])
    except Exception as e:
        print(f"Header read failed: {e}")
    
    # Fallback: read field 0
    print("Reading field 0 to determine NSIDE...")
    m = hp.read_map(filepath, field=0, verbose=False)
    nside = hp.get_nside(m)
    del m
    gc.collect()
    return nside

def resolve_input_file(pattern):
    files = glob.glob(pattern)
    files = [f for f in files if f.lower().endswith('.fits')]
    if not files:
        print(f"Error: No FITS files found matching '{pattern}'")
        sys.exit(1)
    
    # Priority: SEVEM > SMICA > Others
    sevem_files = [f for f in files if 'sevem' in f.lower()]
    smica_files = [f for f in files if 'smica' in f.lower()]
    
    if sevem_files:
        selected = sevem_files[0]
        print(f"Auto-selected SEVEM file: {selected}")
    elif smica_files:
        selected = smica_files[0]
        print(f"Auto-selected SMICA file: {selected}")
    else:
        selected = files[0]
        print(f"Using first available FITS: {selected}")
            
    return selected

def main():
    print("Starting CCC Fractal Analysis (PR4/SEVEM Ready)...")
    
    INPUT_FILE = resolve_input_file(INPUT_PATTERN)
    print(f"Target File: {INPUT_FILE}")

    # 1. Determine NSIDE
    print("Determining Map Properties...")
    nside = get_nside_from_header(INPUT_FILE)
    npix = hp.nside2npix(nside)
    print(f"NSIDE: {nside}, Total Pixels: {npix}")

    # 2. Scan Strategy & Index Calculation
    print("Phase 1: Calculating Interest Indices (Geometric Scan)...")
    nside_scan = 8
    npix_scan = hp.nside2npix(nside_scan)
    centers_theta, centers_phi = hp.pix2ang(nside_scan, np.arange(npix_scan))
    radii_to_check = [2.0, 4.0, 6.0, 10.0]
    
    ring_definitions = [] 
    all_needed_pixels = [] 
    ring_pixel_map = [] 
    
    ring_id = 0
    width_rad = np.radians(1.0)
    
    for r_deg in radii_to_check:
        r_rad = np.radians(r_deg)
        for i in range(npix_scan):
            vec = hp.ang2vec(centers_theta[i], centers_phi[i])
            pixels = get_ring_pixels(nside, vec, r_rad, width_rad=width_rad)
            
            if len(pixels) > 50:
                ring_definitions.append((ring_id, centers_theta[i], centers_phi[i], r_deg))
                ring_pixel_map.append(pixels)
                all_needed_pixels.append(pixels)
                ring_id += 1
    
    print(f"Identified {len(ring_definitions)} candidate rings.")
    
    # 3. Consolidate Unique Pixels
    print("Phase 2: Consolidating Unique Pixels...")
    flat_pixels = np.concatenate(all_needed_pixels)
    unique_pixels = np.unique(flat_pixels)
    num_unique = len(unique_pixels)
    print(f"Unique pixels to load: {num_unique} ({num_unique/npix*100:.2f}% of map)")
    
    del all_needed_pixels, flat_pixels
    gc.collect()
    
    # 4. Smart-Split Loading Strategy
    print("Phase 3: Smart-Split Data Loading (Low RAM)...")
    extracted_data = {}
    
    # Check if we need Split HDU logic (PR4) or Standard (PR3/SMICA)
    try:
        # Attempt standard single-HDU load first (Standard Healpix)
        print("  -> Checking HDU 1 structure...")
        test_maps = hp.read_map(INPUT_FILE, field=None, hdu=1, verbose=False, memmap=True)
        
        is_split_hdu = False
        if isinstance(test_maps, np.ndarray):
            # Only 1 map found in HDU 1? Likely Split structure.
            if test_maps.ndim == 1:
                is_split_hdu = True
                print("  -> Detected Split HDU Structure (PR4 Style).")
            elif test_maps.shape[0] < 3:
                is_split_hdu = True
                print("  -> Detected Partial HDU Structure. Switching to Split Logic.")
        elif len(test_maps) < 3:
            is_split_hdu = True
            print("  -> Detected < 3 maps in HDU 1. Switching to Split Logic.")
            
        del test_maps
        gc.collect()
        
    except Exception as e:
        print(f"  -> HDU 1 check warning: {e}. Assuming Split Logic.")
        is_split_hdu = True

    try:
        if not is_split_hdu:
            # STANDARD LOAD (I, Q, U in HDU 1)
            print("  -> Mode: Standard Combined HDU (SMICA/PR3)")
            maps = hp.read_map(INPUT_FILE, field=None, hdu=1, verbose=False, memmap=True)
            
            # Extract
            print("  -> Extracting I...")
            extracted_data['I'] = np.array(maps[0][unique_pixels])
            print("  -> Extracting Q...")
            extracted_data['Q'] = np.array(maps[1][unique_pixels])
            print("  -> Extracting U...")
            extracted_data['U'] = np.array(maps[2][unique_pixels])
            
            del maps
            gc.collect()
            
        else:
            # SPLIT HDU LOAD (I in HDU 1, Q/U in HDU 2)
            print("  -> Mode: Split HDU (SEVEM/Commander/PR4)")
            
            # 1. Load I from HDU 1 (Field 0)
            print("  -> Loading I (HDU 1, Field 0)...")
            map_I = hp.read_map(INPUT_FILE, field=0, hdu=1, verbose=False, memmap=True)
            extracted_data['I'] = np.array(map_I[unique_pixels])
            del map_I
            gc.collect()
            print("     I extracted.")
            
            # 2. Load Q from HDU 2 (Field 0)
            print("  -> Loading Q (HDU 2, Field 0)...")
            map_Q = hp.read_map(INPUT_FILE, field=0, hdu=2, verbose=False, memmap=True)
            extracted_data['Q'] = np.array(map_Q[unique_pixels])
            del map_Q
            gc.collect()
            print("     Q extracted.")
            
            # 3. Load U from HDU 2 (Field 1)
            print("  -> Loading U (HDU 2, Field 1)...")
            map_U = hp.read_map(INPUT_FILE, field=1, hdu=2, verbose=False, memmap=True)
            extracted_data['U'] = np.array(map_U[unique_pixels])
            del map_U
            gc.collect()
            print("     U extracted.")

        # Validation
        if np.all(extracted_data['Q'] == 0) or np.all(extracted_data['U'] == 0):
             print("CRITICAL ERROR: Polarization data is all zeros. Check file or HDU logic.")
             sys.exit(1)
             
    except Exception as e:
        print(f"FATAL: Data load error: {e}")
        try:
            with fits.open(INPUT_FILE) as hdul:
                print("FITS Headers Diagnosis:")
                hdul.info()
        except: pass
        sys.exit(1)
            
    # 5. Prepare Tasks
    print("Phase 4: Remapping Indices...")
    tasks = []
    for i, (rid, theta, phi, rad) in enumerate(ring_definitions):
        original_pixels = ring_pixel_map[i]
        reduced_idxs = np.searchsorted(unique_pixels, original_pixels)
        tasks.append((rid, theta, phi, rad, reduced_idxs))
    
    del ring_pixel_map, unique_pixels
    gc.collect()
    
    # 6. Parallel Execution
    print(f"Phase 5: Processing {len(tasks)} rings...")
    num_processes = max(1, mp.cpu_count() - 1)
    
    results = []
    with mp.Pool(processes=num_processes, initializer=init_worker, 
                 initargs=(extracted_data['I'], extracted_data['Q'], extracted_data['U'])) as pool:
        for res in tqdm(pool.imap_unordered(process_ring_memory_optimized, tasks), total=len(tasks)):
            if res:
                results.append(res)
                
    # 7. Save
    print("Phase 6: Saving Results...")
    df = pd.DataFrame(results)
    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
    df.to_csv(OUTPUT_FILE, index=False)
    print(f"Done. Saved to {OUTPUT_FILE}")

if __name__ == "__main__":
    try:
        mp.set_start_method('fork', force=True)
    except RuntimeError: pass
    main()
