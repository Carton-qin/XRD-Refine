"""
XRD Auto-Refine: Peak Search & Candidate Crystal Phase Matcher
Uses Hanawalt 3-strong-peak indexing and Figure-of-Merit (FOM) cross-correlation.
"""

import os
import glob
import numpy as np
from scipy.signal import find_peaks, savgol_filter
from engine import CrystalPhase

def extract_observed_peaks(two_theta, y_obs, min_prominence_ratio=0.03):
    """
    Smooth pattern and detect significant diffraction peaks.
    Returns: list of dicts with '2theta', 'd', 'intensity', 'prominence'.
    """
    two_theta = np.array(two_theta, dtype=float)
    y_obs = np.array(y_obs, dtype=float)
    
    # Smooth data to suppress Poisson counting noise
    window_length = min(15, len(y_obs) if len(y_obs) % 2 != 0 else len(y_obs) - 1)
    if window_length >= 5:
        y_smooth = savgol_filter(y_obs, window_length=window_length, polyorder=2)
    else:
        y_smooth = y_obs
        
    # Baseline estimation
    baseline = np.percentile(y_smooth, 10)
    net_intensity = np.maximum(y_smooth - baseline, 0.0)
    max_peak_val = np.max(net_intensity) if len(net_intensity) > 0 else 1.0
    prominence_thresh = max_peak_val * min_prominence_ratio
    
    # Find peaks
    peaks_idx, props = find_peaks(net_intensity, prominence=prominence_thresh, distance=5)
    
    wavelength = 1.5406  # Cu K-alpha
    detected_peaks = []
    for idx, prom in zip(peaks_idx, props['prominences']):
        tt = two_theta[idx]
        th_rad = np.radians(tt / 2.0)
        d_val = wavelength / (2.0 * np.sin(th_rad)) if np.sin(th_rad) > 0 else 0.0
        detected_peaks.append({
            '2theta': float(tt),
            'd': float(d_val),
            'intensity': float(net_intensity[idx]),
            'prominence': float(prom)
        })
        
    detected_peaks.sort(key=lambda x: x['intensity'], reverse=True)
    return detected_peaks

def match_candidate_cifs(two_theta, y_obs, db_dir, chemical_elements=None, tolerance_deg=0.25):
    """
    Match observed diffraction peaks against all CIF files in db_dir.
    Computes Hanawalt match and continuous FOM cross-correlation.
    """
    obs_peaks = extract_observed_peaks(two_theta, y_obs)
    if not obs_peaks:
        return []
        
    obs_tt = [p['2theta'] for p in obs_peaks]
    obs_d  = [p['d'] for p in obs_peaks]
    obs_I  = [p['intensity'] for p in obs_peaks]
    max_obs_I = max(obs_I) if obs_I else 1.0
    obs_rel_I = [i / max_obs_I * 100.0 for i in obs_I]
    
    cif_files = glob.glob(os.path.join(db_dir, '**', '*.cif'), recursive=True)
    results = []
    
    for cif_path in cif_files:
        pname = os.path.splitext(os.path.basename(cif_path))[0]
        category = os.path.basename(os.path.dirname(cif_path))
        
        try:
            phase = CrystalPhase(pname, cif_path)
        except Exception:
            continue
            
        # Check chemical filter if provided
        comp_elements = set(phase.structure.composition.keys())
        if chemical_elements:
            chem_set = set(chemical_elements)
            # If phase contains foreign metals not in the chemical system, lower its priority
            system_metals = {'Al', 'Si', 'O', 'H'}.union(chem_set)
            foreign_metals = {el for el in comp_elements if el not in system_metals}
            if foreign_metals:
                chem_penalty = 0.5
            else:
                chem_penalty = 1.0
        else:
            chem_penalty = 1.0
            
        # Get theoretical reflections in range
        refs = phase.get_reflections(min_2theta=two_theta[0] - 0.5, max_2theta=two_theta[-1] + 0.5)
        if not refs:
            continue
            
        refs.sort(key=lambda x: x['I_base'], reverse=True)
        max_ref_I = refs[0]['I_base'] if refs else 1.0
        
        # Match top 8 theoretical reflections against observed peaks
        top_refs = refs[:8]
        matched_peaks = []
        total_score = 0.0
        
        for r in top_refs:
            r_tt = r['2theta']
            r_rel_I = r['I_base'] / max_ref_I * 100.0
            
            # Find closest observed peak
            deltas = [abs(r_tt - o_tt) for o_tt in obs_tt]
            min_delta_idx = int(np.argmin(deltas))
            min_delta = deltas[min_delta_idx]
            
            if min_delta <= tolerance_deg:
                o_rel = obs_rel_I[min_delta_idx]
                intensity_weight = min(r_rel_I, o_rel) / max(r_rel_I, o_rel, 1.0)
                proximity_weight = np.exp(-0.5 * (min_delta / (tolerance_deg / 2.0))**2)
                score = (r_rel_I / 100.0) * proximity_weight * (0.5 + 0.5 * intensity_weight)
                total_score += score
                matched_peaks.append({
                    'hkl': r['hkl'],
                    'calc_2theta': round(r_tt, 3),
                    'obs_2theta': round(obs_tt[min_delta_idx], 3),
                    'delta_deg': round(obs_tt[min_delta_idx] - r_tt, 3),
                    'calc_rel_I': round(r_rel_I, 1),
                    'obs_rel_I': round(o_rel, 1)
                })
                
        # Normalize FOM score (0 to 100)
        fom = (total_score / sum(r['I_base'] / max_ref_I for r in top_refs)) * 100.0 * chem_penalty
        fom = min(max(fom, 0.0), 100.0)
        
        # Categorize into Recommendation Tier
        if fom >= 65.0 and len(matched_peaks) >= 3:
            tier = 'Tier A (High Confidence / Primary Candidate)'
        elif fom >= 35.0 and len(matched_peaks) >= 2:
            tier = 'Tier B (Plausible / Intermediate / Secondary Candidate)'
        else:
            tier = 'Tier C (Low Support / Unlikely / Trace)'
            
        results.append({
            'name': pname,
            'category': category,
            'cif_path': cif_path,
            'fom': round(fom, 1),
            'tier': tier,
            'matched_peaks_count': len(matched_peaks),
            'matched_peaks': matched_peaks,
            'formula': phase.structure.composition
        })
        
    results.sort(key=lambda x: x['fom'], reverse=True)
    return results
