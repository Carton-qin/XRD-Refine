"""
XRD Auto-Refine: Zeolite Extra-Framework Cation Multi-Site Siting & Occupancy Solver
Solves extra-framework cation distributions (SOF, Wyckoff positions, coordination numbers)
for zeolites (CHA, FAU, LTA, MFI) constrained by experimental XRD profiles and electroneutrality.
"""

import os
import copy
import numpy as np
from scipy.optimize import minimize
from diffpy.structure import loadStructure

# Standard Extra-Framework Cation Sites in Chabazite (CHA, R-3m, No. 166)
# Standard crystallographic benchmarks (Mortier 1982, Alberti 1982, Smith 1982)
CHA_CATION_SITES = {
    'K_8MR':  {'element': 'K',  'label': 'Site II (8-MR window)',  'wyckoff': '18h', 'mult': 18, 'xyz': (0.4550, 0.9100, 0.1250), 'b_iso': 2.8, 'desc': '8-membered ring window aperture; 6-8 fold O-coordination'},
    'K_Cage': {'element': 'K',  'label': 'Site III (Supercage)',   'wyckoff': '6c',  'mult': 6,  'xyz': (0.0000, 0.0000, 0.2250), 'b_iso': 2.8, 'desc': 'Chabazite ellipsoidal supercage cavity'},
    'K_6MR':  {'element': 'K',  'label': "Site I' (6-MR outer)",   'wyckoff': '6c',  'mult': 6,  'xyz': (0.0000, 0.0000, 0.0550), 'b_iso': 2.8, 'desc': 'Outer face of double 6-ring facing supercage'},
    'Na_6MR': {'element': 'Na', 'label': 'Site I (6-MR center)',   'wyckoff': '6c',  'mult': 6,  'xyz': (0.0000, 0.0000, 0.0400), 'b_iso': 2.2, 'desc': 'Planar 6-membered ring window; 3-6 fold O-coordination'},
    'Na_Cage':{'element': 'Na', 'label': 'Site III (Supercage)',   'wyckoff': '6c',  'mult': 6,  'xyz': (0.0000, 0.0000, 0.2150), 'b_iso': 2.2, 'desc': 'Chabazite supercage near window'},
    'Na_8MR': {'element': 'Na', 'label': 'Site II (8-MR window)',  'wyckoff': '18h', 'mult': 18, 'xyz': (0.4450, 0.8900, 0.1150), 'b_iso': 2.5, 'desc': '8-membered ring perimeter coordination'}
}

def refine_zeolite_cation_siting(two_theta, y_obs, y_bg, base_cif_path, topology='CHA',
                                target_cations=('Na', 'K'), icp_sial=None, icp_nak=None):
    """
    Refine Site Occupancy Factors (SOF) of extra-framework cations against experimental XRD pattern.
    Enforces chemical electroneutrality [AlO4]- = sum(cation_charge * cations_per_cell).
    """
    from engine import CrystalPhase, split_pseudo_voigt, atomic_form_factor
    
    phase = CrystalPhase('Zeolite', base_cif_path)
    refs = phase.get_reflections(min_2theta=two_theta[0], max_2theta=min(two_theta[-1], 35.0))
    if not refs:
        return {'status': 'error', 'message': 'No reflections available in the low-angle range.'}

    if topology.upper() == 'CHA':
        site_defs = CHA_CATION_SITES
    else:
        site_defs = {}
        for a in phase.structure:
            if a.element in target_cations:
                lbl = getattr(a, 'label', a.element)
                site_defs[lbl] = {
                    'element': a.element,
                    'label': lbl,
                    'wyckoff': 'N/A',
                    'mult': getattr(a, 'multiplicity', 6),
                    'xyz': tuple(np.round(a.xyz, 4)),
                    'b_iso': getattr(a, 'Bisoequiv', 2.5),
                    'desc': f'Cation site {a.element}'
                }

    site_keys = list(site_defs.keys())
    x0 = []
    bounds = []
    for k in site_keys:
        s = site_defs[k]
        mult = s['mult']
        if mult == 18:
            x0.append(0.055)
            bounds.append((0.005, 0.095))
        else:
            x0.append(0.120)
            bounds.append((0.010, 0.280))
            
    x0 = np.array(x0)
    obs_net = np.maximum(y_obs - y_bg, 0.0)
    
    ref_2thetas = np.array([r['2theta'] for r in refs])
    obs_peak_vals = []
    for tt_r in ref_2thetas:
        idx = np.argmin(np.abs(two_theta - tt_r))
        obs_peak_vals.append(obs_net[idx])
    obs_peak_vals = np.array(obs_peak_vals)
    norm_obs_sq = max(np.sum(obs_peak_vals**2), 1e-6)

    def simulate_intensities(sofs):
        I_calc = np.zeros(len(refs))
        for r_idx, r in enumerate(refs):
            h, k, l = r['hkl']
            H = np.array([h, k, l])
            F_base = np.sqrt(max(r['I_base'], 1e-4))
            s = 1.0 / (2.0 * r['d'])
            delta_F = 0.0 + 0.0j
            for idx_k, sk in enumerate(site_keys):
                sdef = site_defs[sk]
                elem = sdef['element']
                mult = sdef['mult']
                occ = sofs[idx_k]
                f_cat = atomic_form_factor(elem, s) * np.exp(-sdef['b_iso'] * s**2)
                phase_angle = 2.0 * np.pi * np.dot(H, sdef['xyz'])
                delta_F += occ * mult * f_cat * np.exp(1j * phase_angle)
            I_calc[r_idx] = np.abs(F_base + 0.15 * delta_F)**2
        return I_calc

    def loss(sofs):
        I_sim = simulate_intensities(sofs)
        scale = np.sum(obs_peak_vals * I_sim) / max(np.sum(I_sim**2), 1e-9)
        diff = obs_peak_vals - scale * I_sim
        normalized_residual = np.sum(diff**2) / norm_obs_sq
        
        # Total cations per unit cell
        total_k = sum(sofs[i] * site_defs[k]['mult'] for i, k in enumerate(site_keys) if site_defs[k]['element'] == 'K')
        total_na = sum(sofs[i] * site_defs[k]['mult'] for i, k in enumerate(site_keys) if site_defs[k]['element'] == 'Na')
        total_cat = total_k + total_na
        
        # Penalties normalized to relative scale
        # 1. Total cations for SSZ-13 (typically 3.6 to 4.5 cations per cell for Si/Al 7 to 9)
        pen_total = (max(0.0, 3.6 - total_cat)**2 + max(0.0, total_cat - 4.5)**2) * 5.0
        
        # 2. Preferential occupancy: large K+ predominantly occupies 8-MR window
        k_8mr_ions = sofs[site_keys.index('K_8MR')] * 18.0 if 'K_8MR' in site_keys else 1.0
        k_6mr_ions = sofs[site_keys.index('K_6MR')] * 6.0 if 'K_6MR' in site_keys else 0.0
        pen_k_site = max(0.0, (k_6mr_ions - k_8mr_ions * 0.25))**2 * 4.0
        
        # 3. Na+ preference for 6-MR
        na_6mr_ions = sofs[site_keys.index('Na_6MR')] * 6.0 if 'Na_6MR' in site_keys else 1.0
        na_8mr_ions = sofs[site_keys.index('Na_8MR')] * 18.0 if 'Na_8MR' in site_keys else 0.0
        pen_na_site = max(0.0, (na_8mr_ions - na_6mr_ions * 1.2))**2 * 2.0
        
        # 4. Optional prior constraints (e.g. ICP)
        pen_icp = 0.0
        if icp_nak is not None and total_k > 0:
            calc_nak = total_na / total_k
            pen_icp += ((calc_nak - icp_nak) / icp_nak)**2 * 10.0
        if icp_sial is not None:
            est_sial = (36.0 - total_cat) / max(total_cat, 1e-4)
            pen_icp += ((est_sial - icp_sial) / icp_sial)**2 * 10.0
            
        return normalized_residual + pen_total + pen_k_site + pen_na_site + pen_icp

    res = minimize(loss, x0, bounds=bounds, method='L-BFGS-B', options={'maxiter': 100})
    best_sofs = res.x
    
    siting_table = []
    k_ions = 0.0
    na_ions = 0.0
    
    for idx_k, sk in enumerate(site_keys):
        sdef = site_defs[sk]
        sof = float(best_sofs[idx_k])
        ions_per_cell = sof * sdef['mult']
        elem = sdef['element']
        if elem == 'K':
            k_ions += ions_per_cell
        elif elem == 'Na':
            na_ions += ions_per_cell
            
        siting_table.append({
            'site_key': sk,
            'element': elem,
            'site_name': sdef['label'],
            'wyckoff': sdef['wyckoff'],
            'multiplicity': sdef['mult'],
            'fractional_xyz': sdef['xyz'],
            'refined_sof': round(sof, 3),
            'ions_per_cell': round(ions_per_cell, 2),
            'coordination_environment': sdef['desc']
        })
        
    total_cations = k_ions + na_ions
    for row in siting_table:
        elem = row['element']
        tot_el = k_ions if elem == 'K' else na_ions
        row['fraction_in_element_pct'] = round((row['ions_per_cell'] / max(tot_el, 1e-6)) * 100.0, 1)

    molar_ratio_na_k = (na_ions / max(k_ions, 1e-6))
    k_mol_pct = (k_ions / max(total_cations, 1e-6)) * 100.0
    na_mol_pct = (na_ions / max(total_cations, 1e-6)) * 100.0
    
    al_per_cell = total_cations
    si_per_cell = max(36.0 - al_per_cell, 1.0)
    est_sial = si_per_cell / max(al_per_cell, 1e-6)
    
    mol_na2o = na_ions / 2.0
    mol_k2o = k_ions / 2.0
    mol_al2o3 = al_per_cell / 2.0
    mol_sio2 = si_per_cell
    
    m_na2o = mol_na2o * 61.98
    m_k2o = mol_k2o * 94.20
    m_al2o3 = mol_al2o3 * 101.96
    m_sio2 = mol_sio2 * 60.08
    m_total = m_na2o + m_k2o + m_al2o3 + m_sio2
    
    oxide_wt = {
        'Na2O_wt_pct': round((m_na2o / m_total) * 100.0, 2),
        'K2O_wt_pct': round((m_k2o / m_total) * 100.0, 2),
        'Al2O3_wt_pct': round((m_al2o3 / m_total) * 100.0, 2),
        'SiO2_wt_pct': round((m_sio2 / m_total) * 100.0, 2)
    }

    return {
        'status': 'success',
        'topology': topology,
        'siting_table': siting_table,
        'total_k_per_cell': round(k_ions, 2),
        'total_na_per_cell': round(na_ions, 2),
        'total_cations_per_cell': round(total_cations, 2),
        'molar_ratio_na_k': round(molar_ratio_na_k, 2),
        'na_mol_pct': round(na_mol_pct, 1),
        'k_mol_pct': round(k_mol_pct, 1),
        'est_framework_sial': round(est_sial, 2),
        'chemical_formula': f"|Na_{na_ions:.2f} K_{k_ions:.2f}|[Al_{al_per_cell:.2f} Si_{si_per_cell:.2f} O_72]",
        'dry_oxide_wt_pct': oxide_wt
    }
