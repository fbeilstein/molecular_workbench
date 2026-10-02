import os

def _pack_orbital_group(zf, cube_dir, orbitals, group_name, color, label_fn):
    group = {'name': group_name, 'items': []}
    for info in orbitals:
        base_file = info['file'].replace('.cube', '.json')
        for sign in ['pos', 'neg']:
            sign_file = info['file'].replace('.cube', f'_{sign}.json')
            if os.path.exists(os.path.join(cube_dir, sign_file)):
                zf.write(os.path.join(cube_dir, sign_file), f'molecules/{sign_file}')
        group['items'].append({
            'label': label_fn(info),
            'file': f'molecules/{base_file}',
            'color': color,
        })
    return group

def package_orbitals(orb_manifest, cube_dir, zf):
    orbital_groups = []

    # 1. Frontier Orbitals
    frontier_group = {'name': 'Frontier', 'items': []}
    for lbl in ['homo', 'lumo']:
        if lbl in orb_manifest['canonical']:
            info = orb_manifest['canonical'][lbl]
            base_file = info['file'].replace('.cube', '.json')
            for sign in ['pos', 'neg']:
                sign_file = info['file'].replace('.cube', f'_{sign}.json')
                if os.path.exists(os.path.join(cube_dir, sign_file)):
                    zf.write(os.path.join(cube_dir, sign_file), f'molecules/{sign_file}')
            frontier_group['items'].append({
                'label': lbl.upper(),
                'file': f'molecules/{base_file}',
                'energy_ev': info.get('energy_ev'),
                'color': '#ff5c5c' if lbl == 'homo' else '#55aaff'
            })
    if frontier_group['items']:
        orbital_groups.append(frontier_group)

    # 2. Canonical π System
    pi_system = orb_manifest['canonical'].get('pi_system', [])
    if pi_system:
        pi_group = {'name': 'π System (Canonical)', 'items': []}
        PI_COLORS = ['#e040fb', '#ab47bc', '#7b1fa2', '#4a148c', '#ea80fc']
        for ci, info in enumerate(pi_system):
            base_file = info['file'].replace('.cube', '.json')
            for sign in ['pos', 'neg']:
                sign_file = info['file'].replace('.cube', f'_{sign}.json')
                if os.path.exists(os.path.join(cube_dir, sign_file)):
                    zf.write(os.path.join(cube_dir, sign_file), f'molecules/{sign_file}')
            pi_group['items'].append({
                'label': info['homo_label'],
                'file': f'molecules/{base_file}',
                'energy_ev': info.get('energy_ev'),
                'color': PI_COLORS[ci % len(PI_COLORS)]
            })
        orbital_groups.append(pi_group)

    # 3. Localized Orbitals
    loc = orb_manifest.get('localized', {})
    
    def _sort_bonds(bonds):
        return sorted(bonds, key=lambda x: (min(x['atoms']), max(x['atoms'])))
        
    if 'sigma' in loc: loc['sigma'] = _sort_bonds(loc['sigma'])
    if 'sigma_star' in loc: loc['sigma_star'] = _sort_bonds(loc['sigma_star'])
    if 'pi' in loc: loc['pi'] = _sort_bonds(loc['pi'])
    if 'pi_star' in loc: loc['pi_star'] = _sort_bonds(loc['pi_star'])
    
    if loc.get('sigma'):
        def sigma_label(info):
            if len(info['atoms']) > 2: return "deloc-σ(" + ",".join(info['atoms']) + ")"
            return f"σ({info['atoms'][0]}–{info['atoms'][1]})"
        grp = _pack_orbital_group(zf, cube_dir, loc['sigma'], 'σ Bonds', '#44cc77', sigma_label)
        orbital_groups.append(grp)
        
    if loc.get('pi'):
        def pi_label(info):
            if 'canonical_label' in info: return info['canonical_label']
            if len(info['atoms']) > 2: return "deloc-π(" + ",".join(info['atoms']) + ")"
            return f"π({info['atoms'][0]}={info['atoms'][1]})"
        grp = _pack_orbital_group(zf, cube_dir, loc['pi'], 'π Bonds', '#cc44bb', pi_label)
        orbital_groups.append(grp)
        
    if loc.get('sigma_star'):
        def sigmastar_label(info):
            if len(info['atoms']) > 2: return "σ*(" + ",".join(info['atoms']) + ")"
            return f"σ*({info['atoms'][0]}–{info['atoms'][1]})"
        grp = _pack_orbital_group(zf, cube_dir, loc['sigma_star'], 'σ* Anti-Bonds', '#ffaa00', sigmastar_label)
        orbital_groups.append(grp)
        
    if loc.get('pi_star'):
        def pistar_label(info):
            if 'canonical_label' in info: return info['canonical_label']
            if len(info['atoms']) > 2: return "π*(" + ",".join(info['atoms']) + ")"
            return f"π*({info['atoms'][0]}={info['atoms'][1]})"
        grp = _pack_orbital_group(zf, cube_dir, loc['pi_star'], 'π* Anti-Bonds', '#ff00aa', pistar_label)
        orbital_groups.append(grp)
        
    if loc.get('lone_pairs'):
        def lp_label(info):
            return f"LP({info['atoms'][0]})"
        grp = _pack_orbital_group(zf, cube_dir, loc['lone_pairs'], 'Lone Pairs', '#eeee22', lp_label)
        orbital_groups.append(grp)

    # 4. ESP surface if available
    esp_surface = None
    esp_info = orb_manifest.get('esp_surface')
    if esp_info:
        esp_file = esp_info['file']
        esp_path = os.path.join(cube_dir, esp_file)
        if os.path.exists(esp_path):
            zf.write(esp_path, f'molecules/{esp_file}')
            esp_surface = {
                'file': f'molecules/{esp_file}',
                'esp_min': esp_info['esp_min'],
                'esp_max': esp_info['esp_max'],
            }

    return orbital_groups, esp_surface
