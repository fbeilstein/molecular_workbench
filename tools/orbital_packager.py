import os
import re

BONDING_COLOR = '#38bdf8'
ANTIBONDING_COLOR = '#fb923c'

_ATOM_RE = re.compile(r'^([A-Za-z]+)(\d+)$')


def _atom_sort_key(atom_id: str):
    m = _ATOM_RE.match(str(atom_id))
    if m:
        return (int(m.group(2)), m.group(1))
    return (999999, str(atom_id))


def _with_frontier_tag(base_label: str, info: dict) -> str:
    tag = info.get('frontier_tag')
    if tag and f"({tag})" not in base_label:
        return f"{base_label} ({tag})"
    return base_label


def _pack_orbital_group(zf, cube_dir, orbitals, group_name, color, label_fn):
    group = {'name': group_name, 'items': []}
    for info in orbitals:
        base_file = info['file'].replace('.cube', '.json')
        for sign in ['pos', 'neg']:
            sign_file = info['file'].replace('.cube', f'_{sign}.json')
            if os.path.exists(os.path.join(cube_dir, sign_file)):
                zf.write(os.path.join(cube_dir, sign_file), f'molecules/{sign_file}')
        item = {
            'label': _with_frontier_tag(label_fn(info), info),
            'file': f'molecules/{base_file}',
            'color': color,
        }
        if info.get('energy_ev') is not None:
            item['energy_ev'] = round(float(info['energy_ev']), 2)
        group['items'].append(item)
    return group


def package_orbitals(orb_manifest, cube_dir, zf):
    orbital_groups = []

    # 1. Frontier Orbitals (supports default HOMO/LUMO as well as expanded HOMO-k .. LUMO+k)
    frontier_group = {'name': 'Frontier', 'items': []}
    canonical = orb_manifest.get('canonical', {})
    frontier_entries = canonical.get('frontier')
    if not frontier_entries:
        frontier_entries = []
        for lbl in ['homo', 'lumo']:
            if lbl in canonical:
                entry = dict(canonical[lbl])
                entry.setdefault('canonical_label', lbl.upper())
                frontier_entries.append(entry)

    for info in frontier_entries:
        base_file = info['file'].replace('.cube', '.json')
        for sign in ['pos', 'neg']:
            sign_file = info['file'].replace('.cube', f'_{sign}.json')
            if os.path.exists(os.path.join(cube_dir, sign_file)):
                zf.write(os.path.join(cube_dir, sign_file), f'molecules/{sign_file}')
        lbl = info.get('canonical_label') or info.get('frontier_key', '').upper()
        is_anti = ('LUMO' in lbl.upper()) or ('lumo' in info.get('file', '').lower())
        item = {
            'label': lbl,
            'file': f'molecules/{base_file}',
            'color': ANTIBONDING_COLOR if is_anti else BONDING_COLOR,
        }
        if info.get('energy_ev') is not None:
            item['energy_ev'] = round(float(info['energy_ev']), 2)
        frontier_group['items'].append(item)

    if frontier_group['items']:
        orbital_groups.append(frontier_group)

    # 2. Legacy Canonical π System (if present)
    pi_system = canonical.get('pi_system', [])
    if pi_system:
        pi_group = {'name': 'π System (Canonical)', 'items': []}
        for info in pi_system:
            base_file = info['file'].replace('.cube', '.json')
            for sign in ['pos', 'neg']:
                sign_file = info['file'].replace('.cube', f'_{sign}.json')
                if os.path.exists(os.path.join(cube_dir, sign_file)):
                    zf.write(os.path.join(cube_dir, sign_file), f'molecules/{sign_file}')
            pi_group['items'].append({
                'label': _with_frontier_tag(info['homo_label'], info),
                'file': f'molecules/{base_file}',
                'energy_ev': info.get('energy_ev'),
                'color': BONDING_COLOR,
            })
        orbital_groups.append(pi_group)

    # 3. Chemical Orbitals (Localized σ/π/LP + Delocalized Planar π/π*)
    loc = orb_manifest.get('localized', {})

    def _sort_sigma(bonds):
        return sorted(
            bonds,
            key=lambda x: (
                tuple(sorted(_atom_sort_key(a) for a in x.get('atoms', []))),
                float(x.get('energy_ev') if x.get('energy_ev') is not None else 0.0),
            ),
        )

    def _sort_pi(bonds):
        return sorted(
            bonds,
            key=lambda x: (
                0 if x.get('canonical_label') else 1,
                tuple(_atom_sort_key(a) for a in x.get('atoms', [])) if not x.get('canonical_label') else (tuple(_atom_sort_key(a) for a in x.get('atoms', [])),),
                float(x.get('energy_ev') if x.get('energy_ev') is not None else 0.0),
            ),
        )

    def _sort_lps(lps):
        return sorted(
            lps,
            key=lambda x: (
                tuple(_atom_sort_key(a) for a in x.get('atoms', [])),
                -float(x.get('energy_ev') if x.get('energy_ev') is not None else 0.0),
            ),
        )

    if 'sigma' in loc:
        loc['sigma'] = _sort_sigma(loc['sigma'])
    if 'sigma_star' in loc:
        loc['sigma_star'] = _sort_sigma(loc['sigma_star'])
    if 'pi' in loc:
        loc['pi'] = _sort_pi(loc['pi'])
    if 'pi_star' in loc:
        loc['pi_star'] = _sort_pi(loc['pi_star'])
    if 'lone_pairs' in loc:
        loc['lone_pairs'] = _sort_lps(loc['lone_pairs'])

    if loc.get('sigma'):
        def sigma_label(info):
            if len(info['atoms']) > 2:
                return "deloc-σ(" + ",".join(info['atoms']) + ")"
            return f"σ({info['atoms'][0]}–{info['atoms'][1]})"
        grp = _pack_orbital_group(zf, cube_dir, loc['sigma'], 'σ Bonds', BONDING_COLOR, sigma_label)
        orbital_groups.append(grp)

    if loc.get('pi'):
        def pi_label(info):
            if 'canonical_label' in info:
                return info['canonical_label']
            if len(info['atoms']) > 2:
                return "deloc-π(" + ",".join(info['atoms']) + ")"
            return f"π({info['atoms'][0]}={info['atoms'][1]})"
        grp = _pack_orbital_group(zf, cube_dir, loc['pi'], 'π Bonds', BONDING_COLOR, pi_label)
        orbital_groups.append(grp)

    if loc.get('sigma_star'):
        def sigmastar_label(info):
            if len(info['atoms']) > 2:
                return "σ*(" + ",".join(info['atoms']) + ")"
            return f"σ*({info['atoms'][0]}–{info['atoms'][1]})"
        grp = _pack_orbital_group(zf, cube_dir, loc['sigma_star'], 'σ* Anti-Bonds', ANTIBONDING_COLOR, sigmastar_label)
        orbital_groups.append(grp)

    if loc.get('pi_star'):
        def pistar_label(info):
            if 'canonical_label' in info:
                return info['canonical_label']
            if len(info['atoms']) > 2:
                return "π*(" + ",".join(info['atoms']) + ")"
            return f"π*({info['atoms'][0]}={info['atoms'][1]})"
        grp = _pack_orbital_group(zf, cube_dir, loc['pi_star'], 'π* Anti-Bonds', ANTIBONDING_COLOR, pistar_label)
        orbital_groups.append(grp)

    if loc.get('lone_pairs'):
        def lp_label(info):
            return f"LP({info['atoms'][0]})"
        grp = _pack_orbital_group(zf, cube_dir, loc['lone_pairs'], 'Lone Pairs', BONDING_COLOR, lp_label)
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
