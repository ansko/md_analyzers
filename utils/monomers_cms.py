def monomers_cms(dfc, masses, monomer_idx):
    lx = dfc.xhi - dfc.xlo
    ly = dfc.yhi - dfc.ylo
    lz = dfc.zhi - dfc.zlo

    result = {}
    monomers = {}
    for atom in dfc.atoms:
        idx = monomer_idx(atom['atom_id'])
        if idx is None:
            continue
        try:
            monomers[idx].append(atom)
        except KeyError:
            monomers[idx] = [atom]

    for k, v in monomers.items():
        x_cm = y_cm = z_cm = M = 0
        for atom in v:
            m = masses[atom['atom_type_id']]
            M += m
            x_cm += (atom['x'] + atom['nx'] * lx) * m
            y_cm += (atom['y'] + atom['ny'] * ly) * m
            z_cm += (atom['z'] + atom['nz'] * lz) * m
        x_cm /= M
        y_cm /= M
        z_cm /= M
        result[k] = {'x_cm': x_cm, 'y_cm': y_cm, 'z_cm': z_cm}

    return result
