from utils.single_molecule_gyration_tensor import single_molecule_gyration_tensor


def gyration_composite(dfc, **kwargs):
    try:
        multip = kwargs['multip']
    except KeyError:
        multip = 1.

    lx = dfc.xhi - dfc.xlo
    ly = dfc.yhi - dfc.ylo
    lz = dfc.zhi - dfc.zlo
    result = {}

    # Distribute atoms into their molecules
    molecules = {}
    for atom in dfc.atoms:
        k = kwargs['poly_mol_idx'](atom['atom_id'])
        try:
            molecules[k].append(atom)
        except KeyError:
            if k is None:
                continue
            molecules[k] = [atom]

    # Find clay ranges
    clay_zlo = dfc.zhi
    clay_zhi = dfc.zlo
    for atom in dfc.atoms:
        if atom['atom_id'] in kwargs['clay_idcs']:
            clay_zlo = min(clay_zlo, atom['z'])
            clay_zhi = max(clay_zhi, atom['z'])

    # If clay crosses the box in z direction, clay ranges are incorrect
    clay_thickness_margined = 12
    if clay_zhi - clay_zlo > clay_thickness_margined:
        clay_zlo = dfc.zhi
        clay_zhi = dfc.zlo
        for atom in dfc.atoms:
            if atom['atom_id'] in kwargs['clay_idcs']:
                if abs(atom['z'] - clay_zlo) < clay_thickness_margined:
                    clay_zlo = min(clay_zlo, atom['z'])
                if abs(atom['z'] - clay_zhi) < clay_thickness_margined:
                    clay_zhi = max(clay_zhi, atom['z'])

    for k, v in molecules.items():
        tmp_result = single_molecule_gyration_tensor(dfc, v, kwargs['masses'])
        # distance to the clay surface:
        z = tmp_result['z_cm']
        z = min([abs(z - clay_zlo), abs(z - clay_zhi),
                 abs(z + lz - clay_zlo), abs(z + lz - clay_zhi),
                 abs(z - lz - clay_zlo), abs(z - lz - clay_zhi)])
        z = round(z*multip, 0)
        try:
            result[z]['rx'] += tmp_result['x']
            result[z]['ry'] += tmp_result['y']
            result[z]['rz'] += tmp_result['z']
            result[z]['count'] += 1
        except KeyError:
            result[z] = {
                'rx': tmp_result['x'],
                'ry': tmp_result['y'],
                'rz': tmp_result['z'],
                'count': 1
            }

    for v in result.values():
        v['rx'] /= v['count']
        v['ry'] /= v['count']
        v['rz'] /= v['count']

    return result
