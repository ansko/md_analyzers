def single_z_profile(dfc, **kwargs):
    '''
    '''
    try:
        multip = kwargs['multip']
    except KeyError:
        multip = 1.0
    try:
        masses = kwargs['masses']
    except:
        masses = None
    try:
        replicate = kwargs['replicate']
    except KeyError:
        replicate = False
    try:
        clay_idcs = kwargs['clay_idcs']
    except KeyError:
        clay_idcs = set()
    try:
        ignored_bond_types = kwargs['ignored_bond_types']
    except KeyError:
        ignored_bond_types = []

    # Find clay ranges
    clay_zlo = dfc.zhi
    clay_zhi = dfc.zlo
    for atom in dfc.atoms:
        if atom['atom_id'] in clay_idcs:
            clay_zlo = min(clay_zlo, atom['z'])
            clay_zhi = max(clay_zhi, atom['z'])

    # If clay crosses the box in z direction, clay ranges are incorrect
    clay_thickness_margined = 12
    if clay_zhi - clay_zlo > clay_thickness_margined:
        clay_zlo = dfc.zhi
        clay_zhi = dfc.zlo
        for atom in dfc.atoms:
            if atom['atom_id'] in clay_idcs:
                if abs(atom['z'] - clay_zlo) < clay_thickness_margined:
                    clay_zlo = min(clay_zlo, atom['z'])
                if abs(atom['z'] - clay_zhi) < clay_thickness_margined:
                    clay_zhi = max(clay_zhi, atom['z'])

    lx = dfc.xhi - dfc.xlo
    ly = dfc.yhi - dfc.ylo
    lz = dfc.zhi - dfc.zlo
    dz = 2  # since key is distance to the closest clay surface
    if multip is not None:
        dz /= multip

    result = {}

    atoms = {a['atom_id']: a for a in dfc.atoms}
    bonds = {b['bond_id']: b for b in dfc.bonds}

    for bond in bonds.values():
        if bond['bond_type_id'] in ignored_bond_types:
            continue

        a1 = atoms[bond['atom_one_id']]
        a2 = atoms[bond['atom_two_id']]
        mid_x = mid_y = mid_z = None  # bond center

        dx = a1['x'] - a2['x']
        if abs(dx) > lx - abs(dx):
            dx = lx - abs(dx)
            mid_x = (a1['x'] + a2['x'])/2 + lx/2
        else:
            mid_x = (a1['x'] + a2['x'])/2
        dy = a1['y'] - a2['y']
        if abs(dy) > ly - abs(dy):
            dy = ly - abs(dy)
            mid_y = (a1['y'] + a2['y'])/2 + ly/2
        else:
            mid_y = (a1['y'] + a2['y'])/2
        dz = a1['z'] - a2['z']
        if abs(dz) > lz - abs(dz):
            dz = lz - abs(dz)
            mid_z = (a1['z'] + a2['z'])/2 + lz/2
        else:
            mid_z = (a1['z'] + a2['z'])/2

        cos_theta2 = dz**2 / (dx**2 + dy**2 + dz**2)
        P = 0.5 * (3 * cos_theta2 - 1)


        #mid_z = a1['z']
        # distance to the clay surface:
        z = min([abs(mid_z - clay_zlo), abs(mid_z - clay_zhi),
                 abs(mid_z + lz - clay_zlo), abs(mid_z + lz - clay_zhi),
                 abs(mid_z - lz - clay_zlo), abs(mid_z - lz - clay_zhi)])
        z = round(z*multip, 0)

        try:
            result[z][0] += P
            result[z][1] += 1
        except KeyError:
            result[z] = [P, 1]

    return result
