def atomic_mobility_profile(dfc_prev, dfc_next, **kwargs):
    try:
        multip = kwargs['multip']
    except KeyError:
        multip = 1

    profile = {}

    atoms_prev = {atom['atom_id'] : atom for atom in dfc_prev.atoms}
    atoms_next = {atom['atom_id'] : atom for atom in dfc_next.atoms}
    xlo = (dfc_prev.xlo + dfc_next.xlo)/2
    xhi = (dfc_prev.xhi + dfc_next.xhi)/2
    ylo = (dfc_prev.ylo + dfc_next.ylo)/2
    yhi = (dfc_prev.yhi + dfc_next.yhi)/2
    zlo = (dfc_prev.zlo + dfc_next.zlo)/2
    zhi = (dfc_prev.zhi + dfc_next.zhi)/2
    lx = xhi - xlo
    ly = yhi - ylo
    lz = zhi - zlo

    # Find clay ranges
    clay_zlo = zhi
    clay_zhi = zlo

    for k, v in atoms_prev.items():
        if k in kwargs['clay_idcs']:
            clay_zlo = min(clay_zlo, v['z'])
            clay_zhi = max(clay_zhi, v['z'])

    # If clay crosses the box in z direction, clay ranges are incorrect
    clay_thickness_margined = 12
    if clay_zhi - clay_zlo > clay_thickness_margined:
        clay_zlo = dfc.zhi
        clay_zhi = dfc.zlo
        for k, v in atoms_prev.items():
            if k in clay_idcs:
                if abs(v['z'] - clay_zlo) < clay_thickness_margined:
                    clay_zlo = min(clay_zlo, v['z'])
                if abs(v['z'] - clay_zhi) < clay_thickness_margined:
                    clay_zhi = max(clay_zhi, v['z'])

    for k in atoms_prev.keys():
        a1 = atoms_prev[k]
        a2 = atoms_next[k]

        mid_x = mid_y = mid_z = None  # dr center

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

        dr = (dx**2 + dy**2 + dz**2)**0.5

        # distance to the clay surface:
        z = min([abs(mid_z - clay_zlo), abs(mid_z - clay_zhi),
                 abs(mid_z + lz - clay_zlo), abs(mid_z + lz - clay_zhi),
                 abs(mid_z - lz - clay_zlo), abs(mid_z - lz - clay_zhi)])
        z = round(z*multip, 0)

        try:
            profile[z][0] += dr
            profile[z][1] += 1
        except KeyError:
            profile[z] = [dr, 1]

    return profile
