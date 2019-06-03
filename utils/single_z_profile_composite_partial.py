from collections import Counter


def single_z_profile_partial(dfc, **kwargs):
    '''
    Get density profile from a single file of a composite
    possible kwargs are:
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

    for atom in dfc.atoms:
        if atom['atom_id'] in clay_idcs:
            phase = 'clay'
        elif atom['atom_id'] in kwargs['mod_idcs']:
            phase = 'modifier'
        else:
            phase = 'polymer'

        # mass_value would be in g/cm3 while masses[...] are in a.m.u.
        mass_value = masses[atom['atom_type_id']] / lx / ly / dz * 1.66053886

        # distance to the clay surface:
        #z = min([abs(atom['z'] - clay_zlo), abs(atom['z'] - clay_zhi),
        #    abs(atom['z'] + lz - clay_zlo), abs(atom['z'] + lz - clay_zhi),
        #    abs(atom['z'] - lz - clay_zlo), abs(atom['z'] - lz - clay_zhi)])

        z = round(atom['z'] * multip, 0)

        try:
            result[z][phase] += mass_value
        except KeyError:
            try:
                result[z][phase] = mass_value
            except:
                result[z] = {'clay': 0, 'modifier': 0, 'polymer': 0}
                result[z][phase] += mass_value

        if 'replicate' in kwargs.keys() and kwargs['replicate']:
            z_minus = round((atom['z'] - lz) * multip, 0)
            z_plus = round((atom['z'] + lz) * multip, 0)

            try:
                result[z_minus][phase] += mass_value
            except KeyError:
                try:
                    result[z_minus][phase] = mass_value
                except:
                    result[z_minus] = {'clay': 0, 'modifier': 0, 'polymer': 0}
                    result[z_minus][phase] += mass_value

            try:
                result[z_plus][phase] += mass_value
            except KeyError:
                try:
                    result[z_plus][phase] = mass_value
                except:
                    result[z_plus] = {'clay': 0, 'modifier': 0, 'polymer': 0}
                    result[z_plus][phase] += mass_value

    return result

