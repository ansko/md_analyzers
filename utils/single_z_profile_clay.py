from collections import Counter


def single_z_profile(dfc, **kwargs):
    '''
    Get density profile from a single file of a clay.
    Possible kwargs are:
        multip    -- more points at x axis
        masses    -- calculate mass density if set
        replicate -- replicate the structure in z direction (nz in[1, -1])
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

    lx = dfc.xhi - dfc.xlo
    ly = dfc.yhi - dfc.ylo
    lz = dfc.zhi - dfc.zlo
    dz = 1.0 / multip
    result = Counter()

    for atom in dfc.atoms:
        if masses is None:
            mass_value = 1  # just counting
            print('a')
        else:
            # mass_value would be in g/cm3 while masses[...] are in a.m.u.
            mass_value = masses[atom['atom_type_id']] / lx / ly / dz * 1.66053886
        result[round(atom['z'] * multip, 0)] += mass_value
        if replicate:
            result[round((atom['z'] + lz) * multip, 0)] += mass_value
            result[round((atom['z'] - lz) * multip, 0)] += mass_value

    return result
