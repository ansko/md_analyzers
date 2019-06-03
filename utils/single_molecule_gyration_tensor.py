import copy


def single_molecule_gyration_tensor(dfc, atoms, masses, output=False):
    xlo = dfc.xlo;    xhi = dfc.xhi;    lx = xhi - xlo
    ylo = dfc.ylo;    yhi = dfc.yhi;    ly = yhi - ylo
    zlo = dfc.zlo;    zhi = dfc.zhi;    lz = zhi - zlo

    # get center-of-mass
    x_cm = y_cm = z_cm = m = 0
    for atom in atoms:
        x_cm += (atom['x'] + atom['nx']*lx) * masses[atom['atom_type_id']]
        y_cm += (atom['y'] + atom['ny']*ly) * masses[atom['atom_type_id']]
        z_cm += (atom['z'] + atom['nz']*lz) * masses[atom['atom_type_id']]
        m += masses[atom['atom_type_id']]
    x_cm /= m
    y_cm /= m
    z_cm /= m

    # get rgyrs
    # R = sqrt{1/N sum(r**2)}
    rx = ry = rz = 0
    for atom in atoms:
        dx = abs(atom['x'] + atom['nx']*lx - x_cm)
        dy = abs(atom['y'] + atom['ny']*ly - y_cm)
        dz = abs(atom['z'] + atom['nz']*lz - z_cm)

        rx += dx**2
        ry += dy**2
        rz += dz**2

    rx /= len(atoms)
    ry /= len(atoms)
    rz /= len(atoms)

    if z_cm < zlo:
        z_cm += lz
    elif z_cm > zhi:
        z_cm -= lz

    return {
        'x': rx**0.5,
        'y': ry**0.5,
        'z': rz**0.5,
        'x_cm': x_cm, 'y_cm': y_cm, 'z_cm': z_cm
    }
