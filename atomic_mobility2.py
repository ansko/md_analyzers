from collections import Counter
import copy
import math
import os
import sys

from md_work.structural_tools.datafile_content import DatafileContent
from utils.atomic_mobility_profile import atomic_mobility_profile

if __name__ == '__main__':
    multip = 0.25

    system_name = sys.argv[1]

    fnames = clay_idcs = plot_kwargs = None

    if system_name == 'mmt':
        raise NotImplementedError()
    elif system_name == 'L':
        from utils.monomers_cms import monomers_cms
        from plotters.plot_data import plot_data
        data_dir = ('/media/anton/Seagate Expansion Drive/Backups/'
                    'Cluster calculations for article/BiggerSystems/Comp/'
                    '10chains/2.2 - More relaxation 500 (wiggle)/'
                    '1785842 - relaxation')
        fnames = ['{0}/co.{1}.data'.format(data_dir, idx * 50000)
            for idx in range(1, 51)]
        clay_idcs = set()
        for idx in range(9):
            clay_idcs.update(set(range(1 + idx*5380, 721 + idx*5380)))
    elif system_name == 'S':  # Too small data
        raise NotImplementedError()
    elif system_name == 'mix':
        raise NotImplementedError()
    elif system_name == 'seg':
        raise NotImplementedError()
    elif system_name == 'pa6x20':
        from utils.monomers_cms import monomers_cms
        from plotters.plot_data import plot_data
        data_dir = ('/media/anton/Seagate Expansion Drive/'
                    'md_new_data/pa6x20/3507904_500K_pz100/datafiles')
        fnames = ['{0}/pa6x20.{1}.data'.format(data_dir, idx * 50000)
            for idx in range(1, 51)]
        clay_idcs = set()
    else:
        print('Unknown system name:', sys.argv[1])
        sys.exit()

    dfc_start = DatafileContent(fnames[0])
    dfc_end = DatafileContent(fnames[-1])
    lx = dfc_start.xhi - dfc_start.xlo
    ly = dfc_start.yhi - dfc_start.ylo
    lz = dfc_start.zhi - dfc_start.zlo
    start_positions = { 
        a['atom_id']: {
            'x': a['x'] + a['nx']*lx,
            'y': a['y'] + a['ny']*ly,
            'z': a['z'] + a['nz']*lz
        } for a in dfc_start.atoms
    }
    end_positions = { 
        a['atom_id']: {
            'x': a['x'] + a['nx']*lx,
            'y': a['y'] + a['ny']*ly,
            'z': a['z'] + a['nz']*lz
        } for a in dfc_end.atoms
    }

    # Find clay ranges
    clay_zlo = dfc_start.zhi
    clay_zhi = dfc_start.zlo

    for k, v in start_positions.items():
        if k in clay_idcs:
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

    profile = {}
    for k, v in start_positions.items():
        a1 = start_positions[k]
        a2 = end_positions[k]

        dx = a1['x'] - a2['x']# + (a1['nx'] - a2['nx']) * lx
        dy = a1['y'] - a2['y']# + (a1['ny'] - a2['ny']) * ly
        dz = a1['z'] - a2['z']# + (a1['nz'] - a2['nz']) * lz
        #if abs(dx) > lx - abs(dx):
        #    dx = lx - abs(dx)
        #dy = a1['y'] - a2['y']
        #if abs(dy) > ly - abs(dy):
        #    dy = ly - abs(dy)
        #dz = a1['z'] - a2['z']
        #if abs(dz) > lz - abs(dz):
        #    dz = lz - abs(dz)

        dr = (dx**2 + dy**2 + dz**2)**0.5

        # distance to the clay surface:
        z = min([abs(a1['z'] - clay_zlo), abs(a1['z'] - clay_zhi),
                 abs(a1['z'] + lz - clay_zlo), abs(a1['z'] + lz - clay_zhi),
                 abs(a1['z'] - lz - clay_zlo), abs(a1['z'] - lz - clay_zhi)])
        z = round(z*multip, 0)

        try:
            profile[z][0] += dr
            profile[z][1] += 1
        except KeyError:
            profile[z] = [dr, 1]


    if system_name == 'L':
        f = open('logs/atomic_mobility_L', 'w')
        for k in sorted(profile.keys()):
            print(k, profile[k][0] / profile[k][1], file=f)
        xs = sorted(profile.keys())
        ys = [profile[x][0] / profile[x][1] for x in xs]
        xs = [x / multip for x in xs]
        plot_data(xs, ys=ys, out_fname='atomic_mobility_L.pdf')
    elif system_name == 'pa6x20':
        vals = profile.values()
        ave_dr = sum(v[0] for v in vals) / sum(v[1] for v in vals)
        print(ave_dr)
    else:
        raise NotImplementedError
