'''
Gyration radii studies
'''

from collections import Counter
import os
import sys

from md_work.structural_tools.datafile_content import DatafileContent


if __name__ == '__main__':
    multip = 10

    system_name = sys.argv[1]
    if system_name == 'mmt':
        raise NotImplementedError()
    elif system_name == 'L':
        multip = 0.2
        from plotters.plot_data import plot_data
        from utils.gyration_composite import gyration_composite
        from utils.single_z_profile_composite import single_z_profile
        data_dir = ('/media/anton/Seagate Expansion Drive/Backups/'
                    'Cluster calculations for article/BiggerSystems/Comp/'
                    '10chains/2.2 - More relaxation 500 (wiggle)/'
                    '1785842 - relaxation')
        fnames = ['{0}/co.{1}.data'.format(data_dir, idx * 50000)
            for idx in range(1, 51)]
        clay_idcs = set()
        for idx in range(9):
            clay_idcs.update(set(range(1 + idx*5380, 721 + idx*5380)))
        def poly_mol_idx_seg(atom_id):
            mini_system_idx = (atom_id - 1) // 5380
            in_mini_system_id = atom_id - mini_system_idx * 5380
            if in_mini_system_id < 1561:
                return None
            in_mini_system_id -= 1560
            in_mini_system_idx = (in_mini_system_id - 1) // 382
            return 10 * mini_system_idx + in_mini_system_idx
        masses = {
            1:  26.981540, # ao
            2:  24.305000, # mgo
            3:  28.085500, # st
            4:  15.999400, # ob
            5:  15.999400, # oh
            6:  15.999400, # obts
            7:  15.999400, # ohs
            8:   1.007970, # ho
            9:  14.006700, # n4
           10:  12.011150, # c2
           11:  12.011150, # c3
           12:   1.007970, # h
           13:   1.007970, # hn
           14:  12.011150, # c'
           15:  15.999400, # o'
           16:  14.006700, # n2
           17:  14.006700, # n
        }
        kwargs = {
            'poly_mol_idx': poly_mol_idx_seg,
            'multip': multip,
            'clay_idcs': clay_idcs,
            'masses': masses
        }
        plot_kwargs = {
            'out_fname': 'gyrs_L.eps',
            'legends': [r'$R_{xx}$', r'$R_{yy}$', r'$R_{zz}$'],
            'y_limits': [0, 12]
        }
        radii = []
        for idx, fname in enumerate(fnames):
            print(idx + 1, '/', len(fnames))
            radii.append(gyration_composite(DatafileContent(fname), **kwargs))
            print(min(radii[-1].keys()), max(radii[-1].keys()))
        all_keys = set()
        for r in radii:
            all_keys.update(r.keys())
        all_keys = list(sorted(all_keys))
        ave_radii = {k: {'rx': 0, 'ry': 0, 'rz': 0, 'count': 0}
            for k in all_keys}
        for r in radii:
            for k, v in r.items():
                ave_radii[k]['rx'] += v['rx'] * v['count']
                ave_radii[k]['ry'] += v['ry'] * v['count']
                ave_radii[k]['rz'] += v['rz'] * v['count']
                ave_radii[k]['count'] += v['count']
        xs = [x for x in sorted(ave_radii.keys()) if ave_radii[k]['count']]
        ys_rx = [ave_radii[x]['rx'] / ave_radii[x]['count'] for x in xs]
        ys_ry = [ave_radii[x]['ry'] / ave_radii[x]['count'] for x in xs]
        ys_rz = [ave_radii[x]['rz'] / ave_radii[x]['count'] for x in xs]
        xs = [x / multip for x in xs]
        plot_data(xs=xs, yss=[ys_rx, ys_ry, ys_rz], **plot_kwargs)
        with open('logs/gyr_L', 'w') as f:
            for idx in range(len(xs)):
                print(xs[idx], ys_rx[idx], ys_ry[idx], ys_rz[idx], file=f)
    elif system_name == 'S':  # Too small data
        raise NotImplementedError()
    elif system_name == 'mix':
        raise NotImplementedError()
    elif system_name == 'seg':
        raise NotImplementedError()
    elif system_name == 'pa6x20':
        # 144 chains x 382 atoms
        from utils.gyration_polymer import gyration_polymer
        data_dir = ('/media/anton/Seagate Expansion Drive/'
                    'md_new_data/pa6x20/3507904_500K/datafiles')
        fnames = ['{0}/pa6x20.{1}.data'.format(data_dir, idx * 50000)
            for idx in range(1, 11)]
        masses = {
            1: 1.00797,
            2: 12.0112,
            3: 1.00797,
            4: 12.0112,
            5: 15.9994,
            6: 14.0067,
            7: 14.0067,
            8: 12.0112
        }
        kwargs = {
            'poly_mol_idx': lambda atom_id: (atom_id - 1) // 382,
            'masses': masses
        }
        # Get only average components
        radii = []
        for idx, fname in enumerate(fnames):
            print(idx + 1, '/', len(fnames))
            radii.append(gyration_polymer(DatafileContent(fname), **kwargs))
        ave_rx = sum(r['x'] for r in radii) / len(radii)
        ave_ry = sum(r['y'] for r in radii) / len(radii)
        ave_rz = sum(r['z'] for r in radii) / len(radii)
        print('Polymer:')
        print('\trx =', ave_rx)
        print('\try =', ave_ry)
        print('\trz =', ave_rz)
    else:
        print('Unknown system name:', sys.argv[1])
        sys.exit()
