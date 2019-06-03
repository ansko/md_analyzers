import math
import sys

from md_work.structural_tools.datafile_content import DatafileContent

from plotters.plot_data import plot_data


if __name__ == '__main__':
    characteristic = 'bonds'

    if len(sys.argv) < 3 or sys.argv[2] != 'angles':
        from utils.single_z_profile_bonds import single_z_profile
    else:
        from utils.single_z_profile_angles import single_z_profile
        characteristic = 'angles'

    print('Studying', characteristic)
    multip = 5

    system_name = sys.argv[1]
    if system_name == 'mmt':
        raise NotImplementedError()  # only ho-oh bonds that are not constant
    elif system_name == 'L':
        data_dir = ('/media/anton/Seagate Expansion Drive/Backups/'
                    'Cluster calculations for article/BiggerSystems/Comp/'
                    '10chains/2.1 - Slow cooling (big)/1797434 - wiggle1')
        fnames = ['{0}/co.{1}.data'.format(data_dir, idx * 50000)
            for idx in range(1, 51)]
        kwargs = {
            'ignored_bond_types': [1, 2, 6, 9, 10, 15],  # no h-bonds
            'ignored_angle_types': [4, 5, 6, 8, 11, 12, 13, 14, 15, 17, 18,
                23, 24, 27, 29]  # no h-bonds
        }
        try:
            kwargs['multip'] = multip
        except NameError:
            kwargs['multip'] = 1.0
        clay_idcs = set()
        for idx in range(9):
            clay_idcs.update(set(range(1 + idx*5380, 721 + idx*5380)))
        kwargs['clay_idcs'] = clay_idcs

        plot_kwargs = {
            'out_fname': characteristic + '_comp_L.pdf',
            'x_limits': [0, 27],
            'y_limits': [-0.3, 0.3]
        }
    elif system_name == 'S':
        data_dir = ('/media/anton/Seagate Expansion Drive/Backups/'
                    'Cluster calculations for article/BiggerSystems/Comp/5chains/'
                    '2.1 - Slow cooling (small)/1808725 - wiggle3')
        fnames = ['{0}/co.{1}.data'.format(data_dir, idx * 50000)
            for idx in range(1, 51)]

        kwargs = {
            'ignored_bond_types': [1, 2, 6, 9, 10, 15],  # no h-bonds
            'ignored_angle_types': [4, 5, 6, 8, 11, 12, 13, 14, 15, 17, 18,
                23, 24, 27, 29]  # no h-bonds
        }
        try:
            kwargs['multip'] = multip
        except NameError:
            kwargs['multip'] = 1.0
        clay_idcs = set()
        for idx in range(9):
            clay_idcs.update(set(range(1 + idx*3470, 721 + idx*3470)))
        kwargs['clay_idcs'] = clay_idcs

        plot_kwargs = {
            'out_fname': characteristic + '_comp_S.pdf',
            'x_limits': [0, 18],
            'y_limits': [-0.3, 0.3]
        }
    elif system_name == 'mix':
        data_dir = ('/media/anton/Seagate Expansion Drive/Backups/'
                    'Cluster calculations for article/1st conf/mixed ok/300K/'
                    'relaxation + wiggle/6th wiggle cycle (1432418)')
        fnames = ['{0}/co.{1}.data'.format(data_dir, idx * 50000)
            for idx in range(1, 51)]
        kwargs = {
            'ignored_bond_types': [1, 2, 6, 9, 10, 15],  # no h-bonds
            'ignored_angle_types': [4, 5, 6, 8, 11, 12, 13, 14, 15, 17, 18,
                23, 24, 27, 29]  # no h-bonds
        }
        try:
            kwargs['multip'] = multip
        except NameError:
            kwargs['multip'] = 1.0
        clay_idcs = set()
        for idx in range(9):
            clay_idcs.update(set(range(1 + idx*3480, 721 + idx*3480)))
        kwargs['clay_idcs'] = clay_idcs

        plot_kwargs = {
            'out_fname': characteristic + '_comp_mix.pdf',
            'x_limits': [0, 18],
            'y_limits': [-0.3, 0.3]
        }
    elif system_name == 'seg':
        data_dir = ('/media/anton/Seagate Expansion Drive/Backups/'
                    'Cluster calculations for article/1st conf/segregated/300K/'
                    '2nd wiggle cycle (1426324)')
        fnames = ['{0}/co.{1}.data'.format(data_dir, idx * 50000)
            for idx in range(1, 51)]
        kwargs = {
            'ignored_bond_types': [1, 2, 6, 9, 10, 15],  # no h-bonds
            'ignored_angle_types': [4, 5, 6, 8, 11, 12, 13, 14, 15, 17, 18,
                23, 24, 27, 29]  # no h-bonds
        }
        try:
            kwargs['multip'] = multip
        except NameError:
            kwargs['multip'] = 1.0
        clay_idcs = set()
        for idx in range(9):
            clay_idcs.update(set(range(1 + idx*3480, 721 + idx*3480)))
        kwargs['clay_idcs'] = clay_idcs

        plot_kwargs = {
            'out_fname': characteristic + '_comp_seg.pdf',
            'x_limits': [0, 16.5],
            'y_limits': [-0.3, 0.3]
        }
    elif system_name == 'pa6x20':
        data_dir = ('/media/anton/Seagate Expansion Drive/Backups/'
                    'Cluster calculations for article/BiggerSystems/PA6x20/'
                    '1993784')
        fnames = ['{0}/co.{1}.data'.format(data_dir, idx * 50000)
            for idx in range(1, 51)]
        kwargs = {
        }
        try:
            kwargs['multip'] = multip
        except NameError:
            kwargs['multip'] = 1.0
        plot_kwargs = {
            'out_fname': characteristic + '_comp_pa6x20.pdf',
            'x_limits': [0, 45],
            'y_limits': [-0.3, 0.3]
        }
    else:
        print('Unknown system name:', sys.argv[1])
        sys.exit()


    orientation_agg = {}
    for idx, fname in enumerate(fnames):
        print(idx + 1, '/', len(fnames))
        current_orientation = single_z_profile(DatafileContent(fname), **kwargs)
        for k, v in current_orientation.items():
            try:
                orientation_agg[k][0] += v[0]
                orientation_agg[k][1] += v[1]
            except KeyError:
                orientation_agg[k] = [v[0], v[1]]

    xs = sorted(orientation_agg.keys())
    ys = [orientation_agg[x][0] / orientation_agg[x][1] for x in xs]
    xs = [x / multip for x in xs]

    f = open('logs/bonds_orientation_{}'.format(system_name), 'w')
    for idx in range(len(xs)):
        print(xs[idx], ys[idx], file=f)
    
    plot_data(xs, ys=ys, **plot_kwargs)
