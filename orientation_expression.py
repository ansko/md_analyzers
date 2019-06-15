import math
import os
import sys

from md_work.structural_tools.datafile_content import DatafileContent


if __name__ == '__main__':
    characteristic = 'bonds'

    if len(sys.argv) < 2:
        print('Specify property!')
        sys.exit()
    if sys.argv[1] == 'bonds':
        from utils.single_z_profile_bonds import single_z_profile
    elif sys.argv[1] == 'angles':
        from utils.single_z_profile_angles import single_z_profile
        characteristic = 'angles'

    print('Studying', characteristic)
    multip = 5

    data_dir = ('/media/anton/Seagate Expansion Drive/md_new_data/comp_L/300K/')
    subdirs = ['3514713', '3521870', '3529081', '3545491']

    for subdir_idx in range(len(subdirs)):
        full_path = data_dir + subdirs[subdir_idx] + '/datafiles/'
        fname_part = os.listdir(full_path)[0].split('.')[0]
        fnames = ['{0}/{1}.{2}.data'.format(full_path, fname_part, idx * 50000)
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

        orientation_agg = {}
        for idx, fname in enumerate(fnames):
            print(idx + 1, '/', len(fnames))
            current_orientation = single_z_profile(DatafileContent(fname),
                clay_idcs=clay_idcs)
            for k, v in current_orientation.items():
                try:
                    orientation_agg[k][0] += v[0]
                    orientation_agg[k][1] += v[1]
                except KeyError:
                    orientation_agg[k] = [v[0], v[1]]

        xs = sorted(orientation_agg.keys())
        ys = [orientation_agg[x][0] / orientation_agg[x][1] for x in xs]
        xs = [x / multip for x in xs]

        exspressiveness_abs_sum = 0
        exspressiveness_sq_sum = 0
        for y in ys:
            exspressiveness_abs_sum += abs(y)
            exspressiveness_sq_sum += y**2
        exspressiveness_abs_sum /= len(xs)
        exspressiveness_sq_sum /= len(xs)

        print(subdirs[subdir_idx], len(fnames))
        print('exspressiveness_abs_sum', exspressiveness_abs_sum)
        print('exspressiveness_sq_sum', exspressiveness_sq_sum)
