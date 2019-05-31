'''
Get some denstiy profiles for many systems.
'''


import sys
from collections import Counter

from md_work.structural_tools.datafile_content import DatafileContent

from plotters.plot_data import plot_data


if __name__ == '__main__':
    multip = 10

    system_name = sys.argv[1]
    if system_name == 'mmt':
        from utils.single_z_profile_clay import single_z_profile
        data_dir = ('/media/anton/Seagate Expansion Drive/Backups/'
                    'Cluster calculations for article/1st conf/mmt/'
                    'relaxation (1438006)')
        fnames = ['{0}/mmt.{1}.data'.format(data_dir, idx * 50000)
            for idx in range(1, 11)]
        masses = {
            1: 26.981540, # ao
            2: 24.305000, # mgo
            3: 28.085500, # st
            4: 15.999400, # ob
            5: 15.999400, # oh
            6: 15.999400, # obts
            7: 15.999400, # ohs
            8:  1.007970, # ho
            9: 22.990000  # Na
        }
        kwargs = {
            'masses': masses,
            'replicate': True
        }
        try:
            kwargs['multip'] = multip
        except NameError:
            kwargs['multip'] = 1.0
        plot_kwargs = {
            'out_fname': 'dp_mmt.png',
            'x_limits': [-1, 17.5]
        }
    elif system_name == 'L':
        from utils.single_z_profile_composite import single_z_profile
        data_dir = ('/media/anton/Seagate Expansion Drive/Backups/'
                    'Cluster calculations for article/BiggerSystems/Comp/'
                    '10chains/2.1 - Slow cooling (big)/1797434 - wiggle1')
        fnames = ['{0}/co.{1}.data'.format(data_dir, idx * 50000)
            for idx in range(1, 5)]
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
            'masses': masses,
        }
        try:
            kwargs['multip'] = multip
        except NameError:
            kwargs['multip'] = 1.0
        plot_kwargs = {
            'out_fname': 'dp_comp_L.png',
            'x_limits': [0, 50],
            'y_limits': [0, 2]
        }
    else:
        print('Unknown system name:', sys.argv[1])
        sys.exit()

    density_agg = Counter()
    for idx, fname in enumerate(fnames):
        print(idx + 1, '/', len(fnames))
        density_agg += single_z_profile(DatafileContent(fname), **kwargs)

    xs = sorted(density_agg.keys())
    ys = [density_agg[x] for x in xs]
    xs = [x / multip for x in xs]
    
    plot_data(xs, ys=ys, **plot_kwargs)
