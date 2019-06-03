'''
Get some denstiy profiles for many systems:
    mmt
    pa6x20
    L
    S
    -mix,seg,MSd
'''


import sys

from md_work.structural_tools.datafile_content import DatafileContent
from utils.single_z_profile_composite_partial import single_z_profile_partial
from plotters.plot_data import plot_data


if __name__ == '__main__':
    multip = 10

    system_name = sys.argv[1]

    clay_idcs = set()
    mod_idcs = set()
    if system_name == 'L':
        data_dir = ('/media/anton/Seagate Expansion Drive/Backups/'
                    'Cluster calculations for article/BiggerSystems/Comp/'
                    '10chains/2.1 - Slow cooling (big)/1797434 - wiggle1')
        fnames = ['{0}/co.{1}.data'.format(data_dir, idx * 50000)
            for idx in range(1, 5)]
        for idx in range(9):
            clay_idcs.update(set(range(1 + idx*5380, 721 + idx*5380)))
            mod_idcs.update(set(range(721 + idx*5380, 1561 + idx*5380)))
        kwargs = {
            'clay_idcs': clay_idcs,
            'mod_idcs': mod_idcs,
            'replicate': True
        }
        try:
            kwargs['multip'] = multip
        except NameError:
            kwargs['multip'] = 1.0
        plot_kwargs = {
            'out_fname': 'dp_comp_L_partial.pdf',
            'x_limits': [15, 78],
        }
    elif system_name == 'S':
        data_dir = ('/media/anton/Seagate Expansion Drive/Backups/'
                    'Cluster calculations for article/BiggerSystems/Comp/5chains/'
                    '2.1 - Slow cooling (small)/1808725 - wiggle3')
        fnames = ['{0}/co.{1}.data'.format(data_dir, idx * 50000)
            for idx in range(1, 5)]
        for idx in range(9):
            clay_idcs.update(set(range(1 + idx*3470, 721 + idx*3470)))
            mod_idcs.update(set(range(721 + idx*3470, 1561 + idx*3470)))
        kwargs = {
            'clay_idcs': clay_idcs,
            'mod_idcs': mod_idcs,
            'replicate': True
        }
        try:
            kwargs['multip'] = multip
        except NameError:
            kwargs['multip'] = 1.0
        plot_kwargs = {
            'out_fname': 'dp_comp_S_partial.pdf',
            'x_limits': [-45, -5],
        }
    elif system_name == 'mix':
        data_dir = ('/media/anton/Seagate Expansion Drive/Backups/'
                    'Cluster calculations for article/1st conf/mixed ok/300K/'
                    'relaxation + wiggle/6th wiggle cycle (1432418)')
        fnames = ['{0}/co.{1}.data'.format(data_dir, idx * 50000)
            for idx in range(1, 5)]
        for idx in range(9):
            clay_idcs.update(set(range(1 + idx*3480, 721 + idx*3480)))
            mod_idcs.update(set(range(721 + idx*3480, 1561 + idx*3480)))

        kwargs = {
            'clay_idcs': clay_idcs,
            'mod_idcs': mod_idcs,
            'replicate': True
        }
        try:
            kwargs['multip'] = multip
        except NameError:
            kwargs['multip'] = 1.0
        plot_kwargs = {
            'out_fname': 'dp_comp_mix_partial.pdf',
            'x_limits': [4.5, 45],
        }
    elif system_name == 'seg':
        data_dir = ('/media/anton/Seagate Expansion Drive/Backups/'
                    'Cluster calculations for article/1st conf/segregated/300K/'
                    #'1st wiggle cycle (1414047)')
                    '2nd wiggle cycle (1426324)')
        fnames = ['{0}/co.{1}.data'.format(data_dir, idx * 50000)
            for idx in range(1, 5)]
        for idx in range(9):
            clay_idcs.update(set(range(1 + idx*3480, 721 + idx*3480)))
            mod_idcs.update(set(range(721 + idx*3480, 1561 + idx*3480)))
        kwargs = {
            'clay_idcs': clay_idcs,
            'mod_idcs': mod_idcs,
            'replicate': True
        }
        try:
            kwargs['multip'] = multip
        except NameError:
            kwargs['multip'] = 1.0
        plot_kwargs = {
            'out_fname': 'dp_comp_seg_partial.pdf',
            'x_limits': [10, 51],
        }
    else:
        print('Unknown system name:', sys.argv[1])
        sys.exit()

    kwargs['masses'] = {
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

    density_agg = {}
    for idx, fname in enumerate(fnames):
        print(idx + 1, '/', len(fnames))
        new_profile = single_z_profile_partial(DatafileContent(fname), **kwargs)
        for k, v in new_profile.items():
            try:
                density_agg[k]['clay'] += v['clay']
                density_agg[k]['modifier'] += v['modifier']
                density_agg[k]['polymer'] += v['polymer']
            except KeyError:
                try:
                    density_agg[k]['clay'] = v['clay']
                    density_agg[k]['modifier'] = v['modifier']
                    density_agg[k]['polymer'] = v['polymer']
                except KeyError:
                    density_agg[k] = {
                        'clay': v['clay'],
                        'modifier': v['modifier'],
                        'polymer': v['polymer']
                    }

    plot_kwargs['legends'] = ['MMT', 'Модификатор', 'Полимер', 'Все компоненты']
    plot_kwargs['y_limits'] = [0, 6]
    plot_kwargs['legend_location'] = 'upper right'

    xs = sorted(density_agg.keys())
    ys_clay = [density_agg[x]['clay']/len(fnames) for x in xs]
    ys_modifier = [density_agg[x]['modifier']/len(fnames) for x in xs]
    ys_polymer = [density_agg[x]['polymer']/len(fnames) for x in xs]
    ys_all = [(density_agg[x]['clay'] + density_agg[x]['modifier'] +
        density_agg[x]['polymer'])/len(fnames) for x in xs]

    xs = [x / multip - plot_kwargs['x_limits'][0] for x in xs]
    plot_kwargs['x_limits'] = [0,
        plot_kwargs['x_limits'][1] - plot_kwargs['x_limits'][0]] 

    plot_data(xs, yss=[ys_clay, ys_modifier, ys_polymer, ys_all],
        **plot_kwargs)
