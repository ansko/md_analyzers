from collections import Counter
import math
import os
import sys

from md_work.structural_tools.datafile_content import DatafileContent
from utils.atomic_mobility_profile import atomic_mobility_profile

if __name__ == '__main__':
    multip = 0.5

    system_name = sys.argv[1]

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
        kwargs = {
            'clay_idcs': clay_idcs,
            'multip': multip
        }
        plot_kwargs = {
            'out_fname': 'atoms_mob_L.eps',
        }
        profile_agg = {}
        f = open('logs/atomic_mobility_L', 'w')
        dfc_last = DatafileContent(fnames[0])
        dfc_next = DatafileContent(fnames[-1])
        profile = atomic_mobility_profile(dfc_last, dfc_next, **kwargs)
        xs = sorted(profile.keys())
        ys = [profile[x][0] / profile[x][1] for x in xs]
        plot_data(xs, ys=ys, **plot_kwargs)
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
        kwargs = {
            'clay_idcs': set(),
            'multip': multip
        }
        plot_kwargs = {
            'out_fname': 'atoms_mob_pa6x20.eps',
        }
        f = open('logs/atomic_mobility_pa6x20', 'w')
        dfc_last = DatafileContent(fnames[0])
        dfc_next = DatafileContent(fnames[-1])
        profile = atomic_mobility_profile(dfc_last, dfc_next, **kwargs)
        xs = sorted(profile.keys())
        ys = [profile[x][0] / profile[x][1] for x in xs]
        plot_data(xs, ys=ys, **plot_kwargs)
    else:
        print('Unknown system name:', sys.argv[1])
        sys.exit()
