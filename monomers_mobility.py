from collections import Counter
import math
import os
import sys

from md_work.structural_tools.datafile_content import DatafileContent


if __name__ == '__main__':
    multip = 10

    system_name = sys.argv[1]

    if system_name == 'mmt':
        raise NotImplementedError()
    elif system_name == 'L':
        from utils.monomers_cms import monomers_cms
        from plotters.plot_data import plot_data
        #data_dir = ('/media/anton/Seagate Expansion Drive/Backups/'
        #            'Cluster calculations for article/BiggerSystems/Comp/'
        #            '10chains/2.2 - More relaxation 500 (wiggle)/'
        #            '1785842 - relaxation')
        data_dir = ('/media/anton/Seagate Expansion Drive/md_new_data/comp_L/'
                    '300K/3514713/datafiles')
        fnames = ['{0}/comp_L_d_53.{1}.data'.format(data_dir, idx * 50000)
            for idx in range(1, 51)]
        def monomer_idx(atom_id):
            # 0 -- 1800
            mini_system_idx = (atom_id - 1) // 5380
            in_mini_system_id = atom_id - mini_system_idx * 5380
            if in_mini_system_id < 1561:
                return None
            in_mini_system_id -= 1560
            in_mini_system_idx = (in_mini_system_id - 1) // 382
            in_chain_id = in_mini_system_id - in_mini_system_idx*382
            in_chain_idx = None
            if in_chain_id < 22:
                in_chain_idx = 0
            elif in_chain_id > 361:
                in_chain_idx = 19
            else:
                in_chain_idx = 1 + (in_chain_id - 22) // 19
            return 200 * mini_system_idx + in_mini_system_idx * 20 + in_chain_idx
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
        plot_kwargs = {
            'out_fname': 'mon_mob_L.eps',
        }
        plot_log_kwargs = {
            'out_fname': 'mon_mob_L_log.eps',
        }

        starting = monomers_cms(DatafileContent(fnames[0]), masses, monomer_idx)
        f = open('logs/monomers_mobility_L', 'w')
        xs = []
        ys = []
        for idx, fname in enumerate(fnames):
            dfc = DatafileContent(fname)
            lx = dfc.xhi - dfc.xlo
            ly = dfc.yhi - dfc.ylo
            lz = dfc.zhi - dfc.zlo
            time = idx * 50  # ps
            current = monomers_cms(dfc,masses,  monomer_idx)
            delta = 0
            for k in current.keys():
                dx = abs(current[k]['x_cm'] - starting[k]['x_cm'])
                dy = abs(current[k]['y_cm'] - starting[k]['y_cm'])
                dz = abs(current[k]['z_cm'] - starting[k]['z_cm'])
                dx = min(dx, lx - dx)
                dy = min(dy, ly - dy)
                dz = min(dz, lz - dz)
                delta += (dx**2 + dy**2 + dz**2)**0.5
            delta /= len(current)
            xs.append(time)
            ys.append(delta)
            print(idx + 1, '/', delta, len(current))
            print(time, delta, file=f)
        log_xs = [math.log(x) for x in xs[1:]]
        log_ys = [math.log(y) for y in ys[1:]]
        plot_data(xs, ys=ys, **plot_kwargs)
        plot_data(log_xs, ys=log_ys, **plot_log_kwargs)
    elif system_name == 'S':  # Too small data
        raise NotImplementedError()
    elif system_name == 'mix':
        raise NotImplementedError()
    elif system_name == 'seg':
        raise NotImplementedError()
    elif system_name == 'pa6x20':
        from utils.monomers_cms import monomers_cms
        from plotters.plot_data import plot_data
        from utils.gyration_polymer import gyration_polymer
        data_dir = ('/media/anton/Seagate Expansion Drive/'
                    'md_new_data/pa6x20/3507904_500K/datafiles')
        fnames = ['{0}/pa6x20.{1}.data'.format(data_dir, idx * 50000)
            for idx in range(1, 51)]

        def monomer_idx(atom_id):
            # 0 -- 2879
            chain_id = atom_id // 382
            in_chain_id = atom_id - chain_id * 382
            if in_chain_id < 22:
                in_chain_idx = 0
            elif in_chain_id > 361:
                in_chain_idx = 19
            else:
                in_chain_idx = 1 + (in_chain_id - 22) // 19
            return chain_id * 20 + in_chain_idx
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
        plot_kwargs = {
            'out_fname': 'mon_mob_pa6x20.eps',
        }
        plot_log_kwargs = {
            'out_fname': 'mon_mob_pa6x20_log.eps',
        }
        starting = monomers_cms(DatafileContent(fnames[0]), masses, monomer_idx)
        f = open('logs/monomers_mobility_pa6x20', 'w')
        xs = []
        ys = []
        for idx, fname in enumerate(fnames):
            dfc = DatafileContent(fname)
            lx = dfc.xhi - dfc.xlo
            ly = dfc.yhi - dfc.ylo
            lz = dfc.zhi - dfc.zlo
            time = idx * 50  # ps
            current = monomers_cms(dfc,masses,  monomer_idx)
            delta = 0
            for k in current.keys():
                dx = abs(current[k]['x_cm'] - starting[k]['x_cm'])
                dy = abs(current[k]['y_cm'] - starting[k]['y_cm'])
                dz = abs(current[k]['z_cm'] - starting[k]['z_cm'])
                dx = min(dx, lx - dx)
                dy = min(dy, ly - dy)
                dz = min(dz, lz - dz)
                delta += (dx**2 + dy**2 + dz**2)**0.5
            delta /= len(current)
            xs.append(time)
            ys.append(delta)
            print(idx + 1, '/', delta, len(current))
            print(time, delta, file=f)
        log_xs = [math.log(x) for x in xs[1:]]
        log_ys = [math.log(y) for y in ys[1:]]
        plot_data(xs, ys=ys, **plot_kwargs)
        plot_data(log_xs, ys=log_ys, **plot_log_kwargs)
    else:
        print('Unknown system name:', sys.argv[1])
        sys.exit()
