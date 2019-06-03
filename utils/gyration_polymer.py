from collections import Counter

from utils.single_molecule_gyration_tensor import single_molecule_gyration_tensor


def gyration_polymer(dfc, **kwargs):
    result = {'x': 0, 'y': 0, 'z': 0}

    # Distribute atoms into their molecules
    molecules = {}
    for atom in dfc.atoms:
        k = kwargs['poly_mol_idx'](atom['atom_id'])
        try:
            molecules[k].append(atom)
        except KeyError:
            molecules[k] = [atom]

    for k, v in molecules.items():
        tmp_result = single_molecule_gyration_tensor(dfc, v, kwargs['masses'])
        result['x'] += tmp_result['x']
        result['y'] += tmp_result['y']
        result['z'] += tmp_result['z']

    result['x'] /= len(molecules)
    result['y'] /= len(molecules)
    result['z'] /= len(molecules)

    return result
