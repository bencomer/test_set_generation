from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen import Structure, Lattice, MPRester, Molecule
from ase.data import chemical_symbols, reference_states
from ase.build import bulk
from ase.units import Bohr
from ase.visualize import view
from pickle import dump, load
import numpy as np
import json
from collections import OrderedDict
import math
import os


def select_kpts(structure, kpt_density=15.):
    kpts = []
    atoms = AseAtomsAdaptor.get_atoms(structure)
    for i in range(3):
        l = np.linalg.norm(atoms.cell[i])
        kpts.append(math.ceil(kpt_density / l))
    return kpts

def old_find_shift_abinit(kpts):
    shift = ' '
    for kpt in kpts:
        if kpt % 2:
            shift += '0 '
        else:
            shift += '0.5 '
    return shift

def find_shift_abinit(kpts):
    shift = ' 0 0 0'
    return shift

num_GGA = 0

def add_to_dict(structure):
    formula = structure.composition.reduced_formula
    structures[formula+'_bulk'] = structure.as_dict()
    kpts = select_kpts(structure)
    sparc_p = sprc_params.copy()
    esp_p = esp_params.copy()
    sparc_p['KPOINT_GRID'] = kpts
    esp_p['kpts'] = kpts
    abinit_p = abinit_params.copy()
    abinit_p['kptopt'] = 1
    abinit_p['ngkpt'] = ' '.join([str(a) for a in kpts])
    #abinit_p['nshiftk'] = 1
    abinit_p['shiftk'] = find_shift_abinit(kpts)
    syms =structure.symbol_set
    contains_non_metals = 'O' in syms or 'H' in syms or 'S' in syms or 'N' in syms
    if ('O' in syms or formula == 'PSe') and formula != 'CaO':
        sparc_p['EXCHANGE_CORRELATION'] = 'GGA_PBE'
        abinit_p['xc'] = 'PBE'
        esp_p['xc'] = 'PBE'
    if formula == 'MgCl2':
        abinit_p['diemix'] = 0.1
        esp_p['convergence']['mixing'] = 0.1
        sparc_p['MIXING_PARAMETER'] = 0.1


    for site in structure:
        if hasattr(site, 'magmom'):
            if 'Co' not in formula:
                continue
            abinit_p['nsppol'] = 2
            sparc_p['SPIN_TYP'] = 2
            esp_p['spinpol'] = True


    params_dict[formula+'_bulk'] = {}
    params_dict[formula+'_bulk']['sparc'] = sparc_p.copy()
    params_dict[formula+'_bulk']['espresso'] = esp_p.copy()
    params_dict[formula+'_bulk']['abinit'] = abinit_p.copy()



sprc_params = {}
sprc_params['h'] = 0.4 * Bohr
sprc_params['MAXIT_SCF'] = 250
#sprc_params['TOL_SCF_QE'] = 1e-6
sprc_params['TOL_SCF'] = 5e-5
sprc_params['CALC_STRESS'] = 1
sprc_params['BC'] = [True, True, True]
sprc_params['EXCHANGE_CORRELATION'] = 'LDA_PZ'
sprc_params['MIXING_VARIABLE'] = 'density'
sprc_params['MIXING_PRECOND'] = 'kerker'
sprc_params['ELEC_TEMP'] = 1160.45
sprc_params['KPOINT_SHIFT'] = ' 0 0 0'
#params['MIXING_VARIABLE'] = 'potential'


esp_params = dict(
        #xc='PZ',
        kpts=(1, 1, 1), #only need 1 kpt in z-direction
        pw=300,
        dw=3000,
        calculation='scf',
        #spinpol=True,
        tstress=True,
        convergence={'energy':1e-6,
                    #'mixing':0.1,
                    'maxsteps':1000,
                    'diag':'david'},
        nosym=True,
        startingwfc='random',
        smearing='gaussian', #fermi-dirac electron smearing
        sigma=0.1, #smearing width
        dipole={'status':False}, #dipole corrections True turns them on
        outdir ='.'
        )

abinit_params = {'xc':'LDA_PZ',
         'autoparal':4,
         'nsym':1,
         'optcell':0,
         'optstress':1,
         'nstep':500,
         #'pps':'oncv',
         'pps':'pot',
         'ecut': 2721.14,
         'tsmear':0.1,
         'chksymbreak':0,
         'toldfe':1.0e-6,
         'occopt':7,
}


rest = MPRester('XaCJrv4nBIeuy3kd')

params_dict = {}
structures = {}

psps = os.listdir('/home/bcomer3/sparc/ase_sparc/pysparcx/sparc/pseudos/LDA_pseudos/')
available_psps = [a.split('.')[0] for a in psps]



# metals
metals = ['Pt', 'Sc', 'La', 'Na', 'Au', 'Os', 'V', 'Zn', 'Bi', 'Co', 'Pd']
sizes = [(4,2,2), (4,3,2), (3,1,1), (3,2,1), (1,1,7), (2,2,2), (2,1,1),(2,2,1), (1,1,1), (2,1,1), (2,2,1)]

for size, metal in zip(sizes, metals):
    struc_l = rest.get_entries(metal, sort_by_e_above_hull=True)
    struc_l = struc_l[0]
    struct = rest.get_structure_by_material_id(struc_l.entry_id) 
    if metal == 'Pt':
        struct.add_site_property('magmom', [0] * len(struct))
    struct *= size
    if metal in ['Pt', 'Sc']:
        del struct[0]
    if metal == 'Hf':
        del struct[50]
        del struct[35]
        del struct[100]
        #print(len(struct))
    if metal == 'Sc':
        print(len(struct))
    add_to_dict(struct)

# oxides
oxide_formulas = ['CaO', 'ZrO2', 'Pb3O4',
                  'Al2O3', 'Ag2O', 'SnO']
oxide_mps = ['mp-2605', 'mp-2858', 'mp-22633', 
             'mp-353', 'mp-2097']

oxide_sizes = [(2,1,2), (1,1,1), (1,2,1), (1,1,1), (1,1,1), (2,2,2), (1,1,1),
               (1,1,1), (1,1,1), (2,1,1)]

for formula, mp_id, size in zip(oxide_formulas, oxide_mps, oxide_sizes):
    struct = rest.get_structure_by_material_id(mp_id)
    struct *= size
    add_to_dict(struct)


# Alloys
alloy_formulas = ['SiGe', 'Al3Sn', 'BeCu', 'PtAu', 'Zn3Cu', 'IrPt3', 'Mo4Ru', 'SbAs']
alloy_mps = ['mp-1094056','mp-1183200', 'mp-2323', 'mp-1219709',
             'mp-972042', 'mp-1184759', 'mp-1221461', 'mp-1219457']

for formula, mp_id in zip(alloy_formulas, alloy_mps):
    struct = rest.get_structure_by_material_id(mp_id)
    atoms = AseAtomsAdaptor.get_atoms(struct)
    add_to_dict(struct)

# Sulfides
sulfide_formulas = ['AsS', 'CdS', 'NaS', 'Sb2S3', 'RhS2']
sulfide_mps = ['mp-542810', 'mp-672', 'mp-648', 'mp-2809', 'mp-22555']

for formula, mp_id in zip(sulfide_formulas, sulfide_mps):
    struct = rest.get_structure_by_material_id(mp_id)
    add_to_dict(struct)

# other
other_formulas = ['MgS', 'CaN2', 'WC', 'BeF2', 'PSe', 'NbP']
other_mps = ['mp-1315', 'mp-1009657', 'mp-1894', 'mp-15951', 'mp-28885', 'mp-9339']

for formula, mp_id in zip(other_formulas, other_mps):
    struct = rest.get_structure_by_material_id(mp_id)
    atoms = AseAtomsAdaptor.get_atoms(struct)
    if formula == 'MgS':
        struct *= (3,3,3)
    add_to_dict(struct)


lens = []
for n, s in structures.items():
    lens.append(len(Structure.from_dict(s)))
    #print(len(Structure.from_dict(s)))

#print(([a > 30 for a in lens]).count(True))
#print(([a > 100 for a in lens]).count(True))

print(len(structures))

json.dump(structures, open('../../test_set/bulk.json', 'w'))
json.dump(params_dict, open('../../test_set/bulk_parameters.json', 'w'))
