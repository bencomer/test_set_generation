from pymatgen.io.ase import AseAtomsAdaptor
from ase.units import Bohr
from pymatgen.core.surface import generate_all_slabs
from pymatgen import Structure, Lattice, MPRester, Molecule
from ase.data import chemical_symbols, reference_states
from ase.atoms import Atoms
from ase.build import bulk, surface, molecule, add_adsorbate
from ase.build import graphene_nanoribbon
from ase.build.surfaces_with_termination import surfaces_with_termination
from ase.visualize import view
from pickle import dump, load
import numpy as np
import json
from collections import OrderedDict
import math
import os


def select_kpts(structure, kpt_density=20.):
    kpts = []
    atoms = AseAtomsAdaptor.get_atoms(structure)
    for i in range(3):
        if i == 2:
            kpts.append(1)
        else:
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


def add_to_dict(structure, element, miller_index):
    # formula = structure.composition.reduced_formula
    miller_index_str = ''
    for j in miller_index:
        miller_index_str +=  str(j)
        
    slab = structure.as_dict()
    if 'scale_factor' in slab:
        slab['scale_factor'] = slab['scale_factor'].tolist()

    structures[element + '_' + miller_index_str + '_surface'] = slab
    kpts = select_kpts(structure)
    sparc_p = sprc_params.copy()
    esp_p = esp_params.copy()
    abinit_p = abinit_params.copy()
    syms =structure.symbol_set
    if 'O' in syms or 'S' in syms or 'P' in syms or 'C' in syms or 'Ga' in syms:
        if 'O' in syms:
            print('yeet')
        #if element not in ['Ru_NOads', 'Cu_COads', 'TiO2_defect']:
        sparc_p['EXCHANGE_CORRELATION'] = 'GGA_PBE'
        abinit_p['xc'] = 'PBE'
        esp_p['xc'] = 'PBE'
    sparc_p['KPOINT_GRID'] = kpts
    esp_p['kpts'] = kpts
    abinit_p['kptopt'] = 1
    abinit_p['ngkpt'] = ' '.join([str(a) for a in kpts])
    abinit_p['nshiftk'] = 1
    abinit_p['shiftk'] = find_shift_abinit(kpts)
    sparc_p['MIXING_PARAMETER'] = 0.1
    abinit_p['diemix'] = 0.1
    esp_p['convergence']['mixing'] = 0.1




    params_dict[element + '_' + miller_index_str + '_surface'] = {}
    params_dict[element + '_' +  miller_index_str + '_surface']['sparc'] = sparc_p.copy()
    params_dict[element + '_' + miller_index_str + '_surface']['espresso'] = esp_p.copy()
    params_dict[element + '_' + miller_index_str + '_surface']['abinit'] = abinit_p.copy()
        


sprc_params = {}
sprc_params['h'] = 0.4 * Bohr
sprc_params['MAXIT_SCF'] = 1000
#sprc_params['TOL_SCF_QE'] = 1e-6
sprc_params['TOL_SCF'] = 5e-5
sprc_params['BC'] = [True] * 3
sprc_params['EXCHANGE_CORRELATION'] = 'LDA_PZ'
sprc_params['MIXING_VARIABLE'] = 'density'
sprc_params['MIXING_PRECOND'] = 'kerker'
sprc_params['ELEC_TEMP'] = 1160.45
sprc_params['KPOINT_SHIFT'] = ' 0 0 0'
#params['MIXING_VARIABLE'] = 'potential'


esp_params = dict(#xc='PZ',
        kpts=(1, 1, 1), #only need 1 kpt in z-direction
        pw=300,
        dw=3000,
        calculation='scf',
        #spinpol=True,
        convergence={'energy':1e-6,
                    'mixing':0.1,
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
         'optcell':0,
         'nsym':1,
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
metals = ['Al', 'Ru', 'Cu', 'Bi']

for metal in metals:
    struc_l = rest.get_entries(metal, sort_by_e_above_hull=True)
    struc_l = struc_l[0]
    struct = rest.get_structure_by_material_id(struc_l.entry_id)
    atoms = AseAtomsAdaptor.get_atoms(struct) # switch to ASE
    for millers in [(1,1,1)]:
        slab = surface(atoms, millers, 4, vacuum=6)
        slab *= (2,2,1)
        n_metal = metal
        if millers == (1,1,1):
            if metal == 'Cu':
                add_adsorbate(slab, molecule('CO'), 2.5, (0, 1.64))
                n_metal += '_COads'

        structure = AseAtomsAdaptor.get_structure(slab)
        add_to_dict(structure, n_metal, millers)



metals = ['In', 'Rb', 'Sr', 'Te']

for metal in metals:
    struc_l = rest.get_entries(metal, sort_by_e_above_hull=True)
    struc_l = struc_l[0]
    struct = rest.get_structure_by_material_id(struc_l.entry_id)
    atoms = AseAtomsAdaptor.get_atoms(struct) # switch to ASE
    for millers in [(1,0,0)]:
        if metal not in ['Rb']:
            slab = surface(atoms, millers, 2, vacuum=6)
            slab *= (2,2,1)
        else:
            slab = surface(atoms, millers, 1, vacuum=6)
        n_metal = metal
        if millers == (1,1,1):
            if metal == 'Cu':
                add_adsorbate(slab, molecule('CO'), 2.5, (0, 1.64))
                n_metal += '_COads'
        if millers == (1,0,0):
            if metal == 'Fe':
                add_adsorbate(slab, molecule('O2'), 2.8, (0.411, 2.3255))
                n_metal += '_O2ads'
            if metal == 'Ru':
                add_adsorbate(slab, molecule('NO'), 1.5, (2.739, 3.219))
                n_metal += '_NOads'
            if metal == 'Cu':
                slab *= (2,2,1)
                del slab[10]
                print(len(slab))
                n_metal += '_defect'
        structure = AseAtomsAdaptor.get_structure(slab)
        add_to_dict(structure, n_metal, millers)

metals = ['Tc', 'Ga', 'Ba', 'Cs']

for metal in metals:
    struc_l = rest.get_entries(metal, sort_by_e_above_hull=True)
    struc_l = struc_l[0]
    struct = rest.get_structure_by_material_id(struc_l.entry_id)
    atoms = AseAtomsAdaptor.get_atoms(struct) # switch to ASE
    for millers in [(1,1,0)]:
        if metal not in ['Cs', 'Ga']:
            slab = surface(atoms, millers, 4, vacuum=6)
            slab *= (2,2,1)
        elif metal == 'Ga':
            slab = surface(atoms, millers, 2, vacuum=6)
            slab *= (2,2,1)
        else:
            slab = surface(atoms, millers, 1, vacuum=6)
        n_metal = metal
        if millers == (1,1,1):
            if metal == 'Cu':
                add_adsorbate(slab, molecule('CO'), 2.5, (0, 1.64))
                n_metal += '_COads'
        if millers == (1,0,0):
            if metal == 'Fe':
                add_adsorbate(slab, molecule('O2'), 2.8, (0.411, 2.3255))
                n_metal += '_O2ads'
            if metal == 'Ru':
                add_adsorbate(slab, molecule('NO'), 1.5, (2.739, 3.219))
                n_metal += '_NOads'
            if metal == 'Cu':
                slab *= (2,2,1)
                del slab[10]
                print(len(slab))
                n_metal += '_defect'
        structure = AseAtomsAdaptor.get_structure(slab)
        add_to_dict(structure, n_metal, millers)

# Si

struc_l = rest.get_entries('Si', sort_by_e_above_hull=True)
struc_l = struc_l[0]
struct = rest.get_structure_by_material_id(struc_l.entry_id)

atoms = AseAtomsAdaptor.get_atoms(struct) # switch to ASE

for millers in [(1,1,1)]:
    slab = surfaces_with_termination(atoms, millers, 3, vacuum=6)[0]
    structure = AseAtomsAdaptor.get_structure(slab)
    add_to_dict(structure, 'Si', millers)

# Zn

struc_l = rest.get_entries('Zn', sort_by_e_above_hull=True)
struc_l = struc_l[0]
struct = rest.get_structure_by_material_id(struc_l.entry_id)

atoms = AseAtomsAdaptor.get_atoms(struct) # switch to ASE

for millers in [(1,0,0)]:
    slab = surfaces_with_termination(atoms, millers, 3, vacuum=6)[0]
    structure = AseAtomsAdaptor.get_structure(slab)
    add_to_dict(structure, 'Zn', millers)



# Carbon Bulk Materials

# diamond
material = ['mp-66']
name = ['diamond']

bulk = rest.get_structure_by_material_id('mp-66')
bulk = AseAtomsAdaptor.get_atoms(bulk)
for millers in [(1,1,1)]:
    slab = surfaces_with_termination(bulk, millers, 4, vacuum=6)[0]
    structure = AseAtomsAdaptor.get_structure(slab)
    add_to_dict(structure, 'diamond', millers)

# graphite


material = ['mp-48']
name = ['graphite']

bulk = rest.get_structure_by_material_id('mp-48')
bulk = AseAtomsAdaptor.get_atoms(bulk)
# correct spacing between layers
new_cell = bulk.get_cell()  * np.array([1,1,0.85853])
bulk.set_cell(new_cell,
              scale_atoms=True)
for millers in [(0,0,1)]:
    slab = surfaces_with_termination(bulk, millers, 2, vacuum=6)[0]
    structure = AseAtomsAdaptor.get_structure(slab)
    add_to_dict(structure, 'graphite', millers)


# potassium mercury almalgum

material = ['mp-11462']
name = ['KHg']

bulk = rest.get_structure_by_material_id('mp-11462')
bulk = AseAtomsAdaptor.get_atoms(bulk)
# correct spacing between layers
for millers in [(1,0,0)]:
    slab = surfaces_with_termination(bulk, millers, 2, vacuum=6)[0]
    structure = AseAtomsAdaptor.get_structure(slab)
    add_to_dict(structure, 'KHg', millers)


# TiC

material = ['mp-631']
name = ['TiC']

bulk = rest.get_structure_by_material_id('mp-631')
bulk = AseAtomsAdaptor.get_atoms(bulk)
# correct spacing between layers
for millers in [(1,1,0)]:
    slab = surfaces_with_termination(bulk, millers, 4, vacuum=6)[0]
    slab *= (2,2,1)
    structure = AseAtomsAdaptor.get_structure(slab)
    add_to_dict(structure, 'TiC', millers)

# V8N

"""
material = ['mp-1188283']
name = ['VN']

bulk = rest.get_structure_by_material_id('mp-1188283')
bulk = AseAtomsAdaptor.get_atoms(bulk)
# correct spacing between layers
for millers in [(1,0,0), (1,1,0), (1,1,1)]:
    slab = surfaces_with_termination(bulk, millers, 3, vacuum=6)[0]
    structure = AseAtomsAdaptor.get_structure(slab)
    add_to_dict(structure, 'VN', millers)
"""

# KCl
# doesn't converge in abinit/sparc

"""
material = ['mp-23193']
name = ['KCl']

bulk = rest.get_structure_by_material_id('mp-23193')
bulk = AseAtomsAdaptor.get_atoms(bulk)
# correct spacing between layers
for millers in [(1,0,0), (1,1,0), (1,1,1)]:
    slab = surfaces_with_termination(bulk, millers, 3, vacuum=6)[0]
    slab *= (2,2,1)
    structure = AseAtomsAdaptor.get_structure(slab)
    add_to_dict(structure, 'KCl', millers)
"""

"""
material = ['mp-644481']
name = ['ScO']

bulk = rest.get_structure_by_material_id('mp-644481')
bulk = AseAtomsAdaptor.get_atoms(bulk)
# correct spacing between layers

for millers, t in zip([(1,1,1)], (0,0,0)):
    slab = surfaces_with_termination(bulk, millers, 3, vacuum=6)[t]
    #slab *= (2,2,1)
    structure = AseAtomsAdaptor.get_structure(slab)
    add_to_dict(structure, 'ScO', millers)
"""

material = ['mp-1279']
name = ['TaN']

bulk = rest.get_structure_by_material_id('mp-1279')
bulk = AseAtomsAdaptor.get_atoms(bulk)
# correct spacing between layers

for millers in [(1,1,0)]:
    slab = surfaces_with_termination(bulk, millers, 3, vacuum=6)[0]
    #slab *= (2,2,1)
    structure = AseAtomsAdaptor.get_structure(slab)
    add_to_dict(structure, 'TaN', millers)


# Iron Oxide
"""
bulk = rest.get_structure_by_material_id('mp-19306') # Fe3O4
bulk = AseAtomsAdaptor.get_atoms(bulk)
#bulk.positions += np.array([1,1,1]) * 0.5
# this next function is bugged, future versions of ASE might not recreat this stucture
slab = surfaces_with_termination(bulk, (1,1,1), 4, vacuum=6, termination='O')[0]
structure = AseAtomsAdaptor.get_structure(slab)
add_to_dict(structure, 'Fe3O4', (1,1,1))
"""

# rutile oxides


from ase.spacegroup import crystal

# values obtained from materials project

material = ['mp-2657', 'mp-856', 'mp-19094']
names = ['TiO2',  'VO2']
metals = ['Ti',  'V']
lattice_constants = [(4.607, 2.992), (4.829, 3.243), (4.507, 3.044)]


for name, element, lattice in zip(names, metals, lattice_constants):
    a = lattice[0]
    c = lattice[1]
    bulk = crystal([element, 'O'], basis=[(0, 0, 0), (0.3, 0.3, 0.0)],
                   spacegroup=136, cellpar=[a, a, c, 90, 90, 90])
    # this next function is bugged, future versions of ASE might not recreate this stucture
    slab = surfaces_with_termination(bulk, (1,1,0), 4, vacuum=6, termination='O')[1]
    if element == 'Ti':
        slab *= (1,2,1)
        del slab[45]
        name += '_defect'
    if element == 'Sn':
        ads = molecule('NH3')
        ads.rotate(180, v='y')
        slab *= (1,2,1)
        add_adsorbate(slab, ads, .7, (0,3.243))
        name += '_NH3ads'
    if element == 'V':
        ads = molecule('N2')
        #ads.rotate(180, v='y')
        slab *= (1,2,1)
        add_adsorbate(slab, ads, 1.9, (0,3.243))
        name += '_N2ads'

    #print(len(slab))
    structure = AseAtomsAdaptor.get_structure(slab)
    add_to_dict(structure, name, (1,1,0))


# more ruiltes

material = ['mp-2657']
names = ['TiO2']
metals = ['Ti']
lattice_constants = [(4.607, 2.992)]


for name, element, lattice in zip(names, metals, lattice_constants):
    a = lattice[0]
    c = lattice[1]
    bulk = crystal([element, 'O'], basis=[(0, 0, 0), (0.3, 0.3, 0.0)],
                   spacegroup=136, cellpar=[a, a, c, 90, 90, 90])
    slab = surfaces_with_termination(bulk, (1,1,0), 4, vacuum=6, termination='O')[1]
    ads = molecule('NH3')
    ads.rotate(180, v='y')
    slab *= (1,2,1)
    add_adsorbate(slab, ads, 1, (0,3.243))
    name += '_NH3ads'
    structure = AseAtomsAdaptor.get_structure(slab)
    add_to_dict(structure, name, (1,1,0))

# more oxides

# ZnO
"""
material = ['mp-2133']
name = ['ZnO']

bulk = rest.get_structure_by_material_id('mp-2133') # ZnO
bulk = AseAtomsAdaptor.get_atoms(bulk)
#bulk.positions += np.array([1,1,1]) * 0.5
# this next function is bugged, future versions of ASE might not recreat this stucture
for millers in [(1,1,0), (1,0,0)]:
    slab = surfaces_with_termination(bulk, millers, 4, vacuum=6)[0]
    structure = AseAtomsAdaptor.get_structure(slab)
    add_to_dict(structure, 'ZnO', millers)
"""

material = ['mp-1007776']
name = ['TlP']

bulk = rest.get_structure_by_material_id('mp-1007776') 
bulk = AseAtomsAdaptor.get_atoms(bulk)
#bulk.positions += np.array([1,1,1]) * 0.5
# this next function is bugged, future versions of ASE might not recreat this stucture
for millers in [(1,1,0)]:
    slab = surfaces_with_termination(bulk, millers, 4, vacuum=6)[0]
    structure = AseAtomsAdaptor.get_structure(slab)
    add_to_dict(structure, 'TlP', millers)


material = ['mp-7586']
name = ['Ni3C']

bulk = rest.get_structure_by_material_id('mp-7586') 
bulk = AseAtomsAdaptor.get_atoms(bulk)
#bulk.positions += np.array([1,1,1]) * 0.5
# this next function is bugged, future versions of ASE might not recreat this stucture
for millers in [(1,0,0)]:
    slab = surfaces_with_termination(bulk, millers, 2, vacuum=6)[0]
    structure = AseAtomsAdaptor.get_structure(slab)
    add_to_dict(structure, 'Ni3C', millers)


material = ['mp-2662']
name = ['MnP']

bulk = rest.get_structure_by_material_id('mp-2662') 
bulk = AseAtomsAdaptor.get_atoms(bulk)
#bulk.positions += np.array([1,1,1]) * 0.5
# this next function is bugged, future versions of ASE might not recreat this stucture
for millers in [(1,1,1)]:
    slab = surfaces_with_termination(bulk, millers, 2, vacuum=6)[0]
    structure = AseAtomsAdaptor.get_structure(slab)
    add_to_dict(structure, 'MnP', millers)


material = ['mp-1206445']
name = ['ReW3']

bulk = rest.get_structure_by_material_id('mp-1206445')
bulk = AseAtomsAdaptor.get_atoms(bulk)
#bulk.positions += np.array([1,1,1]) * 0.5
# this next function is bugged, future versions of ASE might not recreat this stucture
for millers in [(0,0,1)]:
    slab = surfaces_with_termination(bulk, millers, 2, vacuum=6)[0]
    structure = AseAtomsAdaptor.get_structure(slab)
    add_to_dict(structure, 'ReW3', millers)


# Al2O3


material = ['mp-1143']
name = ['Al2O3']

bulk = rest.get_structure_by_material_id('mp-1143') # Al2O3
bulk = AseAtomsAdaptor.get_atoms(bulk)
#bulk.positions += np.array([1,1,1]) * 0.5
# this next function is bugged, future versions of ASE might not recreat this stucture
for millers, i in zip([(1,1,1)], [0]):
    slab = surfaces_with_termination(bulk, millers, 3, vacuum=6, termination='O')[i]
    structure = AseAtomsAdaptor.get_structure(slab)
    add_to_dict(structure, 'Al2O3', millers)

material = ['mp-1226211']
name = ['FeCr']

bulk = rest.get_structure_by_material_id('mp-1226211')
bulk = AseAtomsAdaptor.get_atoms(bulk)
#bulk.positions += np.array([1,1,1]) * 0.5
# this next function is bugged, future versions of ASE might not recreat this stucture
for millers, i in zip([(1,0,0)], [0]):
    slab = surfaces_with_termination(bulk, millers, 3, vacuum=6)[0]
    structure = AseAtomsAdaptor.get_structure(slab)
    add_to_dict(structure, 'FeCr', millers)

material = ['mp-30726']
name = ['YHg']

bulk = rest.get_structure_by_material_id('mp-30726') 
bulk = AseAtomsAdaptor.get_atoms(bulk)
#bulk.positions += np.array([1,1,1]) * 0.5
# this next function is bugged, future versions of ASE might not recreat this stucture
for millers, i in zip([(1,1,0)], [0]):
    slab = surfaces_with_termination(bulk, millers, 2, vacuum=6)[0]
    structure = AseAtomsAdaptor.get_structure(slab)
    add_to_dict(structure, 'YHg', millers)


bulk = rest.get_structure_by_material_id('mp-1008918')
bulk = AseAtomsAdaptor.get_atoms(bulk)
#bulk.positions += np.array([1,1,1]) * 0.5
# this next function is bugged, future versions of ASE might not recreat this stucture
for millers, i in zip([(1,1,1)], [0]):
    slab = surfaces_with_termination(bulk, millers, 3, vacuum=6)[0]
    structure = AseAtomsAdaptor.get_structure(slab)
    slab *= (2,1,1)
    add_to_dict(structure, 'CdN', millers)


bulk = rest.get_structure_by_material_id('mp-11460')
bulk = AseAtomsAdaptor.get_atoms(bulk)
#bulk.positions += np.array([1,1,1]) * 0.5
# this next function is bugged, future versions of ASE might not recreat this stucture
for millers, i in zip([(1,0,0)], [0]):
    slab = surfaces_with_termination(bulk, millers, 3, vacuum=6)[0]
    structure = AseAtomsAdaptor.get_structure(slab)
    slab *= (2,1,1)
    add_to_dict(structure, 'HfTc', millers)



# 2D materials

# Graphene
material = ['mp-48']
name = ['graphite']

bulk = rest.get_structure_by_material_id('mp-48')
bulk = AseAtomsAdaptor.get_atoms(bulk)
del bulk[0]
del bulk[1]
new_cell = bulk.get_cell()
new_cell = np.vstack((new_cell[:2],np.array([0,0,12])))
bulk.set_cell(new_cell)
bulk *= (3,3,1)
bulk.center()
# correct spacing between layers
structure = AseAtomsAdaptor.get_structure(bulk)
add_to_dict(structure, 'graphene', millers)

# Pt on Graphene
material = ['mp-48']
name = ['graphite']

bulk = rest.get_structure_by_material_id('mp-48')
bulk = AseAtomsAdaptor.get_atoms(bulk)
del bulk[0]
del bulk[1]
new_cell = bulk.get_cell()
new_cell = np.vstack((new_cell[:2],np.array([0,0,12])))
bulk.set_cell(new_cell)
bulk *= (3,3,1)
bulk.center()
# correct spacing between layers
add_adsorbate(bulk, 'Ir', 2.0, (2.4675, 2.1365))
structure = AseAtomsAdaptor.get_structure(bulk)
add_to_dict(structure, 'Ir_graphene', millers)

# Boro-nitro-graphene

material = ['mp-642462']
name = ['B3C10N3']

bulk = rest.get_structure_by_material_id('mp-642462')
bulk = AseAtomsAdaptor.get_atoms(bulk)
bulk.rotate(90, 'x', rotate_cell=True)
cell = bulk.get_cell()
# get rid of the extra layer
bulk = [a for a in bulk if a.position[2] < 3]
bulk = Atoms(bulk)
# fix the height
cell = cell.array
cell[1,2] = 12
bulk.set_cell(cell)
bulk.center()
structure = AseAtomsAdaptor.get_structure(bulk)
add_to_dict(structure, 'BN_graphene', millers)

# MoS2 sheets

material = ['mp-1143']
name = ['Al2O3']

bulk = rest.get_structure_by_material_id('mp-1018809') # Al2O3
bulk = AseAtomsAdaptor.get_atoms(bulk)
del bulk[4]
del bulk[3]
del bulk[1]

new_cell = bulk.get_cell()
new_cell = np.vstack((new_cell[:2],np.array([0,0,12])))
bulk.set_cell(new_cell)
bulk *= (2,2,1)

structure = AseAtomsAdaptor.get_structure(bulk)
#view(slab)
add_to_dict(structure, 'MoS2', millers)

lens = []
for n, s in structures.items():
    lens.append(len(Structure.from_dict(s)))
    #print(len(Structure.from_dict(s)))

#print(([a > 30 for a in lens]).count(True))
#print(([a > 100 for a in lens]).count(True))
print(len(structures))


#print(structures.keys())

json.dump(params_dict, open('../../test_set/surface_parameters.json', 'w'))
json.dump(structures, open('../../test_set/surface.json', 'w'))
