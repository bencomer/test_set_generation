from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen import Structure, Lattice, MPRester, Molecule
from ase.data.pubchem import pubchem_atoms_search
from ase.units import Bohr
from ase.data import chemical_symbols, reference_states
from ase.build import bulk
from ase.build import nanotube
from ase.build.ribbon import graphene_nanoribbon
from ase.visualize import view
from ase.atom import Atom
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
        if i in [1,2]:
            kpts.append(1)
        else:
            l = np.linalg.norm(atoms.cell[i])
            kpts.append(math.ceil(kpt_density / l))
    return kpts

def find_shift_abinit(kpts):
    shift = ' '
    for kpt in kpts:
        if kpt % 2:
            shift += '0 '
        else:
            shift += '0.5 '
    return shift


def add_to_dict(structure, name):
    # formula = structure.composition.reduced_formula

    slab = structure.as_dict()
    if 'scale_factor' in slab:
        slab['scale_factor'] = slab['scale_factor'].tolist()

    structures[name] = slab
    kpts = select_kpts(structure)
    sparc_p = sprc_params.copy()
    esp_p = esp_params.copy()
    sparc_p['KPOINT_GRID'] = kpts
    esp_p['kpts'] = kpts
    abinit_p = abinit_params.copy()
    abinit_p['kptopt'] = 1
    if 'MoS2' in name or 'polyacetylene' in name: #or 'S' in syms or formula in ['WC', 'PSe']:
        sparc_p['EXCHANGE_CORRELATION'] = 'GGA_PBE'
        abinit_p['xc'] = 'PBE'
        esp_p['xc'] = 'PBE'
    abinit_p['ngkpt'] = ' '.join([str(a) for a in kpts])
    abinit_p['nshiftk'] = 1
    abinit_p['shiftk'] = find_shift_abinit(kpts)
    if 'PVC' in name:
        abinit_p['diemix'] = 0.1
    params_dict[name] = {}
    params_dict[name]['sparc'] = sparc_p
    params_dict[name]['espresso'] = esp_p
    params_dict[name]['abinit'] = abinit_p.copy()


sprc_params = {}
sprc_params['h'] = 0.2 * Bohr
sprc_params['MAXIT_SCF'] = 1000
sprc_params['TOL_SCF'] = 5e-5
#sprc_params['TOL_SCF_QE'] = 1e-6

sprc_params['BC'] = [True] * 3
sprc_params['EXCHANGE_CORRELATION'] = 'LDA_PZ'
sprc_params['MIXING_VARIABLE'] = 'density'
sprc_params['MIXING_PRECOND'] = 'none'
sprc_params['ELEC_TEMP'] = 1160.45
sprc_params['KPOINT_SHIFT'] = ' 0 0 0'
#params['MIXING_VARIABLE'] = 'potential'


esp_params = dict(#xc='PZ',
        kpts=(1, 1, 1), #only need 1 kpt in z-direction
        pw=1360,
        dw=13600,
        calculation='scf',
        #spinpol=True,
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
         'optcell':0,
         'nsym':1,
         #'pps':'oncv',
         'nstep':500,
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


wire = nanotube(6, 0, length=1, vacuum=6)
wire.rotate(90, 'y', rotate_cell=True)
wire.set_cell([wire.cell[2], wire.cell[1], wire.cell[0] * -1])
wire.center()
#view(wire)
add_to_dict(AseAtomsAdaptor.get_structure(wire), 'CNT')


# nano-ribbon

ribbon = graphene_nanoribbon(2, 2, saturated=True, vacuum=6)
ribbon.rotate(90, 'y', rotate_cell=True)
ribbon.set_cell([ribbon.cell[2], ribbon.cell[1], ribbon.cell[0] * -1])
ribbon.center()
add_to_dict(AseAtomsAdaptor.get_structure(ribbon), 'graphene_ribbon')


from ase import Atoms
d = 2.9
L = 12.0
"""
for metal in ['Au', 'Pt', 'Fe']:
    wire = Atoms(metal,
                 positions=[[0, L / 2, L / 2]],
                 cell=[d, L, L],
                 pbc=[1, 0, 0])
    add_to_dict(AseAtomsAdaptor.get_structure(wire), metal + '_Wire')
"""

for metal in ['Au', 'Pt', 'Mo']:
    atoms = bulk(metal, cubic=True)
    if metal == 'Fe':
        atoms *=(1,2,2)
    atoms.set_cell([atoms.cell[0][0]] + [12] * 2)
    atoms.center()
    add_to_dict(AseAtomsAdaptor.get_structure(atoms), metal + '_wire')


# this code is about to get pretty serious if you are trying to follow it
atoms = pubchem_atoms_search('ethane')
vector = atoms[1].position - atoms[0].position
methyl = Atoms([atoms[a] for a in [1,6,5,7]])
del atoms[7]
del atoms[6]
del atoms[5]
del atoms[1]
methyl.rotate(120, v=vector, center=methyl[0].position)
atoms += methyl

del atoms[7]
del atoms[3]
atoms.rotate(180, v=vector)
atoms.positions -= atoms[0].position


bl = np.linalg.norm(vector)

cell = np.array([
                 #[1,0,0],
                 #[bl * (np.cos(np.pi/6) + 1), bl * np.sin(np.pi/6), 0],
                 [bl * (np.sin(np.pi/6)), bl * (np.cos(np.pi/6) + 1), 0],
                 [0,40,0],

                 [0,0,15]])



#atoms.rotate(180, v=vector)
atoms.rotate(90, v = [0,0,1])
n = np.cos(np.pi/6) * 1.512 * 2

#cell *= np.array([bl,n,12])
atoms.set_cell(cell)
atoms *= (1,1,1)
atoms.center()
atoms *= (2,1,1)
add_to_dict(AseAtomsAdaptor.get_structure(atoms), 'polyethylene')



# PVC (syndiatactic)

atoms = pubchem_atoms_search('ethane')
vector = atoms[1].position - atoms[0].position
methyl = Atoms([atoms[a] for a in [1,6,5,7]])
del atoms[7]
del atoms[6]
del atoms[5]
del atoms[1]
methyl.rotate(120, v=vector, center=methyl[0].position)
atoms += methyl

del atoms[7]
del atoms[3]
atoms.rotate(180, v=vector)
atoms.positions -= atoms[0].position


bl = np.linalg.norm(vector)

cell = np.array([
                 #[1,0,0],
                 #[bl * (np.cos(np.pi/6) + 1), bl * np.sin(np.pi/6), 0],
                 [bl * (np.sin(np.pi/6)), bl * (np.cos(np.pi/6) + 1), 0],
                 [0,40,0],

                 [0,0,15]])



#atoms.rotate(180, v=vector)
atoms.rotate(90, v = [0,0,1])
n = np.cos(np.pi/6) * 1.512 * 2

#cell *= np.array([bl,n,12])
atoms.set_cell(cell)
atoms *= (2,1,1)
atoms.center()
atoms[1].symbol = 'Cl'
atoms[8].symbol = 'Cl'
atoms.set_distance(0, 1, 1.712, fix=0)
atoms.set_distance(6, 8, 1.712, fix=0)
atoms.center()
atoms.wrap(pbc = [True, False, False])
#atoms *= (2,1,1)


add_to_dict(AseAtomsAdaptor.get_structure(atoms), 'PVC')
#view(atoms)

# PTFE

atoms = pubchem_atoms_search('ethane')
vector = atoms[1].position - atoms[0].position
methyl = Atoms([atoms[a] for a in [1,6,5,7]])
del atoms[7]
del atoms[6]
del atoms[5]
del atoms[1]
methyl.rotate(120, v=vector, center=methyl[0].position)
atoms += methyl

del atoms[7]
del atoms[3]
atoms.rotate(180, v=vector)
atoms.positions -= atoms[0].position


bl = np.linalg.norm(vector)

cell = np.array([
                 #[1,0,0],
                 #[bl * (np.cos(np.pi/6) + 1), bl * np.sin(np.pi/6), 0],
                 [bl * (np.sin(np.pi/6)), bl * (np.cos(np.pi/6) + 1), 0],
                 [0,40,0],

                 [0,0,15]])



#atoms.rotate(180, v=vector)
atoms.rotate(90, v = [0,0,1])
n = np.cos(np.pi/6) * 1.512 * 2

#cell *= np.array([bl,n,12])
atoms.set_cell(cell)
atoms *= (1,1,1)
atoms.center()
for index in [1,2]:
    atoms[index].symbol = 'F'
    atoms.set_distance(0, index, 1.35, fix=0)

for index in [4,5]:
    atoms[index].symbol = 'F'
    atoms.set_distance(3, index, 1.35, fix=0)
#atoms[8].symbol = 'Cl'
#atoms.set_distance(0, 1, 1.712, fix=0)
#atoms.set_distance(6, 8, 1.712, fix=0)
atoms.wrap(pbc=[True, False, False])
atoms *= (2,1,1)

add_to_dict(AseAtomsAdaptor.get_structure(atoms), 'PTFE')


# polyacetylene

atoms = pubchem_atoms_search('ethene')
vector = atoms[1].position - atoms[0].position

del atoms[4]
del atoms[2]
atoms.rotate(180, v=vector)
atoms.positions -= atoms[0].position


bl = np.linalg.norm(vector)



bl = np.linalg.norm(vector)

cell = np.array([
                 #[1,0,0],
                 [-1 * bl * (np.cos(np.pi/6) + 1), bl * np.sin(np.pi/6), 0],
                 #[bl * (np.sin(np.pi/6)), bl * (np.cos(np.pi/6) + 1), 0],
                 [0,15,0],

                 [0,0,15]])



#atoms.rotate(180, v=vector)
#atoms.rotate(90, v = [0,0,1])
#n = np.cos(np.pi/6) * 1.512 * 2

#cell *= np.array([bl,n,12])
atoms.set_cell(cell)
atoms *= (1,1,1)
atoms.center()
atoms.set_cell([atoms.cell.array[0] * -1,atoms.cell.array[1],atoms.cell.array[2]])
atoms.wrap(pbc=[True, False, False])
atoms *= (2,1,1)
#view(atoms)


add_to_dict(AseAtomsAdaptor.get_structure(atoms), 'polyacetylene')

# MoS2

with open('MoS2') as f:
    txt = f.read()

header, Mo, S = txt.split(':')

MoS2 = Atoms()

for element, block in zip(['Mo', 'S'], [Mo, S]):
    b = block.splitlines()
    for line in b[1:-1]:
        MoS2 += Atom(position=[float(a) * Bohr for a in line.split()], symbol=element)
    
MoS2.set_cell([27,27,6.017129910500000 * Bohr])
MoS2.center()
MoS2.rotate(90, 'y', rotate_cell=True)
#MoS2.set_cell([MoS2.cell[0] * -1, MoS2.cell[1], MoS2.cell[2]])
MoS2.set_cell([MoS2.cell[2], MoS2.cell[1], MoS2.cell[0] * -1])
MoS2.center()
#view(MoS2)
add_to_dict(AseAtomsAdaptor.get_structure(MoS2), 'MoS2_nanotube')

lens = []
for n, s in structures.items():
    lens.append(len(Structure.from_dict(s)))
    #print(len(Structure.from_dict(s)))

print(([a > 30 for a in lens]).count(True))
print(([a > 100 for a in lens]).count(True))



print(len(structures))
json.dump(structures, open('../../test_set/wire.json', 'w'))
json.dump(params_dict, open('../../test_set/wire_parameters.json', 'w'))
