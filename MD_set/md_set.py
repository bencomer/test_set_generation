from ase.build import molecule
from ase.atoms import Atoms
from ase.atom import Atom
from ase.units import Bohr
from ase.io import read
from ase.visualize import view
from pymatgen.io.ase import AseAtomsAdaptor
from lammps_interface.tools import make_box_of_molecules
from ase.data.pubchem import pubchem_atoms_search
import json
from ase.collections import g2



def add_to_dict(atoms, name,  temp):
    structure = AseAtomsAdaptor.get_structure(atoms)
    formula = structure.composition.reduced_formula
    structures[name+'_MD'] = structure.as_dict()
    sparc_p = sprc_params.copy()
    esp_p = esp_params.copy()

    params_dict[name+'_MD'] = {}
    params_dict[name+'_MD']['sparc'] = sparc_p.copy()
    params_dict[name+'_MD']['sparc']['ION_TEMP'] = temp
    params_dict[name+'_MD']['espresso'] = esp_p.copy()
    params_dict[name+'_MD']['espresso']['tempw'] = temp
    params_dict[name+'_MD']['abinit'] = abinit_params.copy()
    params_dict[name+'_MD']['abinit']['ionmov'] = 6
    params_dict[name+'_MD']['abinit']['mdtemp'] = temp
    if name == 'iron_interface':
        params_dict[name+'_MD']['sparc']['TOL_SCF'] = 1e-5
        #del params_dict[name+'_MD']['sparc']['TOL_SCF_QE']
        params_dict[name+'_MD']['sparc']['SPIN_TYP'] = 2


structures = {}
params_dict = {}

sprc_params = {}
sprc_params['h'] = 0.2 * Bohr
sprc_params['MAXIT_SCF'] = 1000
#sprc_params['TOL_SCF_QE'] = 1e-6
sprc_params['TOL_SCF'] = 5e-5
sprc_params['BC'] = [True] * 3
sprc_params['MD_FLAG'] = 1
sprc_params['MD_METHOD'] = 'NVE'
sprc_params['EXCHANGE_CORRELATION'] = 'LDA_PZ'
sprc_params['MIXING_VARIABLE'] = 'density'
sprc_params['MIXING_PRECOND'] = 'none'
sprc_params['ELEC_TEMP'] = 1160.45
sprc_params['KPOINT_SHIFT'] = ' 0 0 0'
#params['SCF_ENERGY_ACC'] = scf_acc * Hartree
#params['MIXING_VARIABLE'] = 'potential'

abinit_params = {'xc':'LDA_PZ',
         'autoparal':4,
         'optcell':0,
         'nstep':500,
         'nsym':1,
         #'pps':'oncv',
         'pps':'pot',
         'ecut': 2721.14,
         'tsmear':0.1,
         'chksymbreak':0,
         'toldfe':1.0e-6,
         'occopt':7,
}


esp_params = dict(#xc='PZ',
        kpts=(1, 1, 1), #only need 1 kpt in z-direction
        pw=1360,
        dw=13600,
        #spinpol=True,
        calculation='md',
        ion_dynamics='verlet',
        nosym=True,
        convergence={'energy':1e-6,
                    #'mixing':0.1,
                    'maxsteps':1000,
                    'diag':'david'},
        startingwfc='random',
        smearing='gaussian', #gaussian electron smearing
        sigma=0.1, #smearing width
        dipole={'status':False}, #dipole corrections True turns them on
        outdir ='.'
        )


# He

# PV=nRT @ 100 bar, 100 C
atoms = make_box_of_molecules(Atoms('He'), 75, [33.811] * 3)
add_to_dict(atoms, 'He', 100 + 273.15)

# water

# liquid water at room temp, 55.555 mol/L

atoms = make_box_of_molecules(molecule('H2O'), 113, [15] * 3)
add_to_dict(atoms, 'water', 300)

# CO2

atoms = make_box_of_molecules(molecule('CO2'), 25, [20] * 3)
print(len(atoms))
add_to_dict(atoms, 'CO2', 70 + 273.15)

# Al

atoms = read('Al.xyz')
add_to_dict(atoms, 'Al', 10000)


# natural gas


CH4 = pubchem_atoms_search('methane')
C2H6 = pubchem_atoms_search('ethane')
C3H8 = pubchem_atoms_search('propane')


atoms = make_box_of_molecules([CH4, C2H6, C3H8], [18, 6, 6], [15] * 3)
add_to_dict(atoms, 'LNG', 273.15 - 162)


# Battery

with open('Li-ion_electrolyte_interface_318.txt', 'r') as f:
    txt = f.read()
    lattice = txt.split('lattice vectors (au):')[1]
    lattice, positions = txt.split('Atom types and positions (au):')

lattice = lattice.splitlines()[1:]
lattice = [[float(b) * Bohr for b in a.split()] for a in lattice]

sym_dict = {'carbon':'C',
            'oxygen':'O',
            'hydrogen':'H',
            'lithium':'Li',
            'fluorine':'F',
            'phosphorus':'P'}

positions = positions.splitlines()[1:]

atoms = Atoms()

for line in positions:
    element, x, y, z = line.split()
    position = [float(a) * Bohr for a in [x, y, z]]
    atoms += Atom(symbol=sym_dict[element], position=position)
atoms.set_cell(lattice)
atoms.set_pbc([True] * 3)

add_to_dict(atoms, 'electrolyte', 273.15 + 23)

# iron interface

#atoms = read('interface.traj')
#add_to_dict(atoms, 'iron_interface', 273.15 + 23)


json.dump(structures, open('../../test_set/MD.json', 'w'))
json.dump(params_dict, open('../../test_set/MD_parameters.json', 'w'))

