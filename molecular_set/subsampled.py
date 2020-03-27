from ase.build import molecule
from ase.units import Bohr
from ase.atom import Atom
from ase.visualize import view
from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.core import Structure
from lammps_interface.tools import make_box_of_molecules
from ase.data.pubchem import pubchem_atoms_search
from ase.units import Bohr
from ase.atoms import Atoms
import json
from ase.collections import g2



def add_to_dict(structure, name=None):
    formula = structure.composition.reduced_formula
    if name is not None:
        formula = name
    structures[formula+'_molecule'] = structure.as_dict()
    sparc_p = sprc_params.copy()
    esp_p = esp_params.copy()
    abinit_p = abinit_params.copy()
    if 'O' in formula:
        print('yeet')
    if formula in ['O2', 'SO2']:
        sparc_p['EXCHANGE_CORRELATION'] = 'GGA_PBE'
        abinit_p['xc'] = 'PBE'
        esp_p['xc'] = 'PBE'
    if formula == 'Pt_cluster_molecule':
        abinit_p['diemix'] = 0.1


    for site in structure:
        if hasattr(site, 'magmom'):
            print(formula)
            abinit_p['nsppol'] = 2
            sparc_p['SPIN_TYP'] = 2
            esp_p['spinpol'] = True

    params_dict[formula+'_molecule'] = {}
    params_dict[formula+'_molecule']['sparc'] = sparc_p.copy()
    if 'box' in formula:
        params_dict[formula+'_molecule']['sparc']['BC'] = [True] * 3
    if formula == 'Pt_cluster':
        params_dict[formula+'_molecule']['sparc']['MIXING_VARIABLE'] = 'density'
        params_dict[formula+'_molecule']['sparc']['MIXING_PRECOND'] = 'kerker'
    params_dict[formula+'_molecule']['espresso'] = esp_p.copy()
    params_dict[formula+'_molecule']['abinit'] = abinit_p.copy()


structures = {}
params_dict = {}

sprc_params = {}
sprc_params['h'] = 0.4 * Bohr
sprc_params['MAXIT_SCF'] = 1000
#sprc_params['TOL_SCF_QE'] = 1e-6
sprc_params['TOL_SCF'] = 5e-5
sprc_params['BC'] = [True] * 3
sprc_params['EXCHANGE_CORRELATION'] = 'LDA_PZ'
#sprc_params['MIXING_VARIABLE'] = 'density'
#sprc_params['MIXING_PRECOND'] = 'none'
sprc_params['ELEC_TEMP'] = 1160.45
sprc_params['KPOINT_SHIFT'] = ' 0 0 0'
#params['SCF_ENERGY_ACC'] = scf_acc * Hartree
#params['MIXING_VARIABLE'] = 'potential'

esp_params = dict(#xc='PZ',
        kpts=(1, 1, 1), #only need 1 kpt in z-direction
        pw=300,
        dw=3000,
        #spinpol=True,
        calculation='scf',
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



boring_list = ['methylenecyclopropane', 'isobutene', 'CH2NHCH2', 'CH2_s3B1d',
               'C3H4_C3v', 'C3H4_C2v', 'C5H5N', 'CH3CH2OCH3', 'C4H4O', 'CH',
               'NH', 'NH2', '2-butyne', 'OCHCHO', 'C4H4NH', 'SiH2_s3B1d',
               'CH3CO', 'H3CNH2', 'CH3CH2OH', 'C2H3', 'CH3CN', 'CH3ONO',
               'C3H6_D3h', 'trans-butane', 'CH2OCH2', 'C6H6', 'CH3CONH2',
               'H2CCHCN', 'butadiene', 'H2CO', 'CH3COOH', 'HCF3', 'SiH2_s1A1d',
               'C4H4S', 'N2H4', 'CH2_s1A1d', 'CH3CH2SH', 'CH3NO2', 'C4H4O', 
               'CH3O', 'CH3OH', 'isobutane', 'CH3CH2O', 'H2CCHF', 'SiO', 
               'C3H4_D2d', 'COF2', '2-butyne', 'C2H5', 'F2O', 'H2CCl2', 'CF3CN',
               'C2H6NH',  'ClO', 'CH3COCl', 'C3H9N', 'C3H6_Cs', 'HCOOCH3', 'CCH',
               'Si2', 'C2H6SO', 'C5H8', 'H2CF2', 'CH2SCH2', 'C2Cl4', 'C3H4_C3v',
               'CH3COCH3', 'H2CCO', 'CH3CH2NH2', 'H2O2', 'C3H4_C2v', 'CH3', 'O3',
               'cyclobutane', 'cyclobutene', 'CH3CH2Cl', 'CH3SH', 'C3H8',
               'CH3SCH3', 'NCCN', 'CH3CHO', 'CH3OCH3', 'CH3COF', 'CH3OCH3', 
               'C3H7Cl', 'H2CCHCl', 'ClNO', 'C2H6CHOH', 'SiCl4', 'NF3',
               'Si2H6',  'C2H6', 'P2', 'CCl4', 'BF3', 'F2', 'PF3', 'CO2', 'CH4',
               'H2O', 'SiF4', 'HOCl', 'Cl2', 'HCCL3', 'CS2', 'PH3', 'CF4', 'CF2',
               'H2C', 'ClF','Cl2', 'CSO', 'CF2', 'CS','HC','BCl3','H2O2',
               'H3CCl','HCl', 'H2CO2', 'HCCl3','H2O2','HC','SiH6C', 'N2O', 'HClO',
               'OCS', 'CH3SiH3', 'ClF3', 'C2F4', 'CO', 'NH3', 'C2H2', 'HF', 'H2O2',
               'SO2', 'C2H4', 'SH2', 'C2H3', 'O2']

molecule_list = []
mol_names = []

for mol in g2.names:
    atoms = molecule(mol)

    # remove the single atoms
    if len(atoms) == 1:
        continue
    # remove the spin polarized systems
    if atoms.get_initial_magnetic_moments().any() != 0.:
        if mol not in ['NO', 'OH']:
            continue
    # capriciously removing systems with Li, Na, and Al:
    for element in ['Li', 'Na', 'Al']:
        if element in atoms.get_chemical_formula():
            boring_list.append(mol)
    # remove random ones that seem superfluous
    if mol in boring_list:
        continue
    mol_names.append(mol)
    atoms.set_cell([10,10,10])
    atoms.center()

    struct = AseAtomsAdaptor.get_structure(atoms)
    if atoms.get_initial_magnetic_moments().any() != 0.:
        struct.add_site_property('magmom', atoms.get_initial_magnetic_moments())
    add_to_dict(struct)
    molecule_list.append(atoms)

# random molecules to add more elements
atoms = read('Xe.sdf')
atoms.set_cell([10,10,10])
for index in [1,2,3,4]:
    atoms.set_distance(0, index, 1.96, fix=0)
atoms.center()

struct = AseAtomsAdaptor.get_structure(atoms)
add_to_dict(struct)

atoms = molecule('HF')
for atom in atoms:
    if atom.symbol == 'F':
        atom.symbol = 'Br'

atoms.set_distance(0, 1, 1.43)
atoms.set_cell([10,10,10])
atoms.center()

struct = AseAtomsAdaptor.get_structure(atoms)
add_to_dict(struct)

atoms = pubchem_atoms_search(cid=807)  # I2

atoms.set_cell([10,10,10])
atoms.center()

struct = AseAtomsAdaptor.get_structure(atoms)
add_to_dict(struct, 'I2')



# nano-particles
import ase
from ase.cluster.cubic import FaceCenteredCubic

surfaces = [(1, 0, 0), (1, 1, 0), (1, 1, 1)]
layers = [3, 3, 3]
lc = 3.61000
atoms = FaceCenteredCubic('Cu', surfaces, layers, latticeconstant=lc, vacuum=10)
add_to_dict(AseAtomsAdaptor.get_structure(atoms), 'Cu_cluster')


from ase.cluster import Octahedron

atoms = Octahedron(['Pt','Au'], length=3, alloy=True, latticeconstant=4)
atoms.set_cell([15] * 3)
atoms.center()
add_to_dict(AseAtomsAdaptor.get_structure(atoms), 'PtAu_cluster')


from ase.cluster import Icosahedron

atoms = Icosahedron('Pt', noshells=3,)
atoms.set_cell([21] * 3)
atoms.center()
add_to_dict(AseAtomsAdaptor.get_structure(atoms), 'Pt_cluster')


surfaces = [(1, 0, 0), (1, 1, 1), (1, -1, 1)]
layers = [3, 1, -1]
atoms = FaceCenteredCubic('Au', surfaces, layers)
atoms.set_cell([10] * 3)
add_to_dict(AseAtomsAdaptor.get_structure(atoms), 'Au_cluster')
atoms.center()

"""
from ase.cluster import  Decahedron

surfaces = [(1, 0, 0), (1, 1, 1), (1, -1, 1)]
layers = [3, 1, -1]
atoms = Decahedron('Rh', 1, 1, 1)
atoms.set_cell([10] * 3)
add_to_dict(AseAtomsAdaptor.get_structure(atoms), 'Rh_cluster')
atoms.center()
view(atoms)
"""
# Si nano-dot

with open('Si29H36.ion','r') as f:
    txt = f.read()

Si = txt.split('Si')[1]
Si, H = Si.split('H')

atoms = Atoms()

for element, coord_set in zip(['Si', 'H'], [Si, H]):
    for line in coord_set.splitlines()[1:-1]:
        atoms += Atom(symbol=element,
                      position=[float(a) * Bohr for a in line.split()])
atoms.set_cell([21] * 3)
atoms.center()
add_to_dict(AseAtomsAdaptor.get_structure(atoms), 'Si_nanodot')


# PV=nRT @ 100 bar, 100 C
atoms = make_box_of_molecules(Atoms('He'), 75, [33.811] * 3)
add_to_dict(AseAtomsAdaptor.get_structure(atoms), 'He_box')

# water

# liquid water at room temp, 55.555 mol/L

water = molecule('H2O')
water.set_distance(0, 1, 1.0, fix=0)
water.set_distance(0, 2, 1.0, fix=0)
atoms = make_box_of_molecules(water, 113, [15] * 3)

#add_to_dict(AseAtomsAdaptor.get_structure(atoms), 'water_box')

# CO2

atoms = make_box_of_molecules(molecule('CO2'), 22, [20] * 3)
add_to_dict(AseAtomsAdaptor.get_structure(atoms), 'CO2_box')


# natural gas


#CH4 = pubchem_atoms_search('methane')
#C2H6 = pubchem_atoms_search('ethane')
#C3H8 = pubchem_atoms_search('propane')


#atoms = make_box_of_molecules([CH4, C2H6, C3H8], [18, 6, 6], [15] * 3)
LNG = json.load(open('LNG.json', 'r'))
add_to_dict(Structure.from_dict(LNG), 'LNG_box')

lens = []
for n, s in structures.items():
    lens.append(len(Structure.from_dict(s)))
    #print(len(Structure.from_dict(s)))

#print(([a > 30 for a in lens]).count(True))
#print(([a > 100 for a in lens]).count(True))
#print(len(structures))

print(len(structures))
json.dump(structures, open('../../test_set/molecular.json', 'w'))
json.dump(params_dict, open('../../test_set/molecular_parameters.json', 'w'))
