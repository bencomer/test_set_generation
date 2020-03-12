import json

relaxations_set = {}
relaxations_set_params = {}


surfaces = json.load(open('../../test_set/surface.json', 'r'))
surfaces_params = json.load(open('../../test_set/surface_parameters.json', 'r'))

bulks = json.load(open('../../test_set/bulk.json', 'r'))
bulks_params = json.load(open('../../test_set/bulk_parameters.json', 'r'))

wires = json.load(open('../../test_set/wire.json', 'r'))
wires_params = json.load(open('../../test_set/wire_parameters.json', 'r'))

molecular = json.load(open('../../test_set/molecular.json', 'r'))
molecular_params = json.load(open('../../test_set/molecular_parameters.json', 'r'))

MD = json.load(open('../../test_set/MD.json', 'r'))
MD_params = json.load(open('../../test_set/MD_parameters.json', 'r'))

for entry, value in MD_params.items():
    del value['sparc']['MD_FLAG'] 
    del value['sparc']['ION_TEMP']

    del value['espresso']['tempw']
    del value['espresso']['ion_dynamics']
    value['espresso']['calculation'] = 'relax'

    del value['abinit']['ionmov'] 
    del value['abinit']['mdtemp'] 
    MD_params[entry] = value

for entry, value in surfaces_params.items():
    value['espresso']['dipole']['status'] = True
    MD_params[entry] = value


for mol in ['H2', 'I2', 'HCN', 'H3C2']:
    relaxations_set[mol + '_molecule'] = molecular[mol + '_molecule']
    relaxations_set_params[mol + '_molecule'] = molecular_params[mol + '_molecule']
    relaxations_set_params[mol + '_molecule']['sparc']['RELAX_FLAG'] = 1
    relaxations_set_params[mol + '_molecule']['sparc']['TOL_RELAX'] = 2.5e-3
    relaxations_set_params[mol + '_molecule']['abinit']['ionmov'] = 2
    relaxations_set_params[mol + '_molecule']['abinit']['ntime'] = 300
    relaxations_set_params[mol + '_molecule']['espresso']['forc_conv_thr'] = 2.5e-3
    relaxations_set_params[mol + '_molecule']['espresso']['calculation'] = 'relax'
    relaxations_set_params[mol + '_molecule']['espresso']['ion_dynamics'] = 'bfgs'



for bulk in ['WC', 'RhS2', 'CaO']:
    relaxations_set[bulk + '_bulk'] = bulks[bulk + '_bulk']
    relaxations_set_params[bulk + '_bulk'] = bulks_params[bulk + '_bulk']
    relaxations_set_params[bulk + '_bulk']['sparc']['RELAX_FLAG'] = 3
    relaxations_set_params[bulk + '_bulk']['sparc']['TOL_RELAX'] = 2.5e-3
    relaxations_set_params[bulk + '_bulk']['abinit']['optcell'] = 2
    relaxations_set_params[bulk + '_bulk']['abinit']['ionmov'] = 2
    relaxations_set_params[bulk + '_bulk']['abinit']['ecutsm'] = 1e-8
    relaxations_set_params[bulk + '_bulk']['abinit']['ntime'] = 300
    relaxations_set_params[bulk + '_bulk']['espresso']['forc_conv_thr'] = 2.5e-3
    relaxations_set_params[bulk + '_bulk']['espresso']['calculation'] = 'vc-relax'
    relaxations_set_params[bulk + '_bulk']['espresso']['ion_dynamics'] = 'bfgs'


for surface in ['Zn_100', 'MnP_111', 'Si_111']:
    relaxations_set[surface + '_surface'] = surfaces[surface + '_surface']
    relaxations_set_params[surface + '_surface'] = surfaces_params[surface + '_surface']
    relaxations_set_params[surface + '_surface']['sparc']['RELAX_FLAG'] = 1
    relaxations_set_params[surface + '_surface']['sparc']['TOL_RELAX'] = 2.5e-3
    relaxations_set_params[surface + '_surface']['abinit']['ionmov'] = 2
    relaxations_set_params[surface + '_surface']['abinit']['ntime'] = 300
    relaxations_set_params[surface + '_surface']['espresso']['forc_conv_thr'] = 2.5e-3
    relaxations_set_params[surface + '_surface']['espresso']['calculation'] = 'relax'
    relaxations_set_params[surface + '_surface']['espresso']['ion_dynamics'] = 'bfgs'


json.dump(relaxations_set, open('../../test_set/relaxation.json', 'w'))
json.dump(relaxations_set_params, open('../../test_set/relaxation_parameters.json', 'w'))
