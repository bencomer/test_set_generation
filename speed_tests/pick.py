import json

speed_set = {}
speed_set_params = {}


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
    value['espresso']['calculation'] = 'scf'

    del value['abinit']['ionmov'] 
    del value['abinit']['mdtemp'] 
    MD_params[entry] = value

for entry, value in surfaces_params.items():
    value['espresso']['dipole']['status'] = True
    MD_params[entry] = value


# 20-50

speed_set['TiO2_NH3ads_110_surface'] = surfaces['TiO2_NH3ads_110_surface']
speed_set_params['TiO2_NH3ads_110_surface'] = surfaces_params['TiO2_NH3ads_110_surface'] 


speed_set['MoS2_nanotube'] = wires['MoS2_nanotube']
speed_set_params['MoS2_nanotube'] = wires_params['MoS2_nanotube'] 

# 50-100

speed_set['Si_nanodot_molecule'] = molecular['Si_nanodot_molecule']
speed_set_params['Si_nanodot_molecule'] = molecular_params['Si_nanodot_molecule']

#speed_set['Cu_defect_100_surface'] = surfaces['Cu_defect_100_surface']
#speed_set_params['Cu_defect_100_surface'] = surfaces_params['Cu_defect_100_surface']

speed_set['CO2_box_molecule'] = molecular['CO2_box_molecule']
speed_set_params['CO2_box_molecule'] = molecular_params['CO2_box_molecule']

# 100-300

speed_set['Al_MD'] = MD['Al_MD']
speed_set_params['Al_MD'] = MD_params['Al_MD']

speed_set['Na_bulk'] = bulks['Na_bulk']
speed_set_params['Na_bulk'] = bulks_params['Na_bulk']

speed_set['MgS_bulk'] = bulks['MgS_bulk']
speed_set_params['MgS_bulk'] = bulks_params['MgS_bulk']

# 300

#speed_set['electrolyte_MD'] = MD['electrolyte_MD']
#speed_set_params['electrolyte_MD'] = MD_params['electrolyte_MD']


#speed_set['iron_interface_MD'] = MD['iron_interface_MD']
#speed_set_params['iron_interface_MD'] = MD_params['iron_interface_MD']

json.dump(speed_set, open('../../test_set/speed.json', 'w'))
json.dump(speed_set_params, open('../../test_set/speed_parameters.json', 'w'))
