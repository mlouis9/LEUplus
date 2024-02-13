"""
A script created to generate properly formatted cross section files, input them
to a nodal core simulator: DYN3D, run DYN3D, and then evaluate the results in
a way that can be sufficiently arbitrary to allow for optimization of core
design by an optimization algorithm (e.g. simulated annealing)

Author: Matthew Louis
Email: matthewlouis31@gmail.com

"""
from serpentTools.settings import rc
from xsInterface.functions.main import Main
from xsInterface.functions.dyn3d import DYN3D
from pathlib import Path
import yaml
from string import Template
rc['serpentVersion'] = '2.1.32'

# ---------------------------------
#    Setup xs-interface files
# ---------------------------------

# First load the unique assembly keys from the adjacent .yaml file
keysPath = str(Path('./keys.yaml'))
with open(keysPath, 'r') as f:
    assyKeys = yaml.safe_load(f)

# Now make a mapping of assembly types in the serpent database to the unique assemblies in the core model
initialMapping = {}
for assyKey in assyKeys['keys']:
    initialMapping.update({assyKey: 'assyw0i0e3-000'})

# Read the univs file template
univs_template_path = Path('./inputs/univ_files/univs_template')
with open(str(univs_template_path), 'r') as f:
    templateLines = f.read()

# Now set the file names in the corresponding univs files
for assyKey, assyType in initialMapping.items():
    with open(str( univs_template_path.parent / f"univs{assyKey}" ), 'w') as f:
        f.write(Template(templateLines).substitute(assemblyType = assyType))


# Now read controlDict and write XS
controlDict = str(Path("./inputs/controlDict")) # --> see definition of controlDict in separate file
xs = Main(controlDict)  # Read control dict
xs.Read(readTemplate=True) # Read xs data and templates and populate data
xs.Write() # Write data to txt files


# casedir = str(Path("./dyn3d"))
# casefile = "core2D"
# exefile = "RUN_DYN3D.sh"

# reslt = DYN3D(xs, casedir, casefile, exefile)
              
# states = {}
# xs.PopulateCoreData(attributes=None,
#                     states=states,
#                     volManip={'infflx': 'divide'},
#                     sph=None, adf=None)

# reslt.Execute()