# ***********************************************
# Utility script for generating input files for
# 2D fuel assembly XS generation for varying
# enrichments and BP loadings (IFBA and WABA)
# and creating scripts for running the inputs on
# the GT PACE cluster
#
# Written by: Matthew Louis
#
# ***********************************************

import os, sys
import yaml
import numpy as np
from itertools import product, combinations
import copy
from string import Template
import math
from pathlib import Path

# Template file to append material and lattice definitions to
serpentTemplate = str(Path('../templateFile/serpentTemplateFA.txt'))
outputDir = str(Path('../database/no_pert_only_burnup'))

# Absolute path of script, to avoid creating directories in unwanted locations
scriptsDir = os.path.dirname(os.path.abspath(sys.argv[0]))
serpentTemplate = os.path.join(scriptsDir, serpentTemplate)
outputDir = os.path.join(scriptsDir, outputDir)

# Create output directory if it doesn't exist
if not os.path.isdir(outputDir):
    os.mkdir(outputDir)


# -------------------------------------------------
# Calculate the Burnable Poison Loading Patterns
#--------------------------------------------------


# Read assembly parameters
with open(str(Path(f'{scriptsDir}/assemblyParameters.yaml')), 'r') as f:
    assemblyData = yaml.safe_load(f)

identifiers = assemblyData['identifiers']


def pattern_string_to_array(patternStr: str) -> np.chararray:
    """Utility for converting assembly pattern string to
    an array that can be manipulated to create all combinations of BP loadings
    
    Parameters
    ----------
        A string representing a given burnable poison loading

    Returns
    -------
        A np chararray containing the pin types for the entire assembly
    """
    return np.array([ line.split() for line in patternStr.split('\n') if len(line) != 0 ])


def pattern_array_to_string(inputPatternArray: np.chararray) -> str:
    """Utility for converting assembly pattern array to
    a string that can be appended to a Serpent input file
    
    Parameters
    ----------
        A np chararray containing the pin types for the entire assembly

    Returns
    -------
        A string representing a given burnable poison loading
        
    """
    patternArray = inputPatternArray.copy()

    # First find the longest character, so that each pin character can be equally spaced
    length_checker = np.vectorize(len)
    charLen = np.max(length_checker(patternArray))
    
    # Now pad pin characters with spaces
    def space_padder(string):
        return string.ljust(charLen)

    array_padder = np.vectorize(space_padder)
    if np.size(patternArray[np.where(length_checker(patternArray) < charLen)]) != 0: # Results in strange behavior with array_padder
        patternArray[np.where(length_checker(patternArray) < charLen)] = array_padder(patternArray[np.where(length_checker(patternArray) < charLen)])

    strPattern = ''
    for row in patternArray.tolist():
        strPattern += ' '.join(row) + '\n'

    return strPattern


# Create arrays of individual burnable poison (BP) types (that will later be combined to get all possible BP loadings)
assemblies = dict()

for key, identifier in identifiers.items():
    # Get only assembly dicts of a given type
    assemblies[identifier['name']] = [ assembly for assembly in assemblyData['assemblies'] \
                                       if assembly['type'] == identifier['name'] ]

    # convert assembly patterns to np chararrays
    for assembly in assemblies[identifier['name']]:
        assembly['assembly_pattern'] = pattern_string_to_array(assembly['assembly_pattern'])

# Now, generate all combinations of BP loadings
all_assembly_combinations = []
bpTypes = list(assemblies.keys())

for assembly_combination in product(*[assemblies[bpType] for bpType in bpTypes]):
    # replace fuel/gt locations in first assembly with bps of the other two assemblies
    bp_indices = [ assembly['assembly_pattern'] == identifiers[assembly['type']]['symbol'] for assembly in assembly_combination ]

    # Raise error if the bp configurations overlap (can deal with this by trying different orderings, but not neecessary right now)
    # we do this by checking all pairwise combinations of the bp indices and checking that none have like entries

    pairwise_combinations = list(combinations(bp_indices, 2))
    if np.any([ np.logical_and(pair[0], pair[1]) for pair in pairwise_combinations ]):
        raise Exception("""Burnable poison designs overlap, i.e. there are patterns with identical burnable poison locations - this is not currently supported""")

    # The combination assembly will be written to the pattern of the first assembly in the combination
    # Note: deepcopy is needed to recursively copy dictionary objects
    combined_assembly = copy.deepcopy(assembly_combination[0])
    for index, assembly in enumerate(assembly_combination):
        combined_assembly['number_' + assembly['type']] = assembly['number']
        if index != 0: # Spkipping "combined assembly"
            combined_assembly['assembly_pattern'][bp_indices[index]] = identifiers[assembly['type']]['symbol']

            del combined_assembly['type']
            del combined_assembly['number']
    
    # Now, convert chararray to text
    combined_assembly['assembly_pattern'] = pattern_array_to_string(combined_assembly['assembly_pattern'])
    all_assembly_combinations.append(combined_assembly)


# --------------------------------------------------------
# Calculate the Fuel Material as a function of Enrichment
#---------------------------------------------------------

# Create array of unique enrichments
low_enrichment = assemblyData['enrichments']['start']
high_enrichment = assemblyData['enrichments']['stop']
num_enrichments = assemblyData['enrichments']['number']

enrichments = np.linspace(low_enrichment, high_enrichment, num_enrichments)

# Define base parameters for enrichment calculation
thDens = 10.96 # g/cm3
fracTh = 0.95  # I.e. 95% theoretical density
dens = fracTh*thDens
A238 = 238.05078
A235 = 235.043923
afU = 1/3   # i.e. UO2 is 1/3 atom fraction U 
afO = 2/3   # i.e. UO2 is 2/3 atom fraction O

# Define fuel material template (note assuming only O-16)
fuelMatTemplate = \
"""mat fuel -{:1.8E} burn 1
U-235.09c   {:1.8E}
U-238.09c   {:1.8E}
 O-16.09c   0.66666667E+00
"""

def enrichment_calculator(e: float)-> str:
    """This is a utility function for calculating the weight fractions of U-235 and
    U-238 for defining a fuel material in Serpent

    Parameters
    ----------
        e: The enrichment (as weight fraction of U-235 in U)
    
    Returns
    -------
        A string describing the fuel material
    """

    # Calculate the atomic fraction of U-235 in U
    wf235 = e
    af235 = A238*wf235/( (1-wf235)*A235 + wf235*A238 )

    # Calculate the total atomic fractions of U-235 and U-238 in the fuel
    af235tot = afU*af235
    af238tot = afU*(1-af235)

    return fuelMatTemplate.format(dens, af235tot, af238tot)


# ------------------------------------
# Generate the Input Files
# ------------------------------------

latticeHeader = "lat 1 1 0.0 0.0 17 17 1.25984\n"
fileNamingTemplate = "assyw{}i{}e{:.3f}"

# Read template file
with open(serpentTemplate, 'r') as f:
    templateLines = f.read()

for assembly in all_assembly_combinations:
    for enrichment in enrichments:
        number_waba = int(assembly['number_waba'])
        number_ifba = int(assembly['number_ifba'])
        outfileName = fileNamingTemplate.format(number_waba, number_ifba, enrichment*100).replace('.', '-')
        with open(str(Path(f"{outputDir}/{outfileName}")), 'w') as f:
            f.write(Template(templateLines).substitute(root_universe_name = '0')) # Just setting root universe to '0' for now

            # Add fuel material
            f.write(f"% {enrichment*100:.3f}w/o Enriched Fuel\n")
            f.write(enrichment_calculator(enrichment))
            
            # Add assembly lattice
            f.write("\n % Assembly Loading Pattern\n")
            f.write(latticeHeader)
            f.write(assembly['assembly_pattern'])

            # Set case title
            f.write(f"\nset title \"2D 17x17 PWR Assembly XS gen. {number_waba} waba, {number_ifba} ifba, {enrichment*100:.3f} w/o enrichment \"\n")

# ------------------------------
# Now write a mpi4py run script
# ------------------------------            

# Define parameters for running serpent cases on PACE
coresPerAssy = 4
numberOfAssyRuns = len(all_assembly_combinations)*len(enrichments)
coresPerNode = 24
tasksPerNode = math.floor(coresPerNode/coresPerAssy)
numberOfNodes = math.ceil(numberOfAssyRuns/tasksPerNode)
numberOfCores = numberOfAssyRuns*coresPerAssy

# This is a template for creating an MPI python script that runs all of the generated input files
runScriptTemplate = Template("""#-------------------------------
# Script for running generated
# input files on PACE
#-------------------------------

from mpi4py import MPI
import subprocess
import numpy as np
from uniqueAssemblies import assemblies, enrichments
                             
# Get enrichments and assemblies for accessing 
numberOfAssyRuns = len(assemblies)*len(enrichments)

# Make an array of enrichment/assembly pairs so that each MPI rank may be assigned a unique run
assembliesIndices = [(i, j) for i in range(len(assemblies)) for j in range(len(enrichments))]

comm = MPI.COMM_WORLD

# Get parameters of assembly
assemblyIndices = assembliesIndices[comm.rank]
assembly = assemblies[assemblyIndices[0]]
enrichment = enrichments[assemblyIndices[1]]
number_waba = assembly['number_waba']
number_ifba = assembly['number_ifba']
                             
fileNamingTemplate = \"assyw{}i{}e{:.3f}\"

# Format Serpent run command and run the case
serpentCommand = f\"sss2 {fileNamingTemplate.format(number_waba, number_ifba, enrichment*100).replace('.', '-')}\" + \" -omp $coresPerAssy\"
subprocess.run(serpentCommand, shell=True)
""")

# This is a template for storing the unique assembly data (for the generated database) that can easily be read by a python script
# e.g. the script that processes the database and runs fuel cycle optimization
uniqueAssembliesTemplate = Template("""import numpy as np
enrichments = np.array($enrichments)
assemblies = $assemblies""")

# Force numpy array to print on single line
np.set_printoptions(linewidth=np.inf)

# Now pop off the assembly_pattern element of each dictionary because it is no longer needed
for assembly in all_assembly_combinations:
    assembly.pop('assembly_pattern')

uniqueAssembliesName = 'uniqueAssemblies.py'
with open(str(Path(f'{outputDir}/{uniqueAssembliesName}')), 'w') as f:
    f.write(uniqueAssembliesTemplate.substitute(enrichments = np.array2string(enrichments, separator=','),
                                                assemblies = all_assembly_combinations))

# Now format template and write run script
runScriptName = 'run.py'
with open(str(Path(f'{outputDir}/{runScriptName}')), 'w') as f:
    f.write(runScriptTemplate.substitute(coresPerAssy = coresPerAssy))


# This is a slurm submission script template for running the cases using the above run script
# Note must fill in charge code by hand
slurmSubmissionTemplate = Template("""#!/bin/bash

#SBATCH -J serpentXSGen                         # Job name
#SBATCH -A $chargeCode                          # account to which job is charged, ex: GT-gburdell3
#SBATCH -N$numNodes --ntasks-per-node=$tasksPerNode                 # Number of nodes and cores per node
#SBATCH --cpus-per-task=$cpusPerTask
#SBATCH --time=7:30:00                          # Walltime = 7.5h
#SBATCH -qembers                                # QOS Name (embers / inferno)
#SBATCH -oReport-%j.out                         # Combined output and error messages file

echo $$SLURM_SUBMIT_DIR
cd   $$SLURM_SUBMIT_DIR                            # change into directory from where job to be executed (where data is / script is located)

source ~/.bashrc                                # load environment
echo "Started on `/bin/hostname`"               # prints the name of the node job started on
srun python $runScriptName
""")

# Now write the slurm submission script
submissionScriptName = 'run.slurm'
chargeCode = 'gts-dkotlyar6-CODA20'
with open(str(Path(f'{outputDir}/{submissionScriptName}')), 'w') as f:
    f.write(slurmSubmissionTemplate.substitute(numNodes = numberOfNodes,
                                               tasksPerNode = tasksPerNode,
                                               cpusPerTask = coresPerAssy,
                                               runScriptName = runScriptName,
                                               totNumTasks = numberOfAssyRuns,
                                               chargeCode = chargeCode))