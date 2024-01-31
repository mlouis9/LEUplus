#-------------------------------
# Script for running generated
# input files on PACE
#-------------------------------

from mpi4py import MPI
import subprocess
import numpy as np
                             
# Get enrichments and assemblies for accessing 
enrichments = np.array([0.03      ,0.03363636,0.03727273,0.04090909,0.04454545,0.04818182,0.05181818,0.05545455,0.05909091,0.06272727,0.06636364,0.07      ])
assemblies = [{'number_ifba': 0, 'number_waba': 0}, {'number_ifba': 0, 'number_waba': 4}, {'number_ifba': 0, 'number_waba': 24}, {'number_ifba': 0, 'number_waba': 20}, {'number_ifba': 0, 'number_waba': 16}, {'number_ifba': 0, 'number_waba': 12}, {'number_ifba': 0, 'number_waba': 8}, {'number_ifba': 156, 'number_waba': 0}, {'number_ifba': 156, 'number_waba': 4}, {'number_ifba': 156, 'number_waba': 24}, {'number_ifba': 156, 'number_waba': 20}, {'number_ifba': 156, 'number_waba': 16}, {'number_ifba': 156, 'number_waba': 12}, {'number_ifba': 156, 'number_waba': 8}, {'number_ifba': 128, 'number_waba': 0}, {'number_ifba': 128, 'number_waba': 4}, {'number_ifba': 128, 'number_waba': 24}, {'number_ifba': 128, 'number_waba': 20}, {'number_ifba': 128, 'number_waba': 16}, {'number_ifba': 128, 'number_waba': 12}, {'number_ifba': 128, 'number_waba': 8}, {'number_ifba': 104, 'number_waba': 0}, {'number_ifba': 104, 'number_waba': 4}, {'number_ifba': 104, 'number_waba': 24}, {'number_ifba': 104, 'number_waba': 20}, {'number_ifba': 104, 'number_waba': 16}, {'number_ifba': 104, 'number_waba': 12}, {'number_ifba': 104, 'number_waba': 8}, {'number_ifba': 80, 'number_waba': 0}, {'number_ifba': 80, 'number_waba': 4}, {'number_ifba': 80, 'number_waba': 24}, {'number_ifba': 80, 'number_waba': 20}, {'number_ifba': 80, 'number_waba': 16}, {'number_ifba': 80, 'number_waba': 12}, {'number_ifba': 80, 'number_waba': 8}, {'number_ifba': 64, 'number_waba': 0}, {'number_ifba': 64, 'number_waba': 4}, {'number_ifba': 64, 'number_waba': 24}, {'number_ifba': 64, 'number_waba': 20}, {'number_ifba': 64, 'number_waba': 16}, {'number_ifba': 64, 'number_waba': 12}, {'number_ifba': 64, 'number_waba': 8}, {'number_ifba': 48, 'number_waba': 0}, {'number_ifba': 48, 'number_waba': 4}, {'number_ifba': 48, 'number_waba': 24}, {'number_ifba': 48, 'number_waba': 20}, {'number_ifba': 48, 'number_waba': 16}, {'number_ifba': 48, 'number_waba': 12}, {'number_ifba': 48, 'number_waba': 8}, {'number_ifba': 32, 'number_waba': 0}, {'number_ifba': 32, 'number_waba': 4}, {'number_ifba': 32, 'number_waba': 24}, {'number_ifba': 32, 'number_waba': 20}, {'number_ifba': 32, 'number_waba': 16}, {'number_ifba': 32, 'number_waba': 12}, {'number_ifba': 32, 'number_waba': 8}, {'number_ifba': 16, 'number_waba': 0}, {'number_ifba': 16, 'number_waba': 4}, {'number_ifba': 16, 'number_waba': 24}, {'number_ifba': 16, 'number_waba': 20}, {'number_ifba': 16, 'number_waba': 16}, {'number_ifba': 16, 'number_waba': 12}, {'number_ifba': 16, 'number_waba': 8}]
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
                             
fileNamingTemplate = "assyw{}i{}e{:.3f}"

# Format Serpent run command and run the case
serpentCommand = f"sss2 {fileNamingTemplate.format(number_waba, number_ifba, enrichment*100).replace('.', '-')}" + " -omp 4"
subprocess.run(serpentCommand, shell=True)
