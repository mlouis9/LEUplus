**Table of Contents**:
---
- [Introduction](#introduction)
- [Serpent\_Input\_Generator.py](#serpent_input_generatorpy)
- [assemblyParameters.yaml](#assemblyparametersyaml)
- [Notes](#notes)
---

# Introduction
This directory contains the scripts used for generating various cross section databases.

# Serpent_Input_Generator.py
This is a script that is used to generate the Serpent input files (for cross section generation) from [the template file](../templateFile/serpentTemplateFA.txt) for varying fuel enrichments and burnable poison loadings.
- This script also generates corresponding run scripts that will automatically parallelize the cross section generation using MPI via a python script `run.py`.
- Once the database is generated, to begin the calculation, just run `sbatch run.slurm`.
- For the slurm run script (`run.slurm`) to be templated with a valid charge code, it is necessary to include a text file `chargeCode.txt` with your charge code.

# assemblyParameters.yaml
This file contains the parameters specifying how many enrichment/burnable poison combinations need to be calculated in the database. 
- For each burnable poison loading, there is a corresponding loading pattern
- In this file, `assemblies` describes the unique burnable poison loadings, i.e. for a specific burnable poison type (e.g. IFBA) and a specific number of burnable poison rods, what the loading is.
- The database is generated by taking all possible combinations of these assembly loading patterns
- Currently, the individual burnable poison loading patterns for IFBA and WABA are _not_ allowed to overlap. It's nothing to worry about, however, because this will raise an error in the `SerpentInputGeneratorFA2D.py` script.

# Notes
- By default the submission script is templated to run on `qinferno`, the production submission QOS on pace's phoenix cluster (more info can be found [here](https://gatech.service-now.com/home?id=kb_article_view&sysparm_article=KB0041966)); if you'd like to change this, just modify the submission script template in [the input geneation script](SerpentInputGeneratorFA2D.py).