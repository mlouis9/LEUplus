# 2D Core With No Burnable Poisons

This directory contains an example of a full-core calculation with DYN3D where the following assumptions have been made
* No burnable poisons, i.e. only considering the effect of assemblies with different enrichments
* Calculations neglecting thermal hydraulics
* A 2D core (i.e. no axial leakage)
* Fixed shuffling scheme

# Files and Directory Structure
Adjacent to this readme, there are two files:
1.  `genXS.ipynb`: An interactive notebook used for running the DYN3D calculations and viewing the results using xs-interface
2.  `coreData.yaml`: A file containing the unique assembly identifiers, a core map, and an index map
    * For a given core, not all assemblies are unique, and sense cross sections are generated only for unique assemblies it makes sense to consider only those
    * The core map specifies how each of the (25) unique assemblies are laid out in the core
    * The index map gives the indices of each of the core assembly locations (used for indexing within xs-interface)

## DYN3D
The DYN3D folder contains the the following DYN3D input files
1. ``core2D_kin.dat``: Input file for overall calculation and neutronics
2. ``core2D_thy.dat``: Input file for thermal hydraulics
3. ``core2D_wqs.dat``: Input file for the cross section library (which is contained in the `xs` directory)

The rest of the files are binary/text result files that are automatically read with xs-interface.

## Serpent
The serpent directory is a symbolic link to the database containing only enrichment variations, and contains all of the serpent _res.m files needed for generating cross sections

## Inputs
This directory contains the input files for  xs-interface.
* `template`: The file which describes the template for xs data. This is consistent with DYN3D, and the templated xs files can be found in the `dyn3d/xs` directory.
* `controlDictTemplate`: This file is a template for the control dict. This is templated so that the core layout can be changed on the fly by the run script by just templating the volume mapping and channel mapping.
    * `controlDict` is the templated output file that is used by xs-interface
* `univ_files`: Is a directory containing the universe files for each unique assembly, created by templating `univs_template`.
    * Doing this again allows chaning which serpent `_res.m` files correspond to which universes.
    * Note it's currently necessary to create a separate universe file for each unique assembly, because only one output file can be specified in the universe file.