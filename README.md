## A pedigree-based map of crossovers and non-crossovers in aye-ayes (*Daubentonia madagascariensis*)

This repository contains the main workflows generated for identifying recombination events in a pedigree-based aye-aye dataset.

#### Code

The following folders contain the workflows for the three main processes described in `Material&Methods`. Each folder consists of a `Snakefile` (a concatenation of the `Snakemake` rules), a `config` file (including the paths and variables - which should be costumized by the user) and a `.sh` execution file. The scripts used as part of the `Snakefile` are included in the `scripts` folders.

- [Variant filtering](Snakepit/Variant_filtering/)
- [Pedigree approach](Snakepit/Pedigree_approach/)
- [Family approach](Snakepit/Family_approach/)