# README #

This is repository for coupling PFLOTRAN into CLM in DOE sponsored NGEE-Arctic Project. It's keeping updates from PLOTRAN-Dev as soon as possible.

The model coupling aims to provide a full alternative solution for CLM-CN's surface-subsurface C/N biogeochemistry and thermal-hydrology, i.e. PFLOTRAN.

This repository is the PFLOTRAN portion for coupling. The CLM portion is in CLM-PFLOTRAN-NGEE-SCI repository.

### Branches and Versions ###

* The repository contains 3 branches for development and stable release purposes:
  
***(1) default*** - the most updated PFLOTRAN codes (change set:c988fe5264f6dd47a5b410f64306e1d071282336 [c988fe5264f6]) with CLM-PFLOTRAN grid/mesh coupling. NOT suggested to use for production run. updated: 2016-08-05

***(2) default-bgc*** - the current develop version for coupling subsurface BGC between CLM and PFLOTRAN (both soil thermal-hydrology are from CLM45, so NO transport BUT can have diffusion). 

Current Version of CLM: 4.5.35; 

***(3) default-thc*** - the current development version of clm-pflotran, specifically for testing all thermal-hydrology and biogeochemistry portion. 

Target CLM Version: CLM4_5_1_r85. (2015-08-12).
Stable CLM Version: CLM4_5_35. (2016-05-19). 

STATUS: UNDER test for TH mode; STABLE for Richards+BGC. (2016-08-05) 

### How do I get set up? ###

(1) hg clone the repository. 

(2) hg update branch to 'default-hc'.

(3) if not coupled with CLM, this repository should be a stand-alone PFLOTRAN model (the Source code directory: ./Src/pflotran, and build it by issuing command: make PETSC_DIR=$PETSC_DIR pflotran, where $PETSC_DIR is your PETSC_DIR directory).

(4) if coupled with PFLOTRAN, first build a libpflotran.a in ./Src/clm-pflotran - 

First, run the script link_files.sh to copy PFLOTRAN codes

Secondly, build the library: make PETSC_DIR=$PETSC_DIR use_characteristic_curves=TRUE libpflotran.a

Then build CLM with this library. 

UPDATED: 2016-05-27