# README #

This is repository for coupling PFLOTRAN into CLM in DOE sponsored NGEE-Arctic Project. It's keeping updates from PLOTRAN-Dev as soon as possible.

The model coupling aims to provide a full alternative solution for CLM-CN's surface-subsurface C/N biogeochemistry and thermal-hydrology, i.e. PFLOTRAN.

This repository is the PFLOTRAN portion for coupling. The CLM portion is in CLM-PFLOTRAN-NGEE-SCI repository.

### Branches and Versions ###

* The repository contains 3 branches for development and stable release purposes:
  
(1) default - the most updated PFLOTRAN codes with necessary editing for coupling. NOT suggested to use for production run. updated: 2015-05-15

STATUS: UNDER development.

(2) default-bgc - the current develop version for coupling subsurface BGC between CLM and PFLOTRAN (both soil thermal-hydrology are from CLM45, so NO transport BUT can have diffusion). 

Current Version of CLM: 4.5.38; 
4.5_r85 (#define FLEXIBLE_POOLS in Line 2 of clm_pflotran_interfaceMod.F90).
 
PFLOTRAN-dev change set: 7dc944a171993b2b27c4a5de84ac0856404a5194 [7dc944a17199] 

STATUS: UNDER development for flexible decomposing C/N pools.

(3) default-hc - the current stable version for coupling subsurface BGC/HYDROLOGY (reactive-transport) between CLM and PFLOTRAN (including Richards only hydrology). 

Current Version of CLM: 4.5.38; 
4.5_r85 (!#define FLEXIBLE_POOLS in Line 2 of clm_pflotran_interfaceMod.F90).

PFLOTRAN-dev change set: 7dc944a171993b2b27c4a5de84ac0856404a5194 [7dc944a17199]

STATUS: UNDER test of Richards mode; stable for BGC coupling. (81bfd59d8396fd45272856fd96d24791e59b4966 [81bfd59d8396] @July-07-2015).  

(4) default-thc - the current development version of clm-pflotran, specifically for testing all thermal-hydrology and biogeochemistry portion. 

Target CLM Version: CLM4_5_1_r85. (2015-03-03). 

PFLOTRAN-dev change set: 87da2df5ab30a429a3816854a304d7f4823d765f [87da2df5ab30] (2015-05-15)

### How do I get set up? ###

(1) hg clone the repository. 

(2) hg update branch to 'default-hc'.

(3) if not coupled with CLM, this repository should be a stand-alone PFLOTRAN model (the Source code directory: ./Src/pflotran, and build it by issuing command: make PETSC_DIR=$PETSC_DIR pflotran, where $PETSC_DIR is your PETSC_DIR directory).

(4) if coupled with PFLOTRAN, first build a libpflotran.a in ./Src/clm-pflotran (First, run the script link_files.sh to copy PFLOTRAN codes, and then build the library: make PETSC_DIR=$PETSC_DIR libpflotran.a). Then build CLM with this library. 

UPDATED: 2015-08-10