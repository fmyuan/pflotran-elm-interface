# README #

This is repository for coupling PFLOTRAN into CLM in DOE sponsored NGEE-Arctic Project. It's keeping updates from PLOTRAN-Dev as soon as possible.

The model coupling aims to provide a full alternative solution for CLM-CN's surface-subsurface C/N biogeochemistry and thermal-hydrology, i.e. PFLOTRAN.

This repository is the PFLOTRAN portion for coupling. The CLM portion is in CLM-PFLOTRAN-NGEE-SCI repository.

### Branches and Versions ###

* The repository contains 3 branches for development and stable release purposes:
  
(1) default - the most updated PFLOTRAN codes with necessary editing for coupling. NOT suggested to use for production run. updated: 2015-03-02

(2) default-bgc - the current stable version for coupling subsurface BGC between CLM and PFLOTRAN (including Richards only hydrology). Current Version of CLM: 4.5.38; PFLOTRAN-dev change set: 7dc944a171993b2b27c4a5de84ac0856404a5194 [7dc944a17199] 

(3) default-thc - the current development version of clm-pflotran, specifically for testing all thermal-hydrology and biogeochemistry portion. Target CLM Version: CLM4_5_1_r85. (2015-03-03). PFLOTRAN-dev change set: 7a77c0606ca99efa9b58f6fe9d5ba23a82b1f8a1 [7a77c0606ca9]

### How do I get set up? ###

(1) hg clone the repository. 

(2) hg update branch to 'default-bgc'.

(3) if not coupled with CLM, this repository should be a stand-alone PFLOTRAN model (the Source code directory: ./Src/pflotran).

(4) if coupled with PFLOTRAN, first build a libpflotran.a in ./Src/clm-pflotran (First, run the script link_files.sh to copy PFLOTRAN codes, and then build the library). Then build CLM with this library. 

UPDATED: 2015-03-03