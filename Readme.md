# README #

This is repository for coupling PFLOTRAN into ELM in DOE sponsored NGEE-Arctic Project. It's keeping updates from PLOTRAN-Dev as soon as possible.

The model coupling aims to provide a full alternative solution for CLM-CN's surface-subsurface C/N biogeochemistry and thermal-hydrology, i.e. PFLOTRAN.

This repository is the PFLOTRAN portion for coupling. The ELM portion is in E3SM repository, branch /fmyuan/lnd/elm-pflotran.

***UPDATE: (2020-03-26) The coupling with E3SM Land Model (ELM) now is working.***
(E3SM Branches tested: **fmyuan/lnd/elm-pflotran**, or,  **fmyuan/lnd/HumHol_v2**)

### Branches and Versions ###

* The repository contains master and FIVE(5) branches for different versions of PETSc:

***master*** 

 - the most updated PFLOTRAN codes forked from https://bitbucket.org/pflotran/pflotran
 - only a few minor changes. 
 - NOT suggested to use coupling run, rather as a stand-alone PFLOTRAN. 
 - updated: **2020-02-12**

***(1) pflotran-elm-interface*** (**Default**)

 - the current development version of elm-pflotran, specifically for testing biogeochemistry (BGC) portion. 
 - **PFLOTRAN**: It's keeping updating with recent master (f08b131141ddb91dcfe902d37fa49a1dcf0844dc) on **11/13/2019**, and with the most recent ***PETSc-dev***

***(2) rebase2020monthly***
- rebased pflotran-elm-interface with master, but not yet fully tested (***Feb 2020***).

***(3) default-release-v3.8*** and above
- for maintenance purpose, simply saved past versions, which mostly consistent with PETSc version. 
- for BGC coupled PFLOTRAN, with PETSc 3.8.x and above. TH coupling is not available.
- It's working with **ELM v1.1 above**. May not be working with CLM.

***(4) default-release-v3.7***
 - The current STABLE version of clm-pflotran, specifically for testing all thermal-hydrology (TH) and biogeochemistry (BGC) portion. 
 - **PFLOTRAN**: 6571d9ae3d1cdd699bfc993155c958173e76c20f [6571d9a]. It's much older (**Updated: 2017-03-09**) than 'master', and with ***PETSc-dev*** version **3.7.x**.
 - Stable CLM Version: **CLM4_5_35**. (2017-04-25).  *STATUS*: TH & BGC coupling is STABLE. 
 - Target ELM Version: **ELM v1**.  *STATUS*: STABLE for ```BGC```. (2017-06-27); UNDER test for ```TH mode``` 

  

### How do I get set up? ###

**(1)** *git clone the repository*
```
git clone https://github.com/fmyuan/pflotran-elm-interface.git
```

**(2)[OPTIONAL]** *git checkout your specific branch, e.g. 'default-release-v3.7'. The default is 'pflotran-elm-interface'*
```
cd pflotran-interface
git checkout default-release-v3.7
```

**(3a)** *if not coupled with ELM, this repository should be a **stand-alone PFLOTRAN** model.*

- the Source code directory: 
```
cd ./src/pflotran
```

- build it by issuing command:
```
make PETSC_DIR=$PETSC_DIR pflotran
(where $PETSC_DIR is your PETSC_DIR directory)
```

- check your repository (regression test):
```
make PETSC_DIR=$PETSC_DIR test
(where $PETSC_DIR is your PETSC_DIR directory)
```

**(3b)** *if coupling ELM with PFLOTRAN, need to build a library named as **libpflotran.a**.*
```
cd ./src/pflotran-elm-interface
```

***FIRST***, run specific make script to softlink PFLOTRAN source code files (*.F90).
```
make link_common_src

(OR, make clean_common_src 
to unlink PFLOTRAN source code files, and only leave PFLOTRAN-ELM interface codes)
```

***SECONDLY***, build the library (vertically-only mode ON)
```
make PETSC_DIR=$PETSC_DIR column_mode=TRUE libpflotran.a

(OR, make PETSC_DIR=$PETSC_DIR column_mode=TRUE debugbuild=TRUE libpflotran.a
for a library with '-g -O0')

```

***FINALLY***, build CLM (ELM) with this library.

*I.* **Macro (CLM)** or **Macro.make (ELM)** or **Macro.cmake (ELM master since@2019-07)** modified, BY any means, to include $PFLTRAN_INC and to link $PFLOTRAN_LIB -

```
(Macro/Macro.make)


  ifeq ($(MODEL),clm)
    FFLAGS := $(FFLAGS)  $(PFLOTRAN_INC)
  endif
  
  ......
  
  ifeq ($(MODEL),driver)
    ......
    LDFLAGS := $(LDFLAGS)  $(PFLOTRAN_LIB)
  endif

```
**Macro.cmake: *ELM master since@2019-07)***
```
  if("${MODEL}" STREQUAL "clm")
    set(FFLAGS "${FFLAGS}  $(PFLOTRAN_INC)")
  endif()

  ......
  
  if("${MODEL}" STREQUAL "driver")
    ......
    set(LDFLAGS "${LDFLAGS}  $(PFLOTRAN_LIB)")
  endif()

```
**IF don't want to modify Macro.make Macro.cmake** after case.setup as above, please add those 2 into **config_compliers.xml** as following:
```
<FFLAGS>
   ......

  <!-- A NOTE here: $(PFLOTRAN_INC) may be empty, or both PETSC and actual PFLOTRAN include dir, or only PETSC include dir -->
  <append MODEL="clm"> $(PFLOTRAN_INC) </append>
</FFLAGS>

......

<LDFLAGS>
  ......
  <!-- A NOTE here: $(PFLOTRAN_LIB) may be empty, or both PETSC libraries and actual PFLOTRAN libray, or only PETSC libraries -->
  <append MODEL="driver"> $(PFLOTRAN_LIB) </append>
</LDFLAGS>
```


*II.* **config_machines.xml** editing for each supported machine. *OR*: after './case.setup', edit **'env_mach_specific.xml'** to add (OR delete, if don't want to build ELM with pflotran codes.  

*NOTE*: 
PETSC_PATH (or PETSC_DIR), CLM_PFLOTRAN_SOURCE_DIR, can be full paths OR pre-defined in your '.bashrc' or '.bash_profile' in your HOME, as **user defined as in your environmental setting**.
e.g. on ORNL CADES, 
```
export PETSC_DIR=$CCSI_USERTOOLS/petsc-x/openmpi-1.10-gcc-5.3
export PETSC_ARCH=arch-orcondo-openmpi-gcc53-nodebug
export CLM_PFLOTRAN_SOURCE_DIR=/lustre/or-hydra/cades-ccsi/proj-shared/models/pflotran-dev/src/pflotran-elm-interface

```


```
<!-- hack for PFLOTRAN coupling to build ELM model.
     this is a temporary solution, and user must manually edit the following after cloning codes.
     By default, model will build as PFLOTRAN coupled and run as non-coupled.
-->
<environment_variables>
  <!-- for CLM-PFLOTRAN coupling, the $PETSC_DIR and $CLM_PFLOTRAN_SOURCE_DIR must be user-defined specifically upon machines, e.g. in .bashrc -->
  <!-- pflotran codes are pre-built as libpflotran.a in $CLM_PFLOTRAN_SOURCE_DIR -->
  <env name="PETSC_PATH">$ENV{PETSC_DIR}</env>
  <env name="PFLOTRAN_INC"> -I$ENV{CLM_PFLOTRAN_SOURCE_DIR} -I$ENV{PETSC_DIR}/include</env>
  <env name="PFLOTRAN_LIB"> -L$ENV{CLM_PFLOTRAN_SOURCE_DIR} -lpflotran -L$ENV{PETSC_DIR}/lib -lpetsc -lparmetis -lmetis</env>
</environment_variables>

```



*III.* Specifically, for **ELM** build with PFLOTRAN, as ***external module***.
(e.g. https://github.com/fmyuan/E3SM.git, branch 'elm-pflotran-II')
In this way, ELM source codes will have copied all PFLOTRAN source codes ONLY into:
```
      $SRCROOT/components/clm/src/external_models/pflotran-interface/src
```
AND, add the following line in ELM building config file ($SRCROOT/components/clm/bld/configure):
```
       "external_models/sbetr/src/Applications/soil-farm/CENT_ECACNP",
+	     "external_models/pflotran-interface/src/clm-pflotran",
		     "utils", 
		     "cpl" );
```

If setting software environments (specifically PETSc library), the whole ELM building process will automatically setup source codes and dependecies for coupled component, and e3sm.exe usable for either standalone or coupling simulation.



***UPDATED: 2020-03-26***
