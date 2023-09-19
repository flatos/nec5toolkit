# Compiling on a Linux System

## Compiling with the Intel Fortran Compiler and Intel Math Kernel Library (MKL)

Install the Intel Fortran Compiler and Intel Math Kernel Library available [here](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html).

For setting command line options correctly, the [oneAPI Math Kernel Library Link Line Advisor](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-link-line-advisor.html) is extremely useful.


### Building NEC5 with the dynamic MKL library  -->  nec5

Compile:
```bash
ifx -c -I${MKLROOT}/include/intel64/lp64 -I\"${MKLROOT}/include\" NECMP_MOD.f NECSEH_MOD.f NECSP_MOD.f
ifx -c -I${MKLROOT}/include/intel64/lp64 -I\"${MKLROOT}/include\" Datagn.f GASYEH.F GASYP.F NecMPCL.f NECMPFLD.f SOMGEH.f SOMGP.f SOMLIB_PEH.F
```

Link:
```bash
ifx -o nec5   ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl Datagn.o  GASYEH.o  GASYP.o  NecMPCL.o  NECMPFLD.o  NECMP_MOD.o  NECSEH_MOD.o  NECSP_MOD.o  SOMGEH.o  SOMGP.o  SOMLIB_PEH.o
```

As of this writing, the file NECMP_MOD.f contains an extraneous parenthesis on the line starting with `SUBROUTINE ALLOC_PARAM_TAB(MXPARAM)` which needs to be edited out.

### Building NEC5 with the static MKL library  -->  nec5s

Compile:
```bash
ifx -c -I${MKLROOT}/include/intel64/lp64 -I\"${MKLROOT}/include\" NECMP_MOD.f NECSEH_MOD.f NECSP_MOD.f
ifx -c -I${MKLROOT}/include/intel64/lp64 -I\"${MKLROOT}/include\" Datagn.f GASYEH.F GASYP.F NecMPCL.f NECMPFLD.f SOMGEH.f SOMGP.f SOMLIB_PEH.F
```

Link:   Options:  Linux/Intel Fortran Compiler/Intel 64/Static/32-bit integers/Threading=Sequential/LAPACK 95
```bash
ifx -o nec5s -static ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a  Datagn.o GASYEH.o GASYP.o NecMPCL.o NECMPFLD.o NECMP_MOD.o NECSEH_MOD.o NECSP_MOD.o SOMGEH.o SOMGP.o SOMLIB_PEH.o -Wl,--end-group -liomp5 -lpthread -lm -ldl
```

## Compiling with the AMD Optimizing Fortran Compiler (AOCC) and AMD Optimizing CPU Libraries (AOCL)

Install the AMD Optimizing C/C++ and Fortran Compilers (AOCC) available [here](https://www.amd.com/en/developer/aocc.html).
Add to .bashrc:  `source /opt/AMD/aocc-compiler-4.1.0/setenv_AOCC.sh`
Install the AMD Optimizing CPU Libraries (AOCL) available [here](https://www.amd.com/en/developer/aocl.html).
Add to .bashrc:  `source /opt/AMD/aocl/aocl-linux-aocc-4.1.0/aocc/amd-libs.cfg`

Some I/O code in NecMPCL.f needs to be tweaked for compatibility with the AMD tools. Here's my effort -- someone with better FORTRAN skills might have other ideas :-)
```bash
cp NecMPCL.f NecMPCL_original.f
patch NecMPCL.f NecMPCL.patch
```

Compile:
```bash
flang -c -march=native NECMP_MOD.f NECSEH_MOD.f NECSP_MOD.f
flang -c -march=native Datagn.f NecMPCL.f NECMPFLD.f SOMGEH.f SOMGP.f GASYEH.F GASYP.F SOMLIB_PEH.F
```

### Link with the AMD shared libraries  -->  nec5amd
Link:
```bash
flang -static-flang-libs -fuse-ld=lld -L/opt/AMD/aocl/aocl-linux-aocc-4.1.0/aocc/lib_LP64/ NECMP_MOD.o NECSEH_MOD.o NECSP_MOD.o Datagn.o NecMPCL.o NECMPFLD.o SOMGEH.o SOMGP.o GASYEH.o GASYP.o SOMLIB_PEH.o -lblis -lflame -lamdlibm -lm -lflang -o nec5amd
```

### Link with the AMD static libraries  -->  nec5amds
Link:
```bash
flang -static-flang-libs -fuse-ld=lld -L/opt/AMD/aocl/aocl-linux-aocc-4.1.0/aocc/lib_LP64/ $LIBROOT/libblis.a $LIBROOT/libflame.a $LIBROOT/libamdlibm.a NECMP_MOD.o NECSEH_MOD.o NECSP_MOD.o Datagn.o NecMPCL.o NECMPFLD.o SOMGEH.o SOMGP.o GASYEH.o GASYP.o SOMLIB_PEH.o -lm -lstdc++ -o nec5amds
```


