# How to run this code in different architectures

## Build PAMR/AMRD

To obtain PAMR/AMRD library, you need to contact William East 
(weast@perimeterinstitute.ca). Once you clone the repository,
follow the instruction that is suitable for your systme

### On Marylou

Working modules to build the library for Marylou are below
```
  1) gcc/6.4                 4) defaultenv/4(default)   7) hdf5/1.8.16_gnu-5.2.0
  2) gdb/8.1                 5) python/3/4
  3) openmpi/3.0_gcc-6       6) gsl/1.16
```

Once you have these modules, you need to configure it first via

`./configure --prefix=/home/<your username>/local`

Once you didn't see any problmes during configuration, then type
```
make install
```
If you build the library correctly, you should fine `libpamr.a`
and `libamrd.a` in your `/home/<your username>/local/lib` 
directory


### On Comet

Working modules to build the library for Comet are below
```
  1) gnutools/2.69         3) intelmpi/2016.3.210   5) mvapich2_ib/2.1
  2) intel/2016.3.210      4) gsl/2.1               6) hdf5/1.8.14
```
Once you have these modules, you need to configure it first via

`./configure --prefix=/home/<your username>/local CC=mpicc`

We recommend to use mpicc as your c complier because somehow it complains
about cannot find `mpi.h` even you call mpi module

Once you didn't see any problmes during configuration, then type
```
make install
```
If you build the library correctly, you should fine `libpamr.a`
and `libamrd.a` in your `/home/<your username>/local/lib` 
directory

### On Symmetry

Working modules to build the library for Symmetry are below
```
 1) julia/1.0.3(default)   3) hdf5_18/1.8.20        5) gcc/7.2.0
 2) slurm/17.11.8          4) python/3.7(default)   6) openmpi/gcc/64/1.10.7
```
Once you have these modules, you need to configure it first via

`./configure --prefix=/home/<your username>/local`

Once you didn't see any problmes during configuration, then type
```
make install
```
If you build the library correctly, you should fine `libpamr.a`
and `libamrd.a` in your `/home/<your username>/local/lib` 
directory

### On Mac

We highly recommend to use `GNU` complier instead of `clang`. You can easily install `GNU` 
complier via various ways such as `homebrew`, `fink`, or manually. Once you successfully install
gnu complier then install `OpenMPI` with `GNU` complier. During configuration, you must specify
`GNU` complier. 

After all installations, you can simply follow similar way as above. For example

`./configure CC=<your GNU gcc> CXX=<your GNU g++> --prefix=/usr/local`

Once you didn't see any problmes during configuration, then type
```
make install
```
If you build the library correctly, you should fine `libpamr.a`
and `libamrd.a` in your `/usr/local/lib` 
directory

Note that there is a problem with `malloc.h` in `cls_.c`. Since Xcode compliation 
can be different with usual compliers so it stores basic libs in different spots.
Especially, if you are using package management like `homebrew` and `macport`, you may
have multiple directories that save bin and libs. Most of case for Mac, `malloc.h` is 
located  in `/usr/include/malloc` rather than `/usr/local/include`. To fix this problem, 
you can simply add this path during configuration 

`./configure CC=<your GNU gcc> CXX=<your GNU g++> CPPFLAGS=-I/usr/include/malloc --prefix=/usr/local`

Then type `make install` to build.


## Build the code

Once you build the PAMR/AMRD correct, you can simply build the code by
typing
```
make
```
If you would like to use different evolution formalism for the code,
you can type
```
make BSSN=1
```
to enable addtional features

In this directory, different `Makefile` is saved for users.

### On Marylou

Use Makefile.m7 then change to your user name and/or path

### On Comet

Use Makefile.comet then change to your user name and/or path

### On Symmetry

Use Makefile.symmetry then change to your user name and/or path

### On Mac

Use Makefile.Mac. If you follow the instruction above, you do not
need to change anything

## Contact
Other architectures should work if you have similar/same modules as
above. If you have questions, please contact Hyun Lim (hylim1988@gmail.com)
