# How to run this code in different architecture

## Build PAMR/AMRD

### On Marylou

### On Comet
`./configure --prefix=/home/<your username>/local --CC=mpicc`

We recommend to use mpicc as cc because somehow ti complains
about cannot find `mpi.h` even you call mpi module

Once you didn't see anyproblme during configuration, then type
```
make
```
### On Mac

## Build the code

### On Marylou

Use Makefile.m7 then change to your user name

### On Comet

Use Makefile.comet then change to your user name

### On Mac
