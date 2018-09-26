MPILIB = -lmpi
prefix = /home/weast/

LIBS       = -lpamr -lamrd $(MPILIB) -lm -lgfortran -lbbhutil  
LDFLAGS    =  -L../../lib  -L. -L$(prefix)/lib -L/home/weast/software/pamr/lib/ 

CC       = mpicc 
CFLAGS     = -g -O2 
CPPFLAGS = -I.  -I../../include -I$(prefix)/include -I/home/weast/software/pamr/include/ 
CC_COMP  = $(CC) -c $(CFLAGS)

F_COMP = mpif77 -c

SHEN_OBJS = loadTable.o shenInterface.o

OBJS = blackHoleFinder.o fluxCorrection.o gravityElliptic.o initialData.o fluidEvolve.o initialVelocity.o rungeKutta.o TOV.o gravity.o boundaryConditions.o flux.o matrixInverse.o conservedPrimitive.o equationOfState.o brent.o ss_gr_hydro.o $(SHEN_OBJS)

all: ss_gr_hydro



ss_gr_hydro: $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJS) $(LIBS) -o ss_gr_hydro



flux.o: flux.c variableDef.h
	$(CC_COMP) flux.c

matrixInverse.o: matrixInverse.c matrixInverse.h
	$(CC_COMP) matrixInverse.c

conservedPrimitive.o: conservedPrimitive.c variableDef.h conservedPrimitive.h
	$(CC_COMP) conservedPrimitive.c

equationOfState.o: equationOfState.c equationOfState.h
	$(CC_COMP) equationOfState.c

#For Shen EOS
loadTable.o: loadTable.c
	$(CC_COMP) loadTable.c
shenInterface.o: shenInterface.c
	$(CC_COMP) shenInterface.c
#End Shen EOS

brent.o: brent.c brent.h
	$(CC_COMP) brent.c

boundaryConditions.o: boundaryConditions.c boundaryConditions.h
	$(CC_COMP) boundaryConditions.c

gravity.o: gravity.c gravity.h
	$(CC_COMP) gravity.c

gravityElliptic.o: gravityElliptic.c gravityElliptic.h
	$(CC_COMP) gravityElliptic.c

TOV.o: TOV.c TOV.h
	$(CC_COMP) TOV.c

rungeKutta.o: rungeKutta.c rungeKutta.h
	$(CC_COMP) rungeKutta.c

initialVelocity.o: initialVelocity.c initialVelocity.h
	$(CC_COMP) initialVelocity.c

initialData.o: initialData.c initialData.h
	$(CC_COMP) initialData.c

blackHoleFinder.o: blackHoleFinder.c blackHoleFinder.h
	$(CC_COMP) blackHoleFinder.c

fluidEvolve.o: fluidEvolve.c fluidEvolve.h
	$(CC_COMP) fluidEvolve.c

fluxCorrection.o: fluxCorrection.c fluxCorrection.h
	$(CC_COMP) fluxCorrection.c

ss_gr_hydro.o: ss_gr_hydro.c ss_gr_hydro.h
	$(CC) $(CPPFLAGS) -c $(CFLAGS) ss_gr_hydro.c