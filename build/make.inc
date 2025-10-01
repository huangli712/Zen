##
##
## Introduction
## ============
##
## This is the official configuration file for the building system. You
## should modify it to fulfill your requirements. The make.inc file is
## the key component of the building system. If it is not configured
## correctly, the building system won't work correctly as well.
##
## Author
## ======
##
## This building system is designed, created, and maintained by
##
## Li Huang // email: huangli@caep.cn
##
## History
## =======
##
## 05/11/2015 by li huang (created)
## 04/30/2025 by li huang (last modified)
##
##

# Fortran compiler, linker, and archiver
#-------------------------------------------------------------------------
F90    = mpiifx
LINKER = $(F90)
ARCHIVER = ar -rv

# Fortran preprocessor options
#-------------------------------------------------------------------------
MPI    = -DMPI
OMP    = #-qopenmp
FPP    = -fpp
CPP    = $(FPP) $(MPI) $(OMP)

# Machine tuning options
#-------------------------------------------------------------------------
CHECK  = -nogen-interfaces -warn all,noexternal #-check all -traceback -g
MTUNE  = -O3 -xHost

# Flags for fortran compiler and linker
#-------------------------------------------------------------------------
FFLAGS = -c $(CPP) $(CHECK) $(MTUNE)
LFLAGS = $(OMP)

# External linear algebra library
#-------------------------------------------------------------------------
LIBS   = -L$(MKLROOT)/lib -qmkl
