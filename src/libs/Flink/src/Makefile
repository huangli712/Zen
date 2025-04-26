##
##
## Introduction
## ============
##
## Makefile for flink/src directory.
##
## Usage
## =====
##
## 'make' or 'make clean'
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
## 09/18/2014 by li huang (created)
## 12/30/2024 by li huang (last modified)
##
##

.SUFFIXES: .f90

include ../../../../build/make.inc

mods1 = m_constants.o
mods2 = m_linkedlist.o m_stack.o
mods3 = m_face.o m_mpi.o m_parser.o m_sparse.o m_spring.o m_tetra.o
subs1 = s_error.o s_util.o
subs2 = s_fourier.o s_function.o s_integrator.o s_spline.o
subs3 = s_matrix.o s_vector.o
objects = $(mods1) $(mods2) $(mods3) $(subs1) $(subs2) $(subs3)

default: all

all: lib

lib: $(objects)
	$(ARCHIVER) libflink.a $(objects)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f libflink.a

clean-dat:
	rm -f *.dat
	rm -f *.out

clean-all: clean clean-dat
