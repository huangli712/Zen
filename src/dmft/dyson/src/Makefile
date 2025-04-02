##
##
## Introduction
## ============
##
## Makefile for dyson/src directory.
##
## Usage
## =====
##
## 'make', 'make clean', 'make clean-dat' or 'make clean-all'
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
## 02/23/2021 by li huang (created)
## 03/26/2025 by li huang (last modified)
##
##

.SUFFIXES: .f90

include ../../../../build/make.inc

modc = dmft_control.o dmft_context.o
core = dmft_task.o dmft_flow.o dmft_core.o
lev1 = dmft_psichi.o
lev2 = dmft_stream.o dmft_dump.o dmft_print.o
lev3 = dmft_util.o
main = dmft_main.o
mlib = libflink.a

objects = $(modc) $(core) $(lev1) $(lev2) $(lev3) $(main) $(mlib)

default: all

all: exe

flink: flink_lib flink_mod

flink_lib:
	cp $(FLINK)/libflink.a .

flink_mod:
	cp $(FLINK)/constants.mod .
	cp $(FLINK)/mmpi.mod .
	cp $(FLINK)/parser.mod .
	cp $(FLINK)/linkedlist.mod .

exe: flink $(objects)
	$(LINKER) $(objects) -o dmft $(LFLAGS) $(LIBS)

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f dmft
	rm -f libflink.a

clean-dat:
	rm -f *.dat
	rm -f *.out

clean-all: clean clean-dat
