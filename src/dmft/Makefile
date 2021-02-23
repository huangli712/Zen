##
##
## Introduction
## ============
##
## Makefile for zen/src/dmft directory.
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
## Li Huang // email: lihuang.dmft@gmail.com
##
## History
## =======
##
## 02/23/2021 by li huang (created)
## 02/23/2021 by li huang (last modified)
##
##

.SUFFIXES: .f90

include ../../build/make.sys

modc = dmft_control.o dmft_context.o
core = dmft_core.o dmft_task.o
lev1 = dmft_stream.o dmft_util.o
lev2 = dmft_dump.o dmft_print.o
main = dmft_main.o
mlib = libflink.a

objects = $(modc) $(core) $(lev1) $(lev2) $(main) $(mlib)

default: all

all: exe

flink: flink_lib flink_mod

flink_lib:
	cp $(FLINK)/libflink.a .

flink_mod:
	cp $(FLINK)/constants.mod .
	cp $(FLINK)/mmpi.mod .
	cp $(FLINK)/spring.mod .
	cp $(FLINK)/stack.mod .
	cp $(FLINK)/parser.mod .
	cp $(FLINK)/linkedlist.mod .
	cp $(FLINK)/version.mod .

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
