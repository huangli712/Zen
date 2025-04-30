##
##
## Introduction
## ============
##
## This is the top Makefile for the Zen software package. You can use
## it to build all of the executable programs and libraries. Please do
## not modify it by yourself even you are very familiar with Zen.
##
## Usage
## =====
##
## 1. get help
## -----------
##
## 'make help'
##
## 2. quick build
## --------------
##
## 'make all'
##
## 3. build internal library
## -------------------------
##
## 'make libs'
##
## 4. build impurity solver
## ------------------------
##
## 'make solver'
##
## 5. build dmft engine
## --------------------
##
## 'make dmft'
##
## 6. build documentation
## ----------------------
##
## 'make docs'
##
## 7. clean the directories
## ------------------------
##
## 'make clean'
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
## 08/06/2021 by li huang (created)
## 04/30/2025 by li huang (last modified)
##
##

# targets for help
help:
	@echo 'Making The Zen DFT + DMFT Computation Framework'
	@echo
	@echo
	@echo 'targets:'
	@echo '  help          print full help information            '
	@echo
	@echo '  all           build all components and tools         '
	@echo '  clean         clean all directories                  '
	@echo
	@echo '  libs          build internal libraries               '
	@echo '  clean-libs    clean internal libraries               '
	@echo
	@echo '  solver        build quantum impurity solvers         '
	@echo '  clean-solver  clean quantum impurity solvers         '
	@echo
	@echo '  dmft          build dynamical mean-field theory codes'
	@echo '  clean-dmft    clean dynamical mean-field theory codes'
	@echo
	@echo '  docs          build documentation                    '
	@echo '  clean-docs    clean documentation                    '

# all-in-one building target
all: libs solver dmft docs

# all-in-one cleaning target
clean: clean-libs clean-solver clean-dmft clean-docs

# targets for internal libraries
libs:
	@echo 'building libs: flink'
	cd ../src/libs/Flink/src; pwd; make

clean-libs:
	@echo 'cleaning libs: flink'
	cd ../src/libs/Flink/src; pwd; make clean

# targets for quantum impurity solvers
solver:
	@echo 'building solver: ctseg'
	cd ../src/solver/iQIST/src/ctseg; pwd; make

clean-solver:
	@echo 'cleaning solver: ctseg'
	cd ../src/solver/iQIST/src/ctseg; pwd; make clean

# targets for dynamical mean-field theory codes
dmft:
	@echo 'building dmft: dyson'
	cd ../src/dmft/Dyson/src; pwd; make

clean-dmft:
	@echo 'cleaning dmft: dyson'
	cd ../src/dmft/Dyson/src; pwd; make clean

# targets for documentation
docs:
	@echo 'building docs'
	cd ../docs/; pwd; julia make.jl

clean-docs:
	@echo 'cleaning docs'
	cd ../docs/; pwd; rm -fr build
