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
## 3. build deps library
## ---------------------
##
## 'make deps'
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
## Li Huang // email: lihuang.dmft@gmail.com
##
## History
## =======
##
## 08/06/2021 by li huang (created)
## 08/06/2021 by li huang (last modified)
## 01/05/2023 by jia-ming wang (last modified)
##
##

# targets for help
help:
	@echo 'Making The Zen DFT + DMFT Computation Framework'
	@echo
	@echo
	@echo 'targets:'
	@echo '  help          print full help information       '
	@echo
	@echo '  all           build all components and tools    '
	@echo '  clean         clean all directories             '
	@echo
	@echo '  deps          build essential library           '
	@echo '  clean-deps    clean essential library           '
	@echo
	@echo '  solver        build quantum impurity solvers    '
	@echo '  clean-solver  clean quantum impurity solvers    '
	@echo
	@echo '  dmft          build dynamical mean-field engine '
	@echo '  clean-dmft    clean dynamical mean-field engine '
	@echo
	@echo '  docs          build documentations              '
	@echo '  clean-docs    clean documentations              '

# all-in-one building target
all: deps solver dmft docs

# all-in-one cleaning target
clean: clean-deps clean-solver clean-dmft clean-docs

# targets for internal library
deps:
	@echo 'building deps'
	cd ../src/deps/; pwd; make

clean-deps:
	@echo 'cleaning deps'
	cd ../src/deps/; pwd; make clean

# targets for quantum impurity solvers
solver:
	@echo 'building solvers / ct_hyb1'
	cd ../src/solver/ct_hyb1/; pwd; make
	@echo 'building solvers / norg'
	cd ../src/solver/norg/; pwd; make

clean-solver:
	@echo 'cleaning solvers / ct_hyb1'
	cd ../src/solver/ct_hyb1/; pwd; make clean
	@echo 'cleaning solvers / norg'
	cd ../src/solver/norg/; pwd; make clear

# targets for dynamical mean-field engine
dmft:
	@echo 'building dmft'
	cd ../src/dmft/; pwd; make

clean-dmft:
	@echo 'cleaning dmft'
	cd ../src/dmft/; pwd; make clean

# targets for documentations
docs:
	@echo 'building docs'
	cd ../docs/; pwd; julia make.jl

clean-docs:
	@echo 'cleaning docs'
	cd ../docs/; pwd; rm -fr build
