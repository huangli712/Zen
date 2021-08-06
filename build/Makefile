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
##
##

# targets for help
help:
	@echo '  Interacting Quantum Impurity Solver Toolkit Building System'
	@echo
	@echo
	@echo '  targets:'
	@echo '     help              print full usage information          '
	@echo
	@echo '     all               build all components and tools        '
	@echo '     clean             clean all directories                 '
	@echo
	@echo '     solver            build all quantum impurity solvers    '
	@echo '     clean-solver      clean all quantum impurity solvers    '
	@echo
	@echo '     atomic            build atomic eigenvalue solver        '
	@echo '     clean-atomic      clean atomic eigenvalue solver        '
	@echo
	@echo '     ct-hyb1           build ct-hyb1 quantum impurity solver '
	@echo '     clean-ct-hyb1     clean ct-hyb1 quantum impurity solver '
	@echo '     ct-hyb2           build ct-hyb2 quantum impurity solver '
	@echo '     clean-ct-hyb2     clean ct-hyb2 quantum impurity solver '

# all-in-one building target
all: solver atomic

# all-in-one cleaning target
clean: clean-solver clean-atomic

# targets for solver
solver: ct-hyb1 ct-hyb2
clean-solver: clean-ct-hyb1 clean-ct-hyb2

# targets for atomic eigenvalue solver: jasmine
atomic:
	cd ../src/atomic/; pwd; make

clean-atomic:
	cd ../src/atomic/; pwd; make clean

# targets for ct-hyb quantum impurity solver: segment version (narcissus)
ct-hyb1:
	cd ../src/ct_hyb1/; pwd; make

clean-ct-hyb1:
	cd ../src/ct_hyb1/; pwd; make clean

# targets for ct-hyb quantum impurity solver: general version (manjushaka)
ct-hyb2:
	cd ../src/ct_hyb2/; pwd; make

clean-ct-hyb2:
	cd ../src/ct_hyb2/; pwd; make clean
