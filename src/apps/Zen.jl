#
# project :
# program :
# source  :
# type    :
# author  :
# history :
# purpose :
# status  :
# comment :
#

module Zen

#
# using standard library
#
# additional remarks:
#
# the TOML package is included in the standard library since v1.6
# so, please upgrade your julia environment if it is outdated 
#
using Printf
using TOML

#
# const.jl
#
# define some global numerical or string constants
#
include("const.jl")
#
export I32, I64
export F32, F64
export C32, C64
export __libname__
export __version__
export __release__
export __authors__

#
# types.jl
#
# define some dictionaries which contain the configuration parameters
# and some data structures
#
include("types.jl")
#
export PCASE
export PDFT
export PDMFT
export PIMP
export PSOLVER
export IterInfo

#
# config.jl
#
# to extract the configurations from external files or dictionaries 
#
include("config.jl")
#
export parse_toml
export renew_config
export check_config
export _c
export _d
export _m
export _i
export _s

#
# util.jl
#
# to provide some useful utility functions. they can be used tp query
# the environments, print the configurations, and parse the strings, etc. 
#
include("util.jl")
#
export require
export query_args
export query_cars
export query_zen
export query_dft
export view_case
export view_dft
export view_dmft
export view_impurity
export view_solver
export welcome
export goodbye
export sorry
export message
export line_to_array

#
# base.jl
#
# to provide the core functions to control the dft engine, dmft engine,
# and impurity solver
#
include("base.jl")
#
export make_trees
export make_incar
export dft_init
export dft_run
export dft_save

#
# vasp.jl
#
# adaptor for the vasp software package. it provide a lot of functions
# to deal with the vasp-related files
#
include("vasp.jl")
#
export vaspio_poscar
export vaspio_ibzkpt
export vaspio_projcar
export vaspio_locproj
export vaspio_eigenval
export vaspio_chgcar

#
# ir.jl
#
# adaptor for the intermediate representation format
#
include("ir.jl")
#
export irio_lattice
export irio_kmesh
export irio_tetra
export irio_projs
export irio_eigen

end
