#
# project : pansy
# source  : Zen.jl
# author  : Li Huang (lihuang.dmft@gmail.com)
# status  : unstable
# comment :
#
# last modified: 2021/01/04
#

module Zen

"""
    __init__()

This function would be executed immediately after the module is loaded at
runtime for the first time
"""
function __init__() end

#
# using standard library
#
# additional remarks:
#
# the TOML package is included in the standard library since v1.6
# so, please upgrade your julia environment if it is outdated
#
using LinearAlgebra
using Distributed
using Printf
using Dates
using TOML

#
# const.jl
#
# define some global numerical or string constants
#
# summary:
#
# I32, I64    -> numerical types
# F32, F64    -> numerical types
# C32, C64    -> numerical types
# __libname__ -> library's name
# __version__ -> version of zen
# __release__ -> released date of zen
# __authors__ -> authors of zen
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
# and some structs
#
# summary:
#
# PCASE    -> dict for case
# PDFT     -> dict for dft engine
# PDMFT    -> dict for dmft engine
# PIMP     -> dict for quantum impurities
# PSOLVER  -> dict for quantum impurity solver
# IterInfo -> struct for iteration information
# Lattice  -> struct for crystallography information
# PrTrait  -> struct for projectors
# PrGroup  -> struct for groups of projectors
# PrGroupT -> struct for groups of projectors (transformed)
#
include("types.jl")
#
export PCASE
export PDFT
export PDMFT
export PIMP
export PSOLVER
export IterInfo
export Lattice
export PrTrait
export PrGroup
export PrGroupT

#
# config.jl
#
# to extract the configurations from external files or dictionaries
#
# summary:
#
# parse_toml   -> parse case.toml
# renew_config -> update dict (configuration)
# check_config -> check dict (configuration)
# _v           -> verify dict
# _c           -> shortcut for visiting dict (case)
# _d           -> shortcut for visiting dict (dft)
# _m           -> shortcut for visiting dict (dmft)
# _i           -> shortcut for visiting dict (impurity)
# _s           -> shortcut for visiting dict (solver)
#
include("config.jl")
#
export parse_toml
export renew_config
export check_config
export _v
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
# summary:
#
# @cswitch      -> C-style switch
# require       -> check julia envirnoment
# query_args    -> query arguments
# query_inps    -> query input files
# query_zen     -> query home directory of zen
# query_dft     -> query home directory of dft engine
# list_case     -> print dict (case)
# list_dft      -> print dict (dft)
# list_dmft     -> print dict (dmft)
# list_impurity -> print dict (impurity)
# list_solver   -> print dict (solver)
# welcome       -> print welcome message
# overview      -> print overview of zen
# goodbye       -> say goodbye
# sorry         -> say sorry
# message       -> print some message to the screen
# line_to_array -> transform a line to a string array
#
include("util.jl")
#
export @cswitch
export require
export query_args
export query_inps
export query_zen
export query_dft
export list_case
export list_dft
export list_dmft
export list_impurity
export list_solver
export welcome
export overview
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
# summary:
#
# make_trees -> make working directories
# dft_init   -> init dft engine
# dft_run    -> execute dft calculation
# dft_save   -> finalize dft calculation
#
include("base.jl")
#
export make_trees
export adaptor_init
export adaptor_run
export adaptor_save
export dft_init
export dft_run
export dft_save

#
# ir.jl
#
# adaptor for the intermediate representation format
#
# summary:
#
# irio_lattice -> write lattice information
# irio_kmesh   -> write kmesh
# irio_tetra   -> write tetrahedra
# irio_eigen   -> write eigenvalues
# irio_projs   -> write projectors
# irio_fermi   -> write fermi level
# irio_charge  -> write charge density
#
include("ir.jl")
#
export irio_lattice
export irio_kmesh
export irio_tetra
export irio_eigen
export irio_projs
export irio_fermi
export irio_charge

#
# plo.jl
#
# tools for the projection on localized orbitals scheme
#
# summary:
#
# plo_group  -> setup groups of projectors
# plo_rotate -> rotate the projectors
# plo_window -> extract the projectors within a given energy window
# plo_orthog -> orthogonalize the projectors
# plo_diag   -> orthogonalizes a projector defined by a rectangular matrix
# plo_ovlp   -> calculate overlap matrix
# plo_dm     -> calculate density matrix
# plo_hamk   -> calculate local hamiltonian
# plo_dos    -> calculate density of states
# view_ovlp  -> show overlap matrix
# view_dm    -> show density matrix
# view_hamk  -> show local hamiltonian
# view_dos   -> show density of states
#
include("plo.jl")
#
export plo_group
export plo_rotate
export plo_window
export plo_orthog
export plo_diag
export plo_ovlp
export plo_dm
export plo_hamk
export plo_dos
export view_ovlp
export view_dm
export view_hamk
export view_dos

#
# tetra.jl
#
# implementation of analytical tetrahedron method
#
# summary:
#
#
include("tetra.jl")

#
# vasp.jl
#
# adaptor for the vasp software package. it provide a lot of functions
# to deal with the vasp-related files
#
# summary:
#
# vasp_init      -> prepare vasp's input files
# vasp_run       -> execute vasp program
# vasp_save      -> backup vasp's output files
# vasp_incar     -> make essential input file (INCAR)
# vasp_kpoints   -> make essential input file (KPOINTS)
# vasp_files     -> check essential output files
# vaspio_lattice -> read lattice information
# vaspio_kmesh   -> read kmesh
# vaspio_tetra   -> read tetrahedra
# vaspio_eigen   -> read eigenvalues
# vaspio_projs   -> read projectors
# vaspio_fermi   -> read fermi level
# vaspio_charge  -> read charge density
#
include("vasp.jl")
#
export vasp_init
export vasp_run
export vasp_save
export vasp_incar
export vasp_kpoints
export vasp_files
export vaspio_lattice
export vaspio_kmesh
export vaspio_tetra
export vaspio_eigen
export vaspio_projs
export vaspio_fermi
export vaspio_charge

end
