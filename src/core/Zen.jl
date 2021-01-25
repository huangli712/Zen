#
# Project : Pansy
# Source  : Zen.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
# Comment :
#
# Last modified: 2021/01/25
#

module Zen

#
# Using standard libraries
#
using LinearAlgebra
using Distributed
using Printf
using Dates

#
# Using third-party libraries
#
# Remarks:
#
# The TOML package is included in the standard library since v1.6.
# So, please upgrade your julia environment if it is outdated.
#
using TOML

#
# global.jl
#
# Summary:
#
# Define type aliases and some string constants for the Zen framework.
#
# Members:
#
# I32, I64    -> Numerical types (integer)
# F32, F64    -> Numerical types (float)
# C32, C64    -> Numerical types (complex)
# __LIBNAME__ -> Name of this framework
# __VERSION__ -> Version of this framework
# __RELEASE__ -> Released date of this framework
# __AUTHORS__ -> Authors of this framework
#
include("global.jl")
#
export I32, I64
export F32, F64
export C32, C64
export __LIBNAME__
export __VERSION__
export __RELEASE__
export __AUTHORS__

#
# types.jl
#
# Summary:
#
# Define some dicts and structs, which store the config parameters or
# represent some essential data structures.
#
# Members:
#
# PCASE    -> Dict for case
# PDFT     -> Dict for DFT engine
# PDMFT    -> Dict for DMFT engine
# PIMP     -> Dict for quantum impurity problems
# PSOLVER  -> Dict for quantum impurity solvers
# DFTData  -> Dict for storing Kohn-Sham band structure and related data
# Logger   -> Struct for logger
# IterInfo -> Struct for DFT + DMFT iteration information
# Lattice  -> Struct for crystallography information
# PrTrait  -> Struct for projectors
# PrGroup  -> Struct for groups of projectors
# PrGroupT -> Struct for groups of projectors (transformed)
#
include("types.jl")
#
export PCASE
export PDFT
export PDMFT
export PIMP
export PSOLVER
export DFTData
export Logger
export IterInfo
export Lattice
export PrTrait
export PrGroup
export PrGroupT

#
# config.jl
#
# Summary:
#
# To extract the configurations from external files or dictionaries.
#
# Members:
#
# setup    -> Setup parameters
# exhibit  -> Display parameters for reference
# inp_toml -> Parse case.toml, return raw configuration information
# new_dict -> Update dicts for configuration
# chk_dict -> Check dicts for configuration
# _v       -> Verify dict's values
# cat_c    -> Print dict (PCASE dict)
# cat_d    -> Print dict (PDFT dict)
# cat_m    -> Print dict (PDMFT dict)
# cat_i    -> Print dict (PIMP dict)
# cat_s    -> Print dict (PSOLVER dict)
# get_c    -> Extract value from dict (PCASE dict), return original value
# get_d    -> Extract value from dict (PDFT dict), return original value
# get_m    -> Extract value from dict (PDMFT dict), return original value
# get_i    -> Extract value from dict (PIMP dict), return original value
# get_s    -> Extract value from dict (PSOLVER dict), return original value
# str_c    -> Extract value from dict (PCASE dict), return everything as string
# str_d    -> Extract value from dict (PDFT dict), return everything as string
# str_m    -> Extract value from dict (PDMFT dict), return everything as string
# str_i    -> Extract value from dict (PIMP dict), return everything as string
# str_s    -> Extract value from dict (PSOLVER dict), return everything as string
#
include("config.jl")
#
export setup
export exhibit
export inp_toml
export new_dict
export chk_dict
export _v
export cat_c
export cat_d
export cat_m
export cat_i
export cat_s
export get_c
export get_d
export get_m
export get_i
export get_s
export str_c
export str_d
export str_m
export str_i
export str_s

#
# util.jl
#
# Summary:
#
# To provide some useful utility functions. They can be used to query
# the environments and parse the strings, etc.
#
# Members:
#
# @cswitch      -> C-style switch
# @ps1          -> Wrapper for printstyled function
# @ps2          -> Wrapper for printstyled function
# require       -> Check julia envirnoment
# query_args    -> Query program's arguments
# query_case    -> Query case (job's name)
# query_inps    -> Query input files
# query_stop    -> Query case.stop file
# query_zen     -> Query home directory of Zen
# query_dft     -> Query home directory of DFT engine
# welcome       -> Print welcome message
# overview      -> Print overview of Zen
# goodbye       -> Say goodbye
# sorry         -> Say sorry
# prompt        -> Print some messages to the screen
# line_to_array -> Convert a line to a string array
#
include("util.jl")
#
export @cswitch
export @ps1
export @ps2
export require
export query_args
export query_case
export query_inps
export query_stop
export query_zen
export query_dft
export welcome
export overview
export goodbye
export sorry
export prompt
export line_to_array

#
# base.jl
#
# Summary:
#
# To provide the core functions to control the DFT engine, DMFT engine,
# and impurity solvers.
#
# Members:
#
# ready        -> Prepare everything that is essential for DFT + DMFT calculations
# go           -> Dispatcher for DFT + DMFT calculations
# final        -> Finalize the DFT + DMFT calculations
# cycle1       -> Perform DFT + DMFT calculations (one-shot mode)
# cycle2       -> Perform DFT + DMFT calculations (fully self-consistent mode)
# monitor      -> Monitor the DFT + DMFT calculations
# make_trees   -> Make working directories
# rm_trees     -> Remove working directories
# adaptor_init -> Initialize DFT_DMFT adaptor
# adaptor_run  -> Launch DFT_DMFT adaptor
# adaptor_save -> Backup files generated by DFT_DMFT adaptor
# dft_init     -> Initialize DFT engine
# dft_run      -> Launch DFT engine
# dft_save     -> Backup files generated by DFT engine
# dmft_init    -> Initialize DMFT engine
# dmft_run     -> Launch DMFT engine
# dmft_save    -> Backup files generated by DMFT engine
# solver_init  -> Initialize quantum impurity solvers
# solver_run   -> Launch quantum impurity solvers
# solver_save  -> Backup files generated by quantum impurity solvers
#
include("base.jl")
#
export ready
export go
export final
export cycle1
export cycle2
export monitor
export make_trees
export rm_trees
export adaptor_init
export adaptor_run
export adaptor_save
export dft_init
export dft_run
export dft_save
export dmft_init
export dmft_run
export dmft_save
export solver_init
export solver_run
export solver_save

#
# ir.jl
#
# Summary:
#
# Tools for the intermediate representation format (adaptor).
#
# Members:
#
# ir_adaptor   -> Adaptor support
# irio_lattice -> Write lattice information
# irio_kmesh   -> Write kmesh
# irio_tetra   -> Write tetrahedra
# irio_eigen   -> Write eigenvalues
# irio_projs   -> Write projectors
# irio_fermi   -> Write fermi level
# irio_charge  -> Write charge density
#
include("ir.jl")
#
export ir_adaptor
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
# Summary:
#
# Tools for the projection on localized orbitals scheme (adaptor).
#
# Members:
#
# plo_adaptor -> Adaptor support
# plo_group   -> Setup groups of projectors
# plo_rotate  -> Rotate the projectors
# plo_filter  -> Extract the projectors within a given energy window
# plo_orthog  -> Orthogonalize the projectors
# plo_diag    -> Orthogonalizes a projector defined by a rectangular matrix
# plo_ovlp    -> Calculate overlap matrix
# plo_dm      -> Calculate density matrix
# plo_hamk    -> Calculate local hamiltonian
# plo_dos     -> Calculate density of states
# view_ovlp   -> Show overlap matrix
# view_dm     -> Show density matrix
# view_hamk   -> Show local hamiltonian
# view_dos    -> Show density of states
#
include("plo.jl")
#
export plo_adaptor
export plo_group
export plo_rotate
export plo_filter
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
# Summary:
#
# Tools for the analytical tetrahedron method (adaptor).
#
# Members:
#
#
include("tetra.jl")

#
# vasp.jl
#
# Summary:
#
# Tools for the vasp software package (adaptor). It provide a lot of
# functions to deal with the vasp-related files.
#
# Members:
#
# vasp_adaptor   -> Adaptor support
# vasp_init      -> Prepare vasp's input files
# vasp_run       -> Execute vasp program
# vasp_save      -> Backup vasp's output files
# vasp_incar     -> Make essential input file (INCAR)
# vasp_kpoints   -> Make essential input file (KPOINTS)
# vasp_files     -> Check essential output files
# vaspio_lattice -> Read lattice information
# vaspio_kmesh   -> Read kmesh
# vaspio_tetra   -> Read tetrahedra
# vaspio_eigen   -> Read eigenvalues
# vaspio_projs   -> Read projectors
# vaspio_fermi   -> Read fermi level
# vaspio_charge  -> Read charge density
#
include("vasp.jl")
#
export vasp_adaptor
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

"""
    __init__()

This function would be executed immediately after the module is loaded at
runtime for the first time.
"""
function __init__() end

end
