#
# Project : Pansy
# Source  : Zen.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
# Comment :
#
# Last modified: 2021/02/14
#

"""
    Zen

Zen is a modern DFT + DMFT computation framework.
"""
module Zen

#
# Using standard libraries
#

using LinearAlgebra
using Distributed
using Printf
using Dates
using Base.Math: libm

#
# Using third-party libraries
#

#
# Remarks:
#
# The TOML.jl package is included in the standard library since v1.6.
# So, please upgrade your julia environment if it is outdated.
#
# We need this package to parse the configuration file.
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
# I32, I64    -> Numerical types (Integer).
# F32, F64    -> Numerical types (Float).
# C32, C64    -> Numerical types (Complex).
# R32, R64    -> Numerical types (Union of Integer and Float).
# N32, N64    -> Numerical types (Union of Integer, Float, and Complex).
# __LIBNAME__ -> Name of this package.
# __VERSION__ -> Version of this package.
# __RELEASE__ -> Released date of this package.
# __AUTHORS__ -> Authors of this package.
#
include("global.jl")
#
export I32, I64
export F32, F64
export C32, C64
export R32, R64
export N32, N64
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
# PCASE    -> Dict for case.
# PDFT     -> Dict for DFT engine.
# PDMFT    -> Dict for DMFT engine.
# PIMP     -> Dict for quantum impurity problems.
# PSOLVER  -> Dict for quantum impurity solvers.
# Logger   -> Struct for logger.
# IterInfo -> Struct for DFT + DMFT iteration information.
# Lattice  -> Struct for crystallography information.
# PrTrait  -> Struct for projectors.
# PrGroup  -> Struct for groups of projectors.
# PrWindow -> Struct for band window.
#
include("types.jl")
#
export PCASE
export PDFT
export PDMFT
export PIMP
export PSOLVER
export Logger
export IterInfo
export Lattice
export PrTrait
export PrGroup
export PrWindow

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
# @cswitch      -> C-style switch.
# @ps1          -> Wrapper for printstyled function.
# @ps2          -> Wrapper for printstyled function.
# require       -> Check julia envirnoment.
# setup_args    -> Setup ARGS manually.
# query_args    -> Query program's arguments.
# query_case    -> Query case (job's name).
# query_inps    -> Query input files.
# query_stop    -> Query case.stop file.
# query_zen     -> Query home directory of Zen.
# query_dft     -> Query home directory of DFT engine.
# welcome       -> Print welcome message.
# overview      -> Print overview of Zen.
# goodbye       -> Say goodbye.
# sorry         -> Say sorry.
# prompt        -> Print some messages to the device.
# line_to_array -> Convert a line to a string array.
# line_to_cmplx -> Convert a line to a cmplx number.
# erf           -> Gauss error function.
#
include("util.jl")
#
export @cswitch
export @ps1
export @ps2
export require
export setup_args
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
export line_to_cmplx
export erf

#
# config.jl
#
# Summary:
#
# To extract the configurations from external files or dictionaries.
#
# Members:
#
# setup    -> Setup parameters.
# exhibit  -> Display parameters for reference.
# inp_toml -> Parse case.toml, return raw configuration information.
# rev_dict -> Update dicts for configuration.
# chk_dict -> Check dicts for configuration.
# _v       -> Verify dict's values.
# cat_c    -> Print dict (PCASE dict).
# cat_d    -> Print dict (PDFT dict).
# cat_m    -> Print dict (PDMFT dict).
# cat_i    -> Print dict (PIMP dict).
# cat_s    -> Print dict (PSOLVER dict).
# get_c    -> Extract value from dict (PCASE dict), return raw value.
# get_d    -> Extract value from dict (PDFT dict), return raw value.
# get_m    -> Extract value from dict (PDMFT dict), return raw value.
# get_i    -> Extract value from dict (PIMP dict), return raw value.
# get_s    -> Extract value from dict (PSOLVER dict), return raw value.
# str_c    -> Extract value from dict (PCASE dict), return string.
# str_d    -> Extract value from dict (PDFT dict), return string.
# str_m    -> Extract value from dict (PDMFT dict), return string.
# str_i    -> Extract value from dict (PIMP dict), return string.
# str_s    -> Extract value from dict (PSOLVER dict), return string.
#
include("config.jl")
#
export setup
export exhibit
export inp_toml
export rev_dict
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
# base.jl
#
# Summary:
#
# To provide the core functions to control the DFT engine, DMFT engine,
# and quantum impurity solvers.
#
# Members:
#
# ready        -> Prepare runtime environment for DFT + DMFT calculations.
# go           -> Dispatcher for DFT + DMFT calculations.
# final        -> Finalize the DFT + DMFT calculations.
# cycle1       -> Perform DFT + DMFT calculations (one-shot mode).
# cycle2       -> Perform DFT + DMFT calculations (fully self-consistent mode).
# monitor      -> Monitor the DFT + DMFT calculations.
# make_trees   -> Make working directories.
# rm_trees     -> Remove working directories.
# adaptor_init -> Initialize DFT_DMFT adaptor.
# adaptor_run  -> Launch DFT_DMFT adaptor.
# adaptor_save -> Backup files generated by DFT_DMFT adaptor.
# dft_init     -> Initialize DFT engine.
# dft_run      -> Launch DFT engine.
# dft_save     -> Backup files generated by DFT engine.
# dmft_init    -> Initialize DMFT engine.
# dmft_run     -> Launch DMFT engine.
# dmft_save    -> Backup files generated by DMFT engine.
# solver_init  -> Initialize quantum impurity solvers.
# solver_run   -> Launch quantum impurity solvers.
# solver_save  -> Backup files generated by quantum impurity solvers.
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
# ir_adaptor   -> Adaptor support.
# ir_save      -> Save the output files by the adaptor.
# irio_params  -> Write key parameters extracted from Kohn-Sham data.
# irio_lattice -> Write lattice information.
# irio_kmesh   -> Write kmesh.
# irio_tetra   -> Write tetrahedra.
# irio_eigen   -> Write eigenvalues.
# irio_projs   -> Write projectors.
# irio_fermi   -> Write fermi level.
# irio_charge  -> Write charge density.
#
include("ir.jl")
#
export ir_adaptor
export ir_save
export irio_params
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
# plo_adaptor -> Adaptor support.
# plo_fermi   -> Calibrate fermi level for Kohn-Sham eigenvalues.
# plo_group   -> Setup groups of projectors.
# plo_window  -> Setup band window of projectors.
# plo_rotate  -> Rotate the projectors.
# plo_filter  -> Extract the projectors within a given energy window.
# plo_orthog  -> Orthogonalize / normalize the projectors.
# plo_monitor -> Generate some physical quantities using the projectors.
# get_win1    -> Evaluate band window.
# get_win2    -> Evaluate energy window.
# try_blk1    -> Orthogonalize / normalize the projectors group by group.
# try_blk2    -> Orthogonalize / normalize the projectors with each other.
# try_diag    -> Orthogonalizes a projector defined by a rectangular matrix.
# calc_ovlp   -> Calculate overlap matrix.
# calc_dm     -> Calculate density matrix.
# calc_hamk   -> Calculate local hamiltonian.
# calc_dos    -> Calculate density of states.
# view_ovlp   -> Show overlap matrix for debug.
# view_dm     -> Show density matrix for debug.
# view_hamk   -> Show local hamiltonian for debug.
# view_dos    -> Show density of states for debug.
#
include("plo.jl")
#
export plo_adaptor
export plo_fermi
export plo_group
export plo_window
export plo_rotate
export plo_filter
export plo_orthog
export plo_monitor
export get_win1
export get_win2
export try_blk1
export try_blk2
export try_diag
export calc_ovlp
export calc_dm
export calc_hamk
export calc_dos
export view_ovlp
export view_dm
export view_hamk
export view_dos

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
# vasp_adaptor   -> Adaptor support.
# vasp_init      -> Prepare vasp's input files.
# vasp_run       -> Execute vasp program.
# vasp_save      -> Backup vasp's output files.
# vasp_incar     -> Make essential input file (INCAR).
# vasp_kpoints   -> Make essential input file (KPOINTS).
# vasp_files     -> Check essential output files.
# vaspio_nband   -> Determine number of bands.
# vaspio_lattice -> Read lattice information.
# vaspio_valence -> Read number of valence electrons per sorts.
# vaspio_kmesh   -> Read kmesh.
# vaspio_tetra   -> Read tetrahedra.
# vaspio_eigen   -> Read eigenvalues.
# vaspio_projs   -> Read projectors.
# vaspio_fermi   -> Read fermi level.
# vaspio_charge  -> Read charge density.
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
export vaspio_nband
export vaspio_lattice
export vaspio_valence
export vaspio_kmesh
export vaspio_tetra
export vaspio_eigen
export vaspio_projs
export vaspio_fermi
export vaspio_charge

#
# tetra.jl
#
# Summary:
#
# Tools for the analytical tetrahedron method (adaptor).
#
# Members:
#
# TetraWeight  -> Struct for integration weights.
# bzint        -> Compute tetrahedron integrated weights.
# gauss_weight -> Compute integrated weights using Gaussian broadening.
# fermi_weight -> Compute integrated weights using Fermi-Dirac broadening.
# tetra_weight -> Compute integrated weights for a given tetrahedron.
# tetra_p_ek1  -> Blochl tetrahedron integration algorithm, case 1.
# tetra_p_ek12 -> Blochl tetrahedron integration algorithm, case 2.
# tetra_p_ek23 -> Blochl tetrahedron integration algorithm, case 3.
# tetra_p_ek34 -> Blochl tetrahedron integration algorithm, case 4.
# tetra_p_ek4  -> Blochl tetrahedron integration algorithm, case 5.
#
include("tetra.jl")
export TetraWeight
export bzint
export gauss_weight
export fermi_weight
export tetra_weight
export tetra_p_ek1
export tetra_p_ek12
export tetra_p_ek23
export tetra_p_ek34
export tetra_p_ek4

"""
    __init__()

This function would be executed immediately after the module is loaded at
runtime for the first time.
"""
function __init__() end

end
