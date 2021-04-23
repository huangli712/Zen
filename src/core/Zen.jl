#
# Project : Pansy
# Source  : Zen.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/04/23
#

"""
    Zen

Zen is a modern DFT + DMFT computation framework. It can be used to
study the correlated electronic structures of a wide range of strongly
correlated materials. Now this framework is under heavy development.
**PLEASE USE IT AT YOUR OWN RISK**.

Zen supports the following DFT backends:
* `VASP`

Zen supports the following schemes for defining local orbitals:
* `PLO`
* `WANNIER`

Zen supports the following quantum impurity solvers:
* `CT-HYB`
* `HIA`
* `NORG`

For more details about how to install and use the Zen framework, please
visit the following website:

`http://doc.zen`
"""
module Zen

#
# Using standard libraries
#

#
# Remarks 1:
#
# The TOML.jl package is included in the standard library since v1.6.
# So please upgrade your julia environment if it is outdated. We need
# this package to parse the configuration file (TOML format).
#
# Remarks 2:
#
# Here we import `libm` explicitly to provide a callable interface for
# the `erf` function. See util.jl/erf() for more details.
#
# Remarks 3:
#
# Actually, `Reexport` is not included in the standard library. It is
# a third-party library.
#

using TOML
using LinearAlgebra
using Distributed
using Dates
using Base.Math: libm
using Reexport
@reexport using Printf

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
# authors     -> Print the authors of Zen to screen.
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
export authors

#
# util.jl
#
# Summary:
#
# To provide some useful utility macros and functions. They can be used
# to colorize the output strings, query the environments, and parse the
# input strings, etc.
#
# Members:
#
# @cswitch      -> C-style switch.
# @ps1          -> Wrapper for printstyled function.
# @ps2          -> Another wrapper for printstyled function.
# require       -> Check julia envirnoment.
# setup_args    -> Setup ARGS manually.
# query_args    -> Query program's arguments.
# query_case    -> Query case (job's name).
# query_inps    -> Query input files.
# query_stop    -> Query case.stop file.
# query_home    -> Query home directory of Zen.
# query_core    -> Query core directory of Zen (where is Zen.jl).
# query_dft     -> Query home directory of DFT engine.
# query_dmft    -> Query home directory of DMFT engine.
# query_solver  -> Query home directory of quantum impurity solvers.
# welcome       -> Print welcome message.
# overview      -> Print overview of Zen.
# goodbye       -> Say goodbye.
# sorry         -> Say sorry.
# prompt        -> Print some messages or logs to the output devices.
# line_to_array -> Convert a line to a string array.
# line_to_cmplx -> Convert a line to a cmplx number.
# erf           -> Gauss error function.
# subscript     -> Convert a number to subscript.
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
export query_home
export query_core
export query_dft
export query_dmft
export query_solver
export welcome
export overview
export goodbye
export sorry
export prompt
export line_to_array
export line_to_cmplx
export erf
export subscript

#
# tetra.jl
#
# Summary:
#
# Implementation of the analytical tetrahedron method.
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

#
# types.jl
#
# Summary:
#
# Define some dicts and structs, which are used to store the config
# parameters or represent some essential data structures.
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
# Mapping  -> Struct for mapping between impurity problems and projectors.
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
export Mapping
export PrTrait
export PrGroup
export PrWindow

#
# config.jl
#
# Summary:
#
# To extract, parse, verify, and print the configuration parameters.
# They are stored in external files (*.toml) or dictionaries.
#
# Members:
#
# setup    -> Setup parameters.
# exhibit  -> Display parameters for reference.
# inp_toml -> Parse case.toml, return raw configuration information.
# rev_dict -> Update dicts for configuration parameters.
# chk_dict -> Check dicts for configuration parameters.
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
# quantum impurity solvers, Kohn-Sham adaptor, self-energy engine, and
# mixer engine. The DFT + DMFT iteration (one-shot mode or charge fully
# self-consistent mode) is also implemented in this file.
#
# Members:
#
# ready       -> Prepare runtime environment for DFT + DMFT calculations.
# go          -> Dispatcher for DFT + DMFT calculations.
# final       -> Finalize the DFT + DMFT calculations.
# cycle0      -> Perform DFT calculations only (for test purpose).
# cycle1      -> Perform DFT + DMFT calculations (one-shot mode).
# cycle2      -> Perform DFT + DMFT calculations (fully self-consistent mode).
# monitor     -> Monitor the DFT + DMFT calculations.
# make_trees  -> Make working directories.
# rm_trees    -> Remove working directories.
# dft_run     -> Driver for DFT engine.
# dmft_run    -> Driver for DMFT engine.
# solver_run  -> Driver for quantum impurity solvers.
# adaptor_run -> Driver for Kohn-Sham adaptor.
# sigma_core  -> Driver for self-energy engine.
# mixer_core  -> Driver for mixer engine.
#
include("base.jl")
#
export ready
export go
export final
export cycle0
export cycle1
export cycle2
export monitor
export make_trees
export rm_trees
export dft_run
export dmft_run
export solver_run
export adaptor_run
export sigma_core
export mixer_core

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
# vasp_exec      -> Execute vasp program.
# vasp_save      -> Backup vasp's output files.
# vasp_incar     -> Generate essential input file (INCAR).
# vasp_kpoints   -> Generate essential input file (KPOINTS).
# vasp_files     -> Check essential output files.
# vaspio_nband   -> Determine number of bands.
# vaspio_valence -> Read number of valence electrons for each sort.
# vaspio_procar  -> Read PROCAR file.
# vaspio_lattice -> Read lattice information.
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
export vasp_exec
export vasp_save
export vasp_incar
export vasp_kpoints
export vasp_files
export vaspio_nband
export vaspio_valence
export vaspio_procar
export vaspio_lattice
export vaspio_kmesh
export vaspio_tetra
export vaspio_eigen
export vaspio_projs
export vaspio_fermi
export vaspio_charge

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
# plo_map     -> Create connection between projectors and impurity problems.
# plo_fermi   -> Calibrate Kohn-Sham eigenvalues with respect to fermi level.
# plo_group   -> Setup groups of projectors.
# plo_window  -> Setup band windows of projectors.
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
# calc_hamk   -> Calculate local hamiltonian or full hamiltonian.
# calc_dos    -> Calculate density of states.
# view_ovlp   -> Show overlap matrix for debug.
# view_dm     -> Show density matrix for debug.
# view_hamk   -> Show local hamiltonian for debug.
# view_dos    -> Show density of states for debug.
#
include("plo.jl")
#
export plo_adaptor
export plo_map
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
# irio_groups  -> Write PrGroup.
# irio_windows -> Write PrWindow.
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
export irio_groups
export irio_windows
export irio_lattice
export irio_kmesh
export irio_tetra
export irio_eigen
export irio_projs
export irio_fermi
export irio_charge

#
# dmft.jl
#
# Summary:
#
# Wrapper for dynamical mean-field theory engine.
#
# Members:
#
# dmft_init -> Prepare input files for the DMFT engine.
# dmft_exec -> Execute the DMFT engine.
# dmft_save -> Backup output files for the DMFT engine.
#
include("dmft.jl")
#
export dmft_init
export dmft_exec
export dmft_save

#
# solver.jl
#
# Summary:
#
# Wrapper for various quantum impurity solvers. Now only the CT-HYB1,
# CT-HYB2, HIA, and NORG quantum impurity solvers are supported.
#
# Members:
#
# s_qmc1_init -> Prepare input files for the CT-HYB1 impurity solver.
# s_qmc1_exec -> Execute the CT-HYB1 impurity solver.
# s_qmc1_save -> Backup output files for the CT-HYB1 impurity solver.
# s_qmc2_init -> Prepare input files for the CT-HYB2 impurity solver.
# s_qmc2_exec -> Execute the CT-HYB2 impurity solver.
# s_qmc2_save -> Backup output files for the CT-HYB2 impurity solver.
# s_hub1_init -> Prepare input files for the HIA impurity solver.
# s_hub1_exec -> Execute the HIA impurity solver.
# s_hub1_save -> Backup output files for the HIA impurity solver.
# s_norg_init -> Prepare input files for the NORG impurity solver.
# s_norg_exec -> Execute the NORG impurity solver.
# s_norg_save -> Backup output files for the NORG impurity solver.
#
include("solver.jl")
#
export s_qmc1_init
export s_qmc1_exec
export s_qmc1_save
export s_qmc2_init
export s_qmc2_exec
export s_qmc2_save
export s_hub1_init
export s_hub1_exec
export s_hub1_save
export s_norg_init
export s_norg_exec
export s_norg_save

#
# sigma.jl
#
# Summary:
#
# Tools for tackling with the self-energy functions.
#
# Members:
#
# sigma_reset  -> Create initial self-energy functions.
# sigma_dcount -> Calculate double counting term.
# sigma_split  -> Split the self-energy functions or hybridization functions.
# sigma_gather -> Gather and combine the self-energy functions.
# cal_dc_fll   -> Fully localized limit scheme for double counting term.
# cal_dc_amf   -> Around mean-field scheme for double counting term.
# cal_dc_exact -> Exact double counting scheme.
#
include("sigma.jl")
#
export sigma_reset
export sigma_dcount
export sigma_split
export sigma_gather
export cal_dc_fll
export cal_dc_amf
export cal_dc_exact

"""
    __init__()

This function would be executed immediately after the module is loaded
at runtime for the first time.

Here, we will try to precompile the whole Zen package to reduce latency
and speed up the later calculations.
"""
function __init__()
    prompt("ZEN", "Loading...")

    # Get an array of the names exported by the `Zen` module
    nl = names(Zen)

    # Go through each name
    cf = 0 # Counter
    for i in eachindex(nl)
        # Please pay attention to that nl[i] is a Symbol, we need to convert it
        # into string and function, respectively.
        str = string(nl[i])
        fun = eval(nl[i])

        # For methods only (macros must be excluded)
        if fun isa Function && !startswith(str, "@")
            # Increase the counter
            cf = cf + 1

            # Extract the signature of the function
            # Actually, `types` is a Core.SimpleVector.
            types = typeof(fun).name.mt.defs.sig.types

            # Convert `types` from SimpleVector into Tuple
            # If length(types) is 1, the method is without arguments.
            T = ()
            if length(types) > 1
                T = tuple(types[2:end]...)
            end

            # Precompile them one by one
            # println(i, " -> ", str, " -> ", length(types), " -> ", T)
            precompile(fun, T)
            @printf("Function %15s (#%3i) is compiled.\r", str, cf)
        end
    end

    prompt("ZEN", "Well, Zen is compiled and loaded ($cf functions).")
    prompt("ZEN", "We are ready to go!")
    println()
end

end # END OF MODULE
