#
# Project : Pansy
# Source  : base.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/03/29
#

#
# Driver Functions
#

"""
    ready()

Examine whether all the conditions for DFT + DMFT calculations are ready.

See also: [`go`](@ref).
"""
function ready()
    # R1: Check the input files
    query_inps(get_d("engine"))

    # R2: Prepare the working directories
    make_trees()
end

"""
    go()

Dispatcher for DFT + DMFT calculations.

See also: [`ready`](@ref).
"""
function go()
    # Get calculation mode
    mode = get_m("mode")

    # Choose suitable computational driver
    @cswitch mode begin
        # One-shot DFT + DMFT calculations
        @case 1
            cycle1()
            break

        # Fully self-consistent DFT + DMFT calculations
        @case 2
            cycle2()
            break

        # To be implemented
        @default
            sorry()
            break
    end
end

"""
    final()

Finalize the DFT + DMFT calculations.
"""
function final()
    sorry()
end

"""
    cycle1()

Perform one-shot DFT + DMFT calculations. In other words, the charge
density won't be fed back to the DFT engine. The self-consistency is
only achieved at the DMFT level.

See also: [`cycle2`](@ref), [`go`](@ref).
"""
function cycle1()
    # C-1: Create IterInfo struct
    it = IterInfo()

    # C00: Create Logger struct
    lr = Logger(query_case())

#
# Initialization (C01-C04)
#
    prompt("ZEN", "Initialization")

#
# Remarks 1:
#
# We would like to perform two successive DFT runs if get_d("loptim") is
# true. The purpose of the first DFT run is to evaluate the fermi level.
# Then an energy window is determined. We will use this window to generate
# optimal projectors in the second DFT run.
#
# On the other hand, if get_d("loptim") is false, only the first DFT run
# is enough.
#

    # C01: Perform DFT calculation (for the first time)
    dft_run(it, lr)

#
# Remarks 2:
#
# We want better optimal projectors.
#
# In the previous DFT run, initial fermi level = 0 -> wrong energy
# window -> wrong optimial projectors. But at this point, the fermi
# level is updated, so we have to generate the optimal projectors
# again within this new window by doing addition DFT calculation.
#

    # C02: Perform DFT calculation (for the second time)
    if get_d("loptim")
        dft_run(it, lr)
    end

#
# Remarks 3:
#
# The key Kohn-Sham data inclue lattice structures, k-mesh and its weights,
# tetrahedra data, eigenvalues, raw projectors, and fermi level, etc. At
# first, the adaptor will read in these data from the output files of DFT
# engine. And then it will process the raw projectors (such as parsing,
# labeling, grouping, filtering, and rotatation). Finally, the adaptor will
# write down the processed data to some specified files using the IR format.
#

    # C03: To bridge the gap between DFT engine and DMFT engine by adaptor
    adaptor_run(it, lr)

    # C04: Prepare default self-energy functions
    sigma_core(lr, "reset")

#
# Remarks 4:
#
# Now everything is ready. We are going to solve the DMFT self-consistent
# equation iterately.
#

#
# Iterations (C05-C09)
#
    prompt("ZEN", "Iterations")

    for iter = 1:get_m("niter")
        prompt("ZEN", "Cycle $iter")

        # C05: Tackle with the double counting term
        sigma_core(lr, "dcount")
        exit(-1)

        # C06: Perform DMFT calculation with `dmft_mode` = 1
        dmft_run(it, lr, 1)

        # C07: Split and distribute the data (hybridization functions)
        sigma_core(lr, "split")

        # C08: Solve the quantum impurity problems
        solver_run(it, lr)

        # C09: Gather and combine the data (impurity self-functions)
        sigma_core(lr, "gather")

        # C10: Mixer for self-energy functions or hybridization functions
        mixer_core(lr)
    end

    # C98: Close Logger.log
    if isopen(lr.log)
        flush(lr.log)
        close(lr.log)
    end

    # C99: Close Logger.cycle
    if isopen(lr.cycle)
        flush(lr.cycle)
        close(lr.cycle)
    end
end

"""
    cycle2()

Perform fully self-consistent DFT + DMFT calculations. The self-consistency
is achieved at both DFT and DMFT levels.

See also: [`cycle1`](@ref), [`go`](@ref).
"""
function cycle2()
    sorry()
end

#
# Service Functions
#

"""
    monitor(force_exit::Bool = false)

Determine whether we need to terminate the Zen code.

See also: [`query_stop`](@ref).
"""
function monitor(force_exit::Bool = false)

#
# Remarks:
#
# In order to terminate the Zen code, the following two conditions
# should be fulfilled at the same time.
#
# 1. The argument `force_exit` is true.
#
# 2. The case.stop file exists (from query_stop()).
#

    if force_exit && query_stop()
        exit(-1)
    end
end

"""
    make_trees()

Prepare the working directories at advance.

See also: [`rm_trees`](@ref).
"""
function make_trees()

#
# Remarks:
#
# The working directories include dft, dmft1, dmft2, and impurity.i.
# If they exist already, it would be better to remove them at first.
#

    # Build an array for folders
    dir_list = ["dft", "dmft1", "dmft2"]
    for i = 1:get_i("nsite")
        push!(dir_list, "impurity.$i")
    end

    # Go through these folders, create them one by one.
    for i in eachindex(dir_list)
        dir = dir_list[i]
        if isdir(dir)
            rm(dir, force = true, recursive = true)
        end
        mkdir(dir)
    end
end

"""
    rm_trees()

Remove the working directories finally.

See also: [`make_trees`](@ref).
"""
function rm_trees()
    # Build an array for folders
    dir_list = ["dft", "dmft1", "dmft2"]
    for i = 1:get_i("nsite")
        push!(dir_list, "impurity.$i")
    end

    # Go through these folders, remove them one by one.
    for i in eachindex(dir_list)
        dir = dir_list[i]
        if isdir(dir)
            rm(dir, force = true, recursive = true)
        end
    end
end

"""
    dft_run(it::IterInfo, lr::Logger)

Simple driver for DFT engine. It performs three tasks: (1) Examine
the runtime environment for the DFT engine. (2) Launch the DFT engine.
(3) Backup the output files by DFT engine for next iterations.

Now only the VASP engine is supported. If you want to support the other
DFT engine, this function must be adapted.

See also: [`adaptor_run`](@ref), [`dmft_run`](@ref), [`solver_run`](@ref).
"""
function dft_run(it::IterInfo, lr::Logger)
    # Determine the chosen engine 
    engine = get_d("engine")

    # Print the log
    prompt("DFT")
    prompt(lr.log, engine)

    # Enter dft directory
    cd("dft")

    # Activate the chosen DFT engine
    @cswitch engine begin
        # For VASP
        @case "vasp"
            vasp_init(it)
            vasp_exec(it)
            vasp_save(it)
            break

        @default
            sorry()
            break
    end

    # Enter the parent directory
    cd("..")

    # Monitor the status
    monitor(true)
end

"""
    dmft_run(it::IterInfo, lr::Logger, dmft_mode::I64)

Simple driver for DMFT engine. It performs three tasks: (1) Examine
the runtime environment for the DMFT engine. (2) Launch the DMFT engine.
(3) Backup the output files by DMFT engine for next iterations.

See also: [`adaptor_run`](@ref), [`dft_run`](@ref), [`solver_run`](@ref).
"""
function dmft_run(it::IterInfo, lr::Logger, dmft_mode::I64)
    # Examine the argument `dmft_mode`
    @assert dmft_mode === 1 || dmft_mode === 2

    # Print the log
    prompt("DMFT")
    prompt(lr.log, "dmft$dmft_mode")

    # Enter dmft1 or dmft2 directory
    cd("dmft$dmft_mode")

    # Activate the chosen DMFT engine
    @cswitch dmft_mode begin
        # Solve the DMFT self-consistent equation
        @case 1
            # sorry()
            break

        # Generate DMFT correction for DFT charge density
        @case 2
            # sorry()
            break
    end

    # Enter the parent directory
    cd("..")

    # Monitor the status
    monitor(true)
end

"""
    solver_run(it::IterInfo, lr::Logger)

Simple driver for quantum impurity solvers. It performs three tasks: (1)
Examine the runtime environment for quantum impurity solver. (2) Launch
the quantum impurity solver. (3) Backup output files by quantum impurity
solver for next iterations.

Now only the `ct_hyb1`, `ct_hyb2`, `hub1`, and `norg` quantum impurity
solvers are supported. If you want to support the other quantum impurity
solvers, this function must be adapted.

See also: [`adaptor_run`](@ref), [`dft_run`](@ref), [`dmft_run`](@ref).
"""
function solver_run(it::IterInfo, lr::Logger)
    # Loop over each impurity site
    for i = 1:get_i("nsite")

        # Determine the chosen solver
        engine = get_s("engine")

        # Print the log
        prompt("Solvers")
        prompt(lr.log, engine)

        # Enter impurity.i directory
        cd("impurity.$i")

        # Activate the chosen quantum impurity solver
        @cswitch engine begin
            @case "ct_hyb1"
                #sorry()
                break

            @case "ct_hyb2"
                #sorry()
                break

            @case "hub1"
                #sorry()
                break

            @case "norg"
                #sorry()
                break

            @default
                sorry()
                break
        end

        # Enter the parent directory
        cd("..")

    end

    # Monitor the status
    monitor(true)
end

"""
    adaptor_run(it::IterInfo, lr::Logger)

Simple driver for the adaptor. It performs three tasks: (1) Initialize
the adaptor, to check whether the essential files exist. (2) Parse the
Kohn-Sham data output by the DFT engine, try to postprocess them, and
then transform them into IR format. (3) Backup the files by adaptor.

For the first task, only the VASP adaptor is supported. While for the
second task, only the PLO adaptor is supported. If you want to support
more adaptor, please adapt this function.

See also: [`dft_run`](@ref), [`dmft_run`](@ref), [`solver_run`](@ref).
"""
function adaptor_run(it::IterInfo, lr::Logger)
    # Enter dft directory
    cd("dft")

    #
    # A0: Create a dict named DFTData
    #
    # This dictionary is for storing the Kohn-Sham band structure and
    # related data. The key-value pairs would be inserted into this
    # dict dynamically.
    #
    DFTData = Dict{Symbol,Any}()

    #
    # A1: Parse the original Kohn-Sham data
    #
    # Choose suitable driver function according to DFT engine. The
    # Kohn-Sham data will be stored in the DFTData dict.
    #
    engine = get_d("engine")
    prompt("Adaptor")
    prompt(lr.log, "adaptor::$engine")
    @cswitch engine begin
        # For VASP
        @case "vasp"
            vasp_files()
            vasp_adaptor(DFTData)
            break

        @default
            sorry()
            break
    end

    #
    # A2: Process the original Kohn-Sham data
    #
    # Well, now we have the Kohn-Sham data. But they can not be used
    # directly. We have to check and process them carefully. Please
    # pay attention to that the DFTData dict will be modified in
    # the plo_adaptor() function. Here the parameter `debug` (= true)
    # means that we are going to calculating some physical quantities
    # additionally to check the correctness of the Kohn-Sham data.
    #
    projtype = get_d("projtype")
    prompt("Adaptor")
    prompt(lr.log, "adaptor::$projtype")
    @cswitch projtype begin
        # For projected local orbital scheme (Now we disable debug)
        @case "plo"
            plo_adaptor(DFTData, false)
            break

        # For maximally localized wannier function scheme
        @case "wannier"
            sorry()
            break

        @default
            sorry()
            break
    end

    #
    # A3: Output the processed Kohn-Sham data
    #
    # Ok, now the Kohn-Sham data are ready. We would like to write them
    # to some specified files with the IR format. Then these files will
    # be saved immediately.
    #
    prompt("Adaptor")
    prompt(lr.log, "adaptor::ir")
    ir_adaptor(DFTData)
    ir_save(it)

    #
    # A4: Clear the DFTData dict
    #
    empty!(DFTData)
    @assert isempty(DFTData)

    # Enter the parent directory
    cd("..")

    # Monitor the status
    monitor(true)
end

"""
    sigma_core(lr::Logger, task::String = "reset")

Simple driver for functions for processing the self-energy functions
and hybridization functions.

Now it supports four tasks: `reset`, `dcount`, `split`, `gather`.

See also: [`mixer_core`](@ref).
"""
function sigma_core(lr::Logger, task::String = "reset")
    # Check the given task
    @assert task in ("reset", "dcount", "split", "gather")

    # Print the log
    prompt("Sigma")
    prompt(lr.log, "sigma::$task")

    # Launch suitable subroutine
    @cswitch task begin
        # Generate default self-energy functions and store them
        @case "reset"
            sigma_reset()
            break

        # Calculate the double counting term and store it
        @case "dcount"
            sigma_dcount()
            break

        # Split the hybridization functions and store them
        @case "split"
            sigma_split()
            break

        # Collect impurity self-energy functions and combine them
        @case "gather"
            sigma_gather()
            break

        @default
            sorry()
            break
    end

    # Monitor the status
    monitor(true)
end

"""
    mixer_core(lr::Logger)

Simple driver for the mixer. It will try to mix the self-energy functions
or hybridization functions and generate a new one.

See also: [`sigma_core`](@ref).
"""
function mixer_core(lr::Logger)
    # Print the log
    prompt("Mixer")
    prompt(lr.log, "mixer")

    # Monitor the status
    monitor(true)
end
