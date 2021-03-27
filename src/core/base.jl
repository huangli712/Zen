#
# Project : Pansy
# Source  : base.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/03/27
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
end

"""
    cycle1()

Perform one-shot DFT + DMFT calculations. In other words, the charge
density won't be fed back to the DFT engine. The self-consistency is
only achieved at the DMFT level.

See also: [`cycle2`](@ref).
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
    prompt("DFT")
    #
    begin
        # C01.1: Prepare and check essential files for the DFT engine
        dft_init(it, lr)
        #
        # C01.2: Perform a self-consitent calculation at the DFT level
        dft_run(it)
        #
        # C01.3: Backup the output files of the DFT engine
        dft_save(it)
        #
        # C01.4: Monitor the status
        monitor(true)
    end

    # We want better optimal projectors.
    # In the previous DFT run, initial fermi level = 0 -> wrong energy
    # window -> wrong optimial projectors. But at this point, the fermi
    # level is updated, so we have to generate the optimal projectors
    # again within this new window by carrying addition DFT calculation.
    if get_d("loptim")

        # C02: Perform DFT calculation (for the second time)
        prompt("DFT")
        #
        # C02.1: Prepare and check essential files for the DFT engine
        dft_init(it, lr)
        #
        # C02.2: Perform a self-consitent calculation at the DFT level
        dft_run(it)
        #
        # C02.3: Backup the output files of the DFT engine
        dft_save(it)
        #
        # C02.4: Monitor the status
        monitor(true)

    end

#
# Remarks 2:
#
# The key Kohn-Sham data inclue lattice structures, k-mesh and its weights,
# tetrahedra data, eigenvalues, raw projectors, and fermi level, etc. At
# first, the adaptor will read in these data from the output files of DFT
# engine. And then it will process the raw projectors (such as parsing,
# labeling, grouping, filtering, and rotatation). Finally, the adaptor will
# write down the processed data to some specified files using the IR format.
#

    # C03: To bridge the gap between DFT engine and DMFT engine by adaptor
    prompt("Adaptor")
    #
    # C03.1: Prepare and check essential files for the adaptor
    adaptor_init(it, lr)
    #
    # C03.2: Launch the adaptor
    adaptor_run(it)
    #
    # C03.3: Backup the output files of the adaptor
    adaptor_save(it)
    #
    # C03.4: Monitor the status
    monitor(true)

    # C04: Prepare default self-energy functions
    prompt("Sigma")
    #
    # C04.1: Generate default self-energy functions and store them
    sigma_reset(lr)
    #
    # C04.2: Monitor the status
    monitor(true)

#
# Remarks 3:
#
# Now everything is ready. We are going to solve the DMFT equation iterately.
#

#
# Iterations (C05-C09)
#
    prompt("ZEN", "Iterations")

    for iter = 1:get_m("niter")
        prompt("ZEN", "Cycle $iter")

        # C05: Tackle with the double counting term
        prompt("Sigma")
        #
        # C05.1: Calculate the double counting term and store it
        sigma_dcount(lr)
        #
        # C05.2: Monitor the status
        monitor(true)

        # C06: Perform DMFT calculation
        prompt("DMFT")
        #
        # C06.1: Prepare and check essential files for the DMFT engine (dmft1)
        dmft_init(it, lr, 1)
        #
        # C06.2: Launch the DMFT engine (dmft1)
        dmft_run(it, 1)
        #
        # C06.3: Backup the output files of the DMFT engine (dmft1)
        dmft_save(it, 1)
        #
        # C06.4: Monitor the status
        monitor(true)

        # C07: Split and distribute the data
        prompt("Sigma")
        #
        # C07.1: Split the hybridization functions and store them
        sigma_split(lr)
        #
        # C07.2: Monitor the status
        monitor(true)

        # C08: Solve the quantum impurity problems
        prompt("Solvers")
        #
        # C08.1: Prepare and check essential files for the quantum impurity solver
        solver_init(it, lr)
        #
        # C08.2: Launch the quantum impurity solver
        solver_run(it)
        #
        # C08.3: Backup the output files of the quantum impurity solver
        solver_save(it)
        #
        # C08.4: Monitor the status
        monitor(true)

        # C09: Gather and combine the data
        prompt("Sigma")
        #
        # C09.1: Collect impurity self-energy functions and combine them
        sigma_gather(lr)
        #
        # C09.2: Monitor the status
        monitor(true)

        #
        # C09: Mixer
        #
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

See also: [`cycle1`](@ref).
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

function adaptor_core()
end

"""
    adaptor_init(it::IterInfo, lr::Logger)

Initialize the adaptor, to check whether the essential files exist.

See also: [`adaptor_run`](@ref), [`adaptor_save`](@ref).
"""
function adaptor_init(it::IterInfo, lr::Logger)
    # Enter dft directory
    cd("dft")

    # Choose suitable adaptor according to DFT engine
    engine = get_d("engine")
    prompt(lr.log, "adaptor")
    @cswitch engine begin
        @case "vasp"
            vasp_files()
            break

        @default
            sorry()
            break
    end

    # Enter the parent directory
    cd("..")
end

"""
    adaptor_run(it::IterInfo)

Parse the data output by the DFT engine, try to postprocess them, and then
transform them into IR format.

See also: [`adaptor_init`](@ref), [`adaptor_save`](@ref).
"""
function adaptor_run(it::IterInfo)
    # Enter dft directory
    cd("dft")

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
    @cswitch engine begin
        @case "vasp"
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
    @cswitch projtype begin
        # Now we disable debug
        @case "plo"
            plo_adaptor(DFTData, false)
            break

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
    # to some specified files with the IR format.
    #
    ir_adaptor(DFTData)

    #
    # A4: Clear the DFTData dict
    #
    empty!(DFTData)
    @assert isempty(DFTData)

    # Enter the parent directory
    cd("..")
end

"""
    adaptor_save(it::IterInfo)

Backup the output files by adaptor.

See also: [`adaptor_init`](@ref), [`adaptor_run`](@ref).
"""
function adaptor_save(it::IterInfo)
    # Enter dft directory
    cd("dft")

    # Save the essential files
    ir_save(it)

    # Enter the parent directory
    cd("..")
end

function dft_core()
end

"""
    dft_init(it::IterInfo, lr::Logger)

To examine the runtime environment for density functional theory engine.

See also: [`dft_run`](@ref), [`dft_save`](@ref).
"""
function dft_init(it::IterInfo, lr::Logger)
    # Enter dft directory
    cd("dft")

    # Choose suitable DFT engine, then initialize it's input files.
    engine = get_d("engine")
    prompt(lr.log, engine)
    @cswitch engine begin
        @case "vasp"
            vasp_init(it)
            break

        @default
            sorry()
            break
    end

    # Enter the parent directory
    cd("..")
end

"""
    dft_run(it::IterInfo)

Launch the density functional theory engine.

See also: [`dft_init`](@ref), [`dft_save`](@ref).
"""
function dft_run(it::IterInfo)
    # Enter dft directory
    cd("dft")

    # Choose suitable DFT engine, then launch it.
    engine = get_d("engine")
    @cswitch engine begin
        @case "vasp"
            vasp_run(it)
            break

        @default
            sorry()
            break
    end

    # Enter the parent directory
    cd("..")
end

"""
    dft_save(it::IterInfo)

Backup the output files by density functional theory engine
for next iterations.

See also: [`dft_init`](@ref), [`dft_run`](@ref).
"""
function dft_save(it::IterInfo)
    # Enter dft directory
    cd("dft")

    # Choose suitable DFT engine, then backup some essential output files.
    engine = get_d("engine")
    @cswitch engine begin
        @case "vasp"
            vasp_save(it)
            break

        @default
            sorry()
            break
    end

    # Enter the parent directory
    cd("..")
end

function dmft_core()
end

"""
    dmft_init(it::IterInfo, lr::Logger, dmft_mode::I64)

To examine the runtime environment for dynamical mean-field theory engine.

See also [`dmft_run`](@ref), [`dmft_save`])(@ref).
"""
function dmft_init(it::IterInfo, lr::Logger, dmft_mode::I64)
    # Examine the argument `dmft_mode`
    @assert dmft_mode === 1 || dmft_mode === 2

    # Solve the DMFT self-consistent equation
    if dmft_mode === 1
        # Enter dmft1 directory
        cd("dmft1")

        # TODO
        prompt(lr.log, "dmft1")

        # Enter the parent directory
        cd("..")

    # Generate DMFT correction to DFT
    else
        # Enter dmft2 directory
        cd("dmft2")

        # TODO
        prompt(lr.log, "dmft2")

        # Enter the parent directory
        cd("..")
    end
end

"""
    dmft_run(it::IterInfo, dmft_mode::I64)

Launch the dynamical mean-field theory engine.

See also [`dmft_init`](@ref), [`dmft_save`])(@ref).
"""
function dmft_run(it::IterInfo, dmft_mode::I64)
    # Examine the argument `dmft_mode`
    @assert dmft_mode === 1 || dmft_mode === 2

    # Solve the DMFT self-consistent equation
    if dmft_mode === 1
        # Enter dmft1 directory
        cd("dmft1")

        # TODO

        # Enter the parent directory
        cd("..")

    # Generate DMFT correction to DFT
    else
        # Enter dmft2 directory
        cd("dmft2")

        # TODO

        # Enter the parent directory
        cd("..")
    end
end

"""
    dmft_save(it::IterInfo, dmft_mode::I64)

Backup the output files by dynamical mean-field theory engine
for next iterations.

See also [`dmft_init`](@ref), [`dmft_run`])(@ref).
"""
function dmft_save(it::IterInfo, dmft_mode::I64)
    # Examine the argument `dmft_mode`
    @assert dmft_mode === 1 || dmft_mode === 2

    # Solve the DMFT self-consistent equation
    if dmft_mode === 1
        # Enter dmft1 directory
        cd("dmft1")

        # TODO

        # Enter the parent directory
        cd("..")

    # Generate DMFT correction to DFT
    else
        # Enter dmft2 directory
        cd("dmft2")

        # TODO

        # Enter the parent directory
        cd("..")
    end
end

function solver_init()
end

"""
    solver_init(it::IterInfo, lr::Logger)

To examine the runtime environment for quantum impurity solver.

See also: [`solver_run`](@ref), [`solver_save`](@ref).
"""
function solver_init(it::IterInfo, lr::Logger)
    # Loop over each impurity site
    for i = 1:get_i("nsite")

        # Enter impurity.i directory
        cd("impurity.$i")

        # Choose suitable quantum impurity solver
        engine = get_s("engine")
        @cswitch engine begin
            @case "ct_hyb1"
                prompt(lr.log, "ct_hyb1")
                break

            @case "ct_hyb2"
                prompt(lr.log, "ct_hyb2")
                break

            @case "hub1"
                prompt(lr.log, "hub1")
                break

            @case "norg"
                prompt(lr.log, "norg")
                break

            @default
                sorry()
                break
        end

        # Enter the parent directory
        cd("..")

    end
end

"""
    solver_run(it::IterInfo)

Launch the quantum impurity solver.

See also: [`solver_init`](@ref), [`solver_save`](@ref).
"""
function solver_run(it::IterInfo)
    # Loop over each impurity site
    for i = 1:get_i("nsite")

        # Enter impurity.i directory
        cd("impurity.$i")

        # Choose suitable quantum impurity solver
        engine = get_s("engine")
        @cswitch engine begin
            @case "ct_hyb1"
                break

            @case "ct_hyb2"
                break

            @case "hub1"
                break

            @case "norg"
                break

            @default
                sorry()
                break
        end

        # Enter the parent directory
        cd("..")

    end
end

"""
    solver_save(it::IterInfo)

Backup the output files by quantum impurity solver for next iterations.

See also: [`solver_init`](@ref), [`solver_run`](@ref).
"""
function solver_save(it::IterInfo)
    # Loop over each impurity site
    for i = 1:get_i("nsite")

        # Enter impurity.i directory
        cd("impurity.$i")

        # Choose suitable quantum impurity solver
        engine = get_s("engine")
        @cswitch engine begin
            @case "ct_hyb1"
                break

            @case "ct_hyb2"
                break

            @case "hub1"
                break

            @case "norg"
                break

            @default
                sorry()
                break
        end

        # Enter the parent directory
        cd("..")

    end
end
