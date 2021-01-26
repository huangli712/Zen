#
# Project : Pansy
# Source  : base.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
# Comment :
#
# Last modified: 2021/01/26
#

#
# Driver Functions
#

"""
    ready()

Examine whether all the conditions for DFT + DMFT calculations are ready.
"""
function ready()
    # S1: Check the input files
    query_inps()

    # S2: Prepare the working directories
    make_trees()
end

"""
    go()

Dispatcher for the DFT + DMFT calculations.
"""
function go()
    # Get calculation mode
    mode = get_m("mode")

    # Choose suitable computational driver
    @cswitch mode begin
        # One-shot DFT + DMFT calculations
        @case 1
            prompt("ZEN", "Cycling")
            cycle1()
            break

        # Fully self-consistent DFT + DMFT calculations
        @case 2
            prompt("ZEN", "Cycling")
            cycle2()
            break

        # To be implemented
        @default
            prompt("ZEN", "Nothing")
            sorry()
            break
    end
end

"""
    final()

Finalize the DFT + DMFT calculations.
"""
function final()
    # TODO
    sorry()
end

"""
    cycle1()

Perform one-shot DFT + DMFT calculations.
"""
function cycle1()
    # S-1: Create IterInfo
    it = IterInfo()

    # S00: Create Logger
    lr = Logger()
    lr.log = open(query_case() * ".log", "a")
    lr.cycle = open(query_case() * ".cycle", "a")

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

    # S01: Perform DFT calculation (for the first time)
    #
    # S01.1: Prepare and check essential files for the DFT engine
    dft_init(it)
    #
    # S01.2: Perform a self-consitent calculation at the DFT level
    dft_run(it)
    #
    # S01.3: Backup the output files of the DFT engine
    dft_save(it)
    #
    # S01.4: Monitor the status
    monitor(true)

    # We want better optimal projectors.
    # In the previous DFT run, initial fermi level = 0 -> wrong energy
    # window -> wrong optimial projectors. But at this point, the fermi
    # level is updated, so we have to generate the optimal projectors
    # again within this new window by carrying addition DFT calculation.
    if get_d("loptim")

        # S02: Perform DFT calculation (for the second time)
        #
        # S02.1: Prepare and check essential files for the DFT engine
        dft_init(it)
        #
        # S02.2: Perform a self-consitent calculation at the DFT level
        dft_run(it)
        #
        # S02.3: Backup the output files of the DFT engine
        dft_save(it)
        #
        # S02.4: Monitor the status
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

    # S03: To bridge the gap between DFT engine and DMFT engine by adaptor
    #
    # S03.1: Prepare and check essential files for the adaptor
    adaptor_init(it)
    #
    # S03.2: Launch the adaptor
    adaptor_run(it)
    #
    # S03.3: Backup the output files of the adaptor
    adaptor_save(it)
    #
    # S03.4: Monitor the status
    monitor(true)

#
# Remarks 3:
#
# Now everything is ready. We are going to solve the DMFT equation iterately.. 
#
    exit(-1)

    for iter = 1:get_m("niter")
        # S04: Perform DMFT calculation
        #
        # S04.1: Prepare and check essential files for the DMFT engine
        dmft_init(it)
        #
        # S04.2: Launch the DMFT engine (dmft1)
        dmft_run(it)
        #
        # S04.3: Backup the output files of the DMFT engine
        dmft_save(it)
        #
        # S04.4: Monitor the status
        monitor(true)

        # S05: Solve the quantum impurity problems
        #
        # S05.1: Prepare and check essential files for the quantum impurity solver
        solver_init(it)
        #
        # S05.2: Launch the quantum impurity solver
        solver_run(it)
        #
        # S05.3: Backup the output files of the quantum impurity solver
        solver_save(it)
        #
        # S05.4: Monitor the status
        monitor(true)
    end
end

"""
    cycle2()

Perform fully self-consistent DFT + DMFT calculations.
"""
function cycle2()
    # TODO
    sorry()
end

#
# Service Functions
#

"""
    monitor(force_exit::Bool = false)

Determine whether we need to terminate the Zen code.
"""
function monitor(force_exit::Bool = false)

#
# Remarks:
#
# In order to terminate the Zen code, the following two conditions
# should be fulfilled at the same time.
#
# 1. The argument force_exit is true.
#
# 2. The case.stop file exists.
#

    if force_exit && query_stop()
        exit(-1)
    end
end

"""
    make_trees()

Prepare the working directories at advance.
"""
function make_trees()

#
# Remarks:
#
# The working directories include dft, dmft1, dmft2, and impurity.i.
# If they exist already, it would be better to remove them at first.
#

    # For dft
    if isdir("dft")
        rm("dft", force = true, recursive = true)
    end
    mkdir("dft")

    # For dmft1
    if isdir("dmft1")
        rm("dmft1", force = true, recursive = true)
    end
    mkdir("dmft1")

    # For dmft2
    if isdir("dmft2")
        rm("dmft2", force = true, recursive = true)
    end
    mkdir("dmft2")

    # For impurity.i
    for i = 1:get_i("nsite")
        if isdir("impurity.$i")
            rm("impurity.$i", force = true, recursive = true)
        end
        mkdir("impurity.$i")
    end
end

"""
    rm_trees()

Remove the working directories finally.
"""
function rm_trees()
    # For dft
    if isdir("dft")
        rm("dft", force = true, recursive = true)
    end

    # For dmft1
    if isdir("dmft1")
        rm("dmft1", force = true, recursive = true)
    end

    # For dmft2
    if isdir("dmft2")
        rm("dmft2", force = true, recursive = true)
    end

    # For impurity.i
    for i = 1:get_i("nsite")
        if isdir("impurity.$i")
            rm("impurity.$i", force = true, recursive = true)
        end
    end
end

"""
    adaptor_init(it::IterInfo)

Initialize the adaptor, to check whether the essential files exist.
"""
function adaptor_init(it::IterInfo)
    # Enter dft directory
    cd("dft")

    # Choose suitable adaptor according to DFT engine
    engine = get_d("engine")
    @cswitch engine begin
        @case "vasp"
            prompt("Adaptor : VASP")
            println("  Init VASP Adaptor")
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

Parse the data output by DFT engine, postprocess them, and then transform
them into IR format.
"""
function adaptor_run(it::IterInfo)
    # Enter dft directory
    cd("dft")

    # Clear the DFTData dict
    for k in keys(DFTData)
        delete!(DFTData, k)
    end

    #
    # A1: Parse the original Kohn-Sham data
    #
    #
    # Choose suitable driver function according to DFT engine. The
    # Kohn-Sham data will be stored in the DFTData dict.
    #
    engine = get_d("engine")
    @cswitch engine begin
        @case "vasp"
            println("  Launch VASP Adaptor")
            vasp_adaptor()
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
    # the plo_adaptor() function. Here the parameter (debug = )true
    # means that we are going to seeing some interesting quantities
    # to check the correctness of the Kohn-Sham data.
    #
    println("  Launch PLO Adaptor")
    plo_adaptor(true)

    #
    # A3: Output the processed Kohn-Sham data
    #
    # Ok, now the Kohn-Sham data are ready. We would like to write them
    # to some specified files.
    #
    println("  Launch IR Adaptor")
    ir_adaptor()

    # Enter the parent directory
    cd("..")
end

"""
    adaptor_save(it::IterInfo)

Backup the output files by adaptor.
"""
function adaptor_save(it::IterInfo)
    # Enter dft directory
    cd("dft")

    # Save the essential files
    ir_save(it)

    # Enter the parent directory
    cd("..")
end

"""
    dft_init(it::IterInfo)

To examine the runtime environment for density functional theory engine.
"""
function dft_init(it::IterInfo)
    # Enter dft directory
    cd("dft")

    # Choose suitable DFT engine, then initialize it's input files
    engine = get_d("engine")
    @cswitch engine begin
        @case "vasp"
            prompt("DFT : VASP")
            println("  Init VASP")
            vasp_init(it)
            break

        @default
            prompt("DFT : Undef")
            sorry()
            break
    end

    # Enter the parent directory
    cd("..")
end

"""
    dft_run(it::IterInfo)

Launch the density functional theory engine.
"""
function dft_run(it::IterInfo)
    # Enter dft directory
    cd("dft")

    # Choose suitable DFT engine, then launch it
    engine = get_d("engine")
    @cswitch engine begin
        @case "vasp"
            println("  Launch VASP")
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
"""
function dft_save(it::IterInfo)
    # Enter dft directory
    cd("dft")

    # Choose suitable DFT engine, then backup some essential output files
    engine = get_d("engine")
    @cswitch engine begin
        @case "vasp"
            println("  Backup VASP's output data\n")
            vasp_save(it)
            break

        @default
            sorry()
            break
    end

    # Enter the parent directory
    cd("..")
end

"""
    dmft_init(it::IterInfo)

To examine the runtime environment for dynamical mean-field theory engine.
"""
function dmft_init(it::IterInfo)
    # Enter dmft1 directory
    cd("dmft1")

    # TODO

    # Enter the parent directory
    cd("..")
end

"""
    dmft_run(it::IterInfo)

Launch the dynamical mean-field theory engine.
"""
function dmft_run(it::IterInfo)
    # Enter dmft1 directory
    cd("dmft1")

    # TODO

    # Enter the parent directory
    cd("..")
end

"""
    dmft_save(it::IterInfo)

Backup the output files by dynamical mean-field theory engine
for next iterations.
"""
function dmft_save(it::IterInfo)
    # Enter dmft1 directory
    cd("dmft1")

    # TODO

    # Enter the parent directory
    cd("..")
end

"""
    solver_init(it::IterInfo)

To examine the runtime environment for quantum impurity solver.
"""
function solver_init(it::IterInfo)
    # Loop over each impurity site
    for i = 1:get_i("nsite")

        # Enter impurity.i directory
        cd("impurity.$i")

        # Choose suitable quantum impurity solver
        engine = get_s("engine")
        @cswitch engine begin
            @case "ct_hyb1"
                sorry()
                break

            @case "ct_hyb2"
                sorry()
                break

            @case "hub1"
                sorry()
                break

            @case "norg"
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
                sorry()
                break

            @case "ct_hyb2"
                sorry()
                break

            @case "hub1"
                sorry()
                break

            @case "norg"
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
                sorry()
                break

            @case "ct_hyb2"
                sorry()
                break

            @case "hub1"
                sorry()
                break

            @case "norg"
                sorry()
                break
        end

        # Enter the parent directory
        cd("..")

    end
end
