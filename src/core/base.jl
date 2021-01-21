#
# Project : Pansy
# Source  : base.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
# Comment :
#
# Last modified: 2021/01/21
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
# Remarks:
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
    prompt("DFT", "BEGIN")
    #
    # S01.1: Prepare and check essential files for the DFT engine
    dft_init(it)
    #
    # S01.2: Perform a self-consitent calculation at the DFT level
    dft_run(it)
    #
    # S01.3: Backup the output files of the DFT engine
    dft_save(it)

    # We want better optimal projectors.
    # In the previous DFT run, initial fermi level = 0 -> wrong energy
    # window -> wrong optimial projectors. But at this point, the fermi
    # level is updated, so we have to generate the optimal projectors
    # again within this new window by carrying addition DFT calculation.
    if get_d("loptim")

        # S02: Perform DFT calculation (for the second time)
        prompt("ZEN", "DFT")
        #
        # S02.1: Prepare and check essential files for the DFT engine
        println("Initialize everything needed by the DFT engine")
        dft_init(it)
        #
        # S02.2: Perform a self-consitent calculation at the DFT level
        println("Launch the DFT engine")
        dft_run(it)
        #
        # S02.3: backup the output files of the DFT engine
        println("Save the output files\n")
        dft_save(it)

    end

#
# remarks:
#
# the key Kohn-Sham data inclue lattice structures, k-mesh and its weights,
# tetrahedra data, eigenvalues, raw projectors, and fermi level, etc. at
# first the adaptor will read in these data from the output files of dft
# engine. and then it will process the raw projectors (parse, label, group,
# filter, and rotate). finally, the adaptor will write down the processed
# data to some specified files within the IR format.
#

    # S07: To bridge the gap between dft engine and dmft engine by adaptor
    prompt("ZEN", "ADAPTOR")
    #
    # S07.1: prepare and check essential files for the adaptor
    println("Initialize everything needed by the adaptor")
    adaptor_init(it)
    #
    # S07.2: launch the adaptor. it will read the Kohn-Sham data from the
    # dft engine, postprocess them, and then write them to external files
    # with the IR format
    println("Launch the adaptor")
    adaptor_run(it)
    #
    # S07.3: backup the output files of the adaptor
    println("Save the output files\n")
    adaptor_save(it)

    for iter = 1:get_m("niter")
        prompt("ZEN", "ITER : $iter")
        prompt("ZEN", "DMFT1")
        prompt("ZEN", "SOLVER")
    end
end

"""
    cycle2()

Perform fully self-consistent DFT + DMFT calculations.
"""
function cycle2()
    sorry()
end

#
# Service Functions
#

"""
    monitor()
"""
function monitor() end

"""
    make_trees()

Prepare the working directories at advance
"""
function make_trees()
    # the working directories include dft, dmft1, dmft2, and impurity.i
    # if they exist already, we have to remove them at first
    #
    # for dft
    if isdir("dft")
        rm("dft", force = true, recursive = true)
    end
    mkdir("dft")

    # for dmft1
    if isdir("dmft1")
        rm("dmft1", force = true, recursive = true)
    end
    mkdir("dmft1")

    # for dmft2
    if isdir("dmft2")
        rm("dmft2", force = true, recursive = true)
    end
    mkdir("dmft2")

    # for impurity.i
    for i = 1:get_i("nsite")
        if isdir("impurity.$i")
            rm("impurity.$i", force = true, recursive = true)
        end
        mkdir("impurity.$i")
    end
end

"""
    rm_trees()

Remove the working directories finally
"""
function rm_trees()
    # for dft
    if isdir("dft")
        rm("dft", force = true, recursive = true)
    end

    # for dmft1
    if isdir("dmft1")
        rm("dmft1", force = true, recursive = true)
    end

    # for dmft2
    if isdir("dmft2")
        rm("dmft2", force = true, recursive = true)
    end

    # for impurity.i
    for i = 1:get_i("nsite")
        if isdir("impurity.$i")
            rm("impurity.$i", force = true, recursive = true)
        end
    end
end

"""
    adaptor_init(it::IterInfo)

Initialize the adaptor, to check whether the essential files exist
"""
function adaptor_init(it::IterInfo)
    # enter dft directory
    cd("dft")

    # choose suitable driver function according to dft engine
    engine = get_d("engine")
    @cswitch engine begin
        @case "vasp"
            vasp_files()
            break

        @default
            sorry()
            break
    end

    # enter the parent directory
    cd("..")
end

"""
    adaptor_run(it::IterInfo)

Parse the data output by dft engine, postprocess them, and then transform
them into IR format
"""
function adaptor_run(it::IterInfo)
    # enter dft directory
    cd("dft")

    #
    # A1: Parse the original Kohn-Sham data
    #
    #
    # choose suitable driver function according to dft engine
    # the Kohn-Sham data will be stored in the KohnShamData dict
    #
    engine = get_d("engine")
    @cswitch engine begin
        @case "vasp"
            vasp_adaptor()
            break

        @default
            sorry()
            break
    end

    #
    # A2: Process the original Kohn-Sham data
    #
    # well, now we have the Kohn-Sham data. but they can not be used
    # directly. we have to check and process them carefully. please
    # pay attention to that the KohnShamData dict will be modified in
    # the plo_adaptor() function
    #
    plo_adaptor()

    #
    # A3: Output the processed Kohn-Sham data
    #
    # ok, now the Kohn-Sham data are ready. we would like to write them
    # to some specified files. here the parameter (view = )true means
    # that we are going to seeing some interesting quantities to check
    # the correctness of the Kohn-Sham data
    #
    ir_adaptor(true)

    # enter the parent directory
    cd("..")
end

"""
    adaptor_save(it::IterInfo)

Backup the output files by adaptor
"""
function adaptor_save(it::IterInfo)
    # enter dft directory
    cd("dft")

    # TODO

    # enter the parent directory
    cd("..")
end

"""
    dft_init(it::IterInfo)

To examine whether the runtime environment for density functional theory
engine is ready
"""
function dft_init(it::IterInfo)
    # enter dft directory
    cd("dft")

    # choose suitable dft engine
    engine = get_d("engine")
    @cswitch engine begin
        @case "vasp"
            vasp_init(it)
            break

        @default
            sorry()
            break
    end

    # enter the parent directory
    cd("..")
end

"""
    dft_run(it::IterInfo)

Launch the density functional theory engine
"""
function dft_run(it::IterInfo)
    # enter dft directory
    cd("dft")

    # choose suitable dft engine
    engine = get_d("engine")
    @cswitch engine begin
        @case "vasp"
            vasp_run(it)
            break

        @default
            sorry()
            break
    end

    # enter the parent directory
    cd("..")
end

"""
    dft_save(it::IterInfo)

Backup the output files by density functional theory engine
for next iterations
"""
function dft_save(it::IterInfo)
    # enter dft directory
    cd("dft")

    # choose suitable dft engine
    engine = get_d("engine")
    @cswitch engine begin
        @case "vasp"
            vasp_save(it)
            break

        @default
            sorry()
            break
    end

    # enter the parent directory
    cd("..")
end

"""
    dmft_init(it::IterInfo)

To examine whether the runtime environment for dynamical mean-field theory
engine is ready
"""
function dmft_init(it::IterInfo) end

"""
    dmft_run(it::IterInfo)

Launch the dynamical mean-field theory engine
"""
function dmft_run(it::IterInfo) end

"""
    dmft_save(it::IterInfo)

Backup the output files by dynamical mean-field theory engine
for next iterations
"""
function dmft_save(it::IterInfo) end

"""
    solver_init(it::IterInfo)

To examine whether the runtime environment for quantum impurity solver
is ready
"""
function solver_init(it::IterInfo) end

"""
    solver_run(it::IterInfo)

Launch the quantum impurity solver
"""
function solver_run(it::IterInfo) end

"""
    solver_save(it::IterInfo)

Backup the output files by quantum impurity solver for next iterations
"""
function solver_save(it::IterInfo) end
