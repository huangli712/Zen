#
# project : pansy
# source  : base.jl
# author  : Li Huang (lihuang.dmft@gmail.com)
# status  : unstable
# comment :
#
# last modified: 2021/01/19
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
    println("Ensure essential input files")
    query_inps()

    # S2: Prepare the working directories
    println("Create working directories")
    make_trees()
end

"""
    go()
"""
function go()
    if get_m("mode") === 1
        cycle1()
    else
        cycle2()
    end
end

function final() end
function cycle1() end
function cycle2() end

#
# Service Functions
#

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
