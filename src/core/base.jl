#
# project : pansy
# source  : base.jl
# author  : Li Huang (lihuang.dmft@gmail.com)
# status  : unstable
# comment :
#
# last modified: 2021/01/13
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
    for i = 1:_i("nsite")
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
    for i = 1:_i("nsite")
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
    engine = _d("engine")
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
    engine = _d("engine")
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
    # directly. we have to check and process them carefully
    #
    plo_adaptor()

    #
    # A3: Output the processed Kohn-Sham data
    #
    # ok, now the Kohn-Sham data are ready. we would like to write them
    # to some specified files. maybe you want to see some key quantities
    #
    ir_adaptor()

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

To examine whether the dft runtime environment is ready
"""
function dft_init(it::IterInfo)
    # enter dft directory
    cd("dft")

    # choose suitable dft engine
    engine = _d("engine")
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

Execute the desired dft engine parallelly or sequentially
"""
function dft_run(it::IterInfo)
    # enter dft directory
    cd("dft")

    # choose suitable dft engine
    engine = _d("engine")
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

Backup the essential dft calculated results for next iterations
"""
function dft_save(it::IterInfo)
    # enter dft directory
    cd("dft")

    # choose suitable dft engine
    engine = _d("engine")
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
    dmft_init()
"""
function dmft_init() end

"""
    dmft_run()
"""
function dmft_run() end

"""
    dmft_save()
"""
function dmft_save() end

"""
    solver_init()
"""
function solver_init() end

"""
    solver_run()
"""
function solver_run() end

"""
    solver_save()
"""
function solver_save() end
