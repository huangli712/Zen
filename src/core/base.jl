#
# project : pansy
# source  : base.jl
# author  : Li Huang (lihuang.dmft@gmail.com)
# status  : unstable
# comment :
#
# last modified: 2020/12/30
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
    adaptor_init(it::IterInfo)

Initialize the adaptor, to check whether the key files exist 
"""
function adaptor_init(it::IterInfo)
    # enter dft directory
    cd("dft")

    # choose suitable dft engine
    engine = _d("engine")
    @cswitch engine begin
        @case "vasp"
            vasp_files(pwd())
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

Parse the data output by dft engine, and transform them into IR format
"""
function adaptor_run(it::IterInfo)
    # enter dft directory
    cd("dft")

    # choose suitable dft engine
    engine = _d("engine")
    @cswitch engine begin
        @case "vasp"
            latt = vaspio_lattice(pwd())
            kmesh, weight = vaspio_kmesh(pwd())
            volt, itet = vaspio_tetra(pwd())
            enk, occupy = vaspio_eigen(pwd())
            PT, PG, chipsi = vaspio_projs(pwd())
            fermi = vaspio_fermi(pwd())
            break

        @default
            sorry()
            break
    end

    # dump the Kohn-Sham data to files with IR format
    irio_lattice(pwd(), latt)
    irio_kmesh(pwd(), kmesh, weight)
    irio_tetra(pwd(), volt, itet)
    irio_eigen(pwd(), enk, occupy)
    irio_projs(pwd(), chipsi)
    irio_fermi(pwd(), fermi)

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
