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

    # choose suitable driver function according to dft engine
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

    # choose suitable driver function according to dft engine
    engine = _d("engine")
    @cswitch engine begin
        @case "vasp"
            println("< Adaptor: Get Kohn-Sham Data >")

            # read in lattice structure
            println("  Lattice")
            latt = vaspio_lattice(pwd())

            # read in kmesh and the corresponding weights
            println("  Kmesh")
            println("  Weight")
            kmesh, weight = vaspio_kmesh(pwd())

            # read in tetrahedron data if they are available
            if _d("smear") === "tetra"
                println("  Tetrahedron")
                volt, itet = vaspio_tetra(pwd())
            end

            # read in band structure and the corresponding occupancies
            println("  Enk")
            println("  Occupation")
            enk, occupy = vaspio_eigen(pwd())

            # read in raw projectors, traits, and groups
            println("  Projector (Trait and Group)")
            PT, PG, chipsi = vaspio_projs(pwd())

            # read in fermi level
            println("  Fermi Level\n")
            fermi = vaspio_fermi(pwd())
            break

        @default
            sorry()
            break
    end

    # well, now we have the Kohn-Sham data. but they can not be written
    # directly. we have to process them carefully
    println("< Adaptor: Eat Kohn-Sham Data >")
    println("  Grouping")
    plo_group(PG)
    println("  Rotating")
    println("  Leveling")
    println("  Filtering")
    println("  Orthogonalizing\n")

    # dump the Kohn-Sham data to files with IR format
    println("< Adaptor: Put Kohn-Sham Data >")

    # write lattice structure
    println("  Lattice")
    irio_lattice(pwd(), latt)

    # write kmesh and the corresponding weights
    println("  Kmesh")
    println("  Weight")
    irio_kmesh(pwd(), kmesh, weight)

    # write tetrahedron data if they are available
    if _d("smear") === "tetra"
        println("  Tetrahedron")
        irio_tetra(pwd(), volt, itet)
    end

    # write band structure and the corresponding occupancies
    println("  Enk")
    println("  Occupation")
    irio_eigen(pwd(), enk, occupy)

    # write projectors, traits, and groups
    println("  Projector (Trait and Group)")
    irio_projs(pwd(), chipsi)

    # write fermi level
    println("  Fermi Level\n")
    irio_fermi(pwd(), fermi)

    println("< Adaptor: View Overlap Matrix >")
    println("< Adaptor: View Density Matrix >")

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
