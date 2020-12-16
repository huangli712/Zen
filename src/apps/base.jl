#
# project : pansy
# source  : base.jl
# author  : Li Huang (lihuang.dmft@gmail.com)
# status  : unstable
# comment :
#
# last modified: 2020/12/16
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
    dft_init(it::IterInfo)

To examine whether the dft runtime environment is ready
"""
function dft_init(it::IterInfo)
    # enter dft directory
    cd("dft")

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
    dft_home = query_dft()
    mpi_prefix = parse_toml("MPI.toml", "dft", false)

    cd("dft")

    if _d("engine") === "vasp"
        if _d("lspinorb")
            vasp_exec = "$dft_home/vasp_ncl"
        else
            vasp_exec = "$dft_home/vasp_std"
        end

        if isnothing(mpi_prefix)
            vasp_cmd = vasp_exec
        else
            vasp_cmd = split("$mpi_prefix $vasp_exec", " ")
        end
        run(pipeline(`$vasp_cmd`, stdout = "vasp.out"))
    else
        sorry()
    end

    cd("..")
end

"""
    dft_save(it::IterInfo)

Backup the essential dft calculated results for next iterations
"""
function dft_save(it::IterInfo)
    cd("dft")

    if _d("engine") === "vasp"
        if it.dft_dmft_iter == 0
            cp("INCAR", "INCAR.$(it.dft_dmft_iter)")
            cp("CHGCAR", "CHGCAR.$(it.dft_dmft_iter)")
            cp("OUTCAR", "OUTCAR.$(it.dft_dmft_iter)")
            cp("PROJCAR", "PROJCAR.$(it.dft_dmft_iter)")
            cp("LOCPROJ", "LOCPROJ.$(it.dft_dmft_iter)")
            cp("EIGENVAL", "EIGENVAL.$(it.dft_dmft_iter)")
            cp("vasp.out", "vasp.out.$(it.dft_dmft_iter)")
            cp("vasprun.xml", "vasprun.xml.$(it.dft_dmft_iter)")
        else
            sorry()
        end
    else
        sorry()
    end

    cd("..")
end
