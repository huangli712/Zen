"""
    make_trees()

Build the working directories at advance 
"""
function make_trees(d::Dict{String,Any})
    mkdir("dft")
    mkdir("dmft1")
    mkdir("dmft2")
    for i = 1:d["nimp"]
        mkdir("impurity.$i")
    end
end

function dft_init(it::IterInfo, d::Dict{String,Any})
    cd("dft")

    # copy essential input files
    if it.dft_dmft_iter == 0
        if d["engine"] == "vasp"
           cp("../POTCAR", pwd() * "/POTCAR")
           cp("../POSCAR", pwd() * "/POSCAR")
        else
           sorry()
        end
    end

    # generate essential input files
    if it.dft_dmft_iter == 0
        if d["engine"] == "vasp"
            make_incar()
            make_kpoints()
        else
            sorry()
        end
    end

    # check essential input files
    if it.dft_dmft_iter >= 1
        if d["engine"] == "vasp"
            if !isfile("INCAR") || !isfile("KPOINTS") || !isfile("POSCAR") || !isfile("POTCAR")
                error("Please make sure the existence of following files: INCAR, KPOINTS, POSCAR, and POTCAR")
            end
        else
           sorry()
        end
    end

    cd("..")
end

function dft_run()
end

function dft_save()
end

"""
    dft_driver(IterInfo, Dict{String,Any})

Drive the engine of density functional theory to carry out calculations
"""
function dft_driver(it::IterInfo, d::Dict{String,Any})
    cd("dft")

    if it.dft_iter == 0
        cp("../INCAR", pwd() * "/INCAR")
        cp("../POSCAR", pwd() * "/POSCAR")
        cp("../POTCAR", pwd() * "/POTCAR")
        cp("../KPOINTS", pwd() * "/KPOINTS")
    else
    end

    mpi_prefix = "mpiexec"
    mpi_options = "-n" 
    mpi_ncpu = "4"
    vasp_exec = "/home/soft/vasp/vasp.6.1.1/vasp.6.1.1/bin/vasp_std"

    run(`$mpi_prefix $mpi_options $mpi_ncpu $vasp_exec`)

    it.dft_iter += 1

    cd("..")
end

function dmft_driver(d::Dict{String,Any})
end
