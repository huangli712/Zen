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

    mpi_prefix = "mpiexec -n 4"
    vasp_exec = "/home/soft/vasp/vasp.6.1.1/vasp.6.1.1/bin/vasp_std"
    run(`$vasp_exec`)

    it.dft_iter += 1

    cd("..")
end

function dmft_driver(d::Dict{String,Any})
end
