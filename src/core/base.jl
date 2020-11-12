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

function make_incar(case::String, d::Dict{String,Any})
    ios = open("INCAR", "w")

    write(ios,"System   = $case \n")
    write(ios,"PREF     = Accurate \n")
    write(ios,"EDIFF    = 1E-8 \n")
    write(ios,"ALGO     = Normal \n")
    write(ios,"ISMEAR   = 0 \n")
    write(ios,"LASPH    = .TRUE. \n")
    write(ios,"KSPACING = 0.5 \n")

    if d["lspins"]
        write(ios,"ISPIN    = .TRUE. \n")
    end

    if d["lspinorb"] 
        write(ios,"LSORBIT  = .TRUE. \n")
    end

    if d["projector"]["lproj"]
        for p in 1:d["projector"]["nproj"]
            str = d["projector"]["sproj"][p]
            write(ios, "LOCPROJ  = $str")
        end
    end

    close(ios)
end

function dft_init(it::IterInfo, case::String, d::Dict{String,Any})
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
            make_incar(case, d)
        else
            sorry()
        end
    end

    # check essential input files
    if it.dft_dmft_iter >= 1
        if d["engine"] == "vasp"
            if !isfile("INCAR") || !isfile("POSCAR") || !isfile("POTCAR")
                error("Please make sure the existence of following files: INCAR, POSCAR, and POTCAR")
            end
        else
           sorry()
        end
    end

    cd("..")
end

function dft_run(it::IterInfo, d::Dict{String,Any})
    dft_home = check_dft(d)
    mpi_prefix = parse_mpi("dft") 

    cd("dft")

    if d["engine"] == "vasp"
        vasp_exec = ""
        if d["lspinorb"]
            global vasp_exec = "$dft_home/vasp_ncl"
        else
            global vasp_exec = "$dft_home/vasp_std"
        end
        vasp_cmd = split("$mpi_prefix $vasp_exec", " ")
        @show `$vasp_cmd`
        exit(-1)
        run(pipeline(`$vasp_cmd`, stdout = "vasp.out"))
    else
        sorry()
    end

    cd("..")
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
