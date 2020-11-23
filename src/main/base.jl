"""
    make_trees(d:Dict{String,Any})

Build the working directories at advance 
"""
function make_trees(d::Dict{String,Any})
    if isdir("dft")
        rm("dft", force = true, recursive = true)
    end
    mkdir("dft")

    if isdir("dmft1")
        rm("dmft1", force = true, recursive = true)
    end
    mkdir("dmft1")

    if isdir("dmft2")
        rm("dmft2", force = true, recursive = true)
    end
    mkdir("dmft2")

    for i = 1:d["nimp"]
        if isdir("impurity.$i")
            rm("impurity.$i", force = true, recursive = true)
        end
        mkdir("impurity.$i")
    end
end

"""
    make_incar(case::String, d::Dict{String,Any})

Generate an INCAR for vasp, which is only suitable for an initial self-consistent run
"""
function make_incar(case::String, d::Dict{String,Any})
    ios = open("INCAR", "w")

    write(ios,"System   = $case \n")
    write(ios,"PREF     = Accurate \n")
    write(ios,"EDIFF    = 1E-8 \n")
    write(ios,"ALGO     = Normal \n")
    write(ios,"LASPH    = .TRUE. \n")

    # customize your INCAR according to the case.toml
    write(ios,"ISMEAR   = 0 \n")

    if d["kgrid"] == "accurate"
        write(ios,"KSPACING = 0.1 \n")
    elseif d["kgrid"] == "medium"
        write(ios,"KSPACING = 0.2 \n")
    elseif d["kgrid"] == "coarse"
        write(ios,"KSPACING = 0.4 \n")
    else 
        write(ios,"KSPACING = 0.5 \n")
    end

    if d["lspins"]
        write(ios,"ISPIN    = 2 \n")
    end

    if d["lspinorb"] 
        write(ios,"LSORBIT  = .TRUE. \n")
    end

    if d["projector"]["lproj"]
        for p in 1:d["projector"]["nproj"]
            str = d["projector"]["sproj"][p]
            write(ios, "LOCPROJ  = $str\n")
        end
    end

    close(ios)
end

"""
    dft_init(it::IterInfo, case::String, d::Dict{String,Any})

To examine whether the dft runtime environment is ready
"""
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

"""
    dft_run(it::IterInfo, d::Dict{String,Any})

Execute the desired dft engine parallelly or sequentially
"""
function dft_run(it::IterInfo, d::Dict{String,Any})
    dft_home = query_dft(d)
    mpi_prefix = parse_toml("MPI.toml", "dft", false) 

    cd("dft")

    if d["engine"] == "vasp"
        if d["lspinorb"]
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
    dft_save(it::IterInfo, d::Dict{String,Any})

Backup the essential dft calculated results for next iterations
"""
function dft_save(it::IterInfo, d::Dict{String,Any})
    cd("dft")

    if d["engine"] == "vasp"
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
