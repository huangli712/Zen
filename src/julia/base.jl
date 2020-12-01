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
    for i = 1:Param(PIMP, "nsite")
        if isdir("impurity.$i")
            rm("impurity.$i", force = true, recursive = true)
        end
        mkdir("impurity.$i")
    end
end

"""
    make_incar(case::String, d::Dict{String,Any})

Generate an INCAR for vasp, which is suitable for an initial self-consistent run
"""
function make_incar(case::String, d::Dict{String,Any})
    ios = open("INCAR", "w")

    # standard part
    write(ios,"System   = $case \n")
    write(ios,"PREF     = Accurate \n")
    write(ios,"EDIFF    = 1E-8 \n")
    write(ios,"ALGO     = Normal \n")
    write(ios,"LASPH    = .TRUE. \n")

    # customize your INCAR according to the case.toml
    #
    # for smearing
    # ismear ==  2: m-p scheme -> metal
    # ismear ==  1: m-p scheme -> metal (default)
    # ismear ==  0: gauss scheme -> metal, semiconductor or insulator
    # ismear == -5: tetrahedron scheme with correction -> insulator
    if haskey(d, "smear")
        if d["smear"] == "m-p"
            write(ios,"ISMEAR   = 2 \n")
        elseif d["smear"] == "gauss"
            write(ios,"ISMEAR   = 0 \n")
        elseif d["smear"] == "tetra"
            write(ios,"ISMEAR   = -5 \n")
        else
            write(ios,"ISMEAR   = 1 \n")
        end
    else
        write(ios,"ISMEAR   = 1 \n")
    end

    # for k-mesh density
    if haskey(d, "kgrid")
        if d["kgrid"] == "accurate"
            write(ios,"KSPACING = 0.1 \n")
        elseif d["kgrid"] == "medium"
            write(ios,"KSPACING = 0.2 \n")
        elseif d["kgrid"] == "coarse"
            write(ios,"KSPACING = 0.4 \n")
        else # very coarse 
            write(ios,"KSPACING = 0.5 \n")
        end
    else
        write(ios,"KSPACING = 0.5 \n")
    end

    # for magnetic moment
    if haskey(d, "magmom")
        str = d["magmom"]
        write(ios,"MAGMOM   = $str \n")
    end

    # for symmetry
    # isym == 2: turn on  symmetry
    # isym == 0: turn off symmetry
    if haskey(d, "lsymm")
        if d["lsymm"]
            write(ios,"ISYM     = 2 \n")
        else
            write(ios,"ISYM     = 0 \n")
        end
    else
        write(ios,"ISYM     = 2 \n")
    end

    # for spin polarizations
    # if spin-orbit coupling is on, then spins must be polarized
    if haskey(d, "lspins")
        if d["lspins"]
            write(ios,"ISPIN    = 2 \n")
        else
            write(ios,"ISPIN    = 1 \n")
        end
    else
        write(ios,"ISPIN    = 1 \n")
    end

    # for spin-orbit coupling
    if haskey(d, "lspinorb")
        if d["lspinorb"]
            write(ios,"LSORBIT  = .TRUE. \n")
        else
            write(ios,"LSORBIT  = .FALSE. \n")
        end
    else
        write(ios,"LSORBIT  = .FALSE. \n")
    end

    # for optimized projectors
    if haskey(d, "lopt") && haskey(d, "window")
        if d["lopt"]
            write(ios,"LORBIT   = 14 \n")
            emin = d["window"][1]
            write(ios,"EMIN     = $emin \n")
            emax = d["window"][2]
            write(ios,"EMAX     = $emax \n")
        end
    end

    # for local orbitals and projectors
    if haskey(d, "lproj") && haskey(d, "nproj") && haskey(d, "sproj")
        if d["lproj"]
            for p in 1:d["nproj"]
                str = d["sproj"][p]
                write(ios, "LOCPROJ  = $str \n")
            end
        end
    end

    close(ios)
end

"""
    dft_init(it::IterInfo)

To examine whether the dft runtime environment is ready
"""
function dft_init(it::IterInfo)
    # enter dft directory
    cd("dft")

    # copy essential input files
    if it.dft_dmft_iter == 0
        if Param(PDFT, "engine") === "vasp"
           cp("../POTCAR", pwd() * "/POTCAR")
           cp("../POSCAR", pwd() * "/POSCAR")
        else
           sorry()
        end
    end

    # generate essential input files
    if it.dft_dmft_iter == 0
        if Param(PDFT, "engine") === "vasp"
            make_incar()
        else
            sorry()
        end
    end

    # check essential input files
    if it.dft_dmft_iter >= 1
        if Param(PDFT, "engine") === "vasp"
            if !isfile("INCAR") || !isfile("POSCAR") || !isfile("POTCAR")
                error("Please make sure the existence of following files: INCAR, POSCAR, and POTCAR")
            end
        else
           sorry()
        end
    end

    # enter the parent directory
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
