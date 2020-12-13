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
    make_incar()

Generate an INCAR for vasp, which is suitable for an initial self-consistent run
"""
function make_incar()
    ios = open("INCAR", "w")

    # standard part
    case = _c("case")
    write(ios, "System   = $case \n")
    write(ios, "PREF     = Accurate \n")
    write(ios, "EDIFF    = 1E-8 \n")
    write(ios, "ALGO     = Normal \n")
    write(ios, "LASPH    = .TRUE. \n")

    # customize your INCAR according to the case.toml
    #
    # for smearing
    smear = _d("smear")
    if smear === "m-p"
        write(ios, "ISMEAR   = 2 \n")
    elseif smear === "gauss"
        write(ios, "ISMEAR   = 0 \n")
    elseif smear === "tetra"
        write(ios, "ISMEAR   = -5 \n")
    else
        write(ios, "ISMEAR   = 1 \n")
    end

    # for k-mesh density
    kmesh = _d("kmesh")
    if kmesh === "accurate"
        write(ios, "KSPACING = 0.1 \n")
    elseif kmesh === "medium"
        write(ios, "KSPACING = 0.2 \n")
    elseif kmesh === "coarse"
        write(ios, "KSPACING = 0.4 \n")
    else # very coarse 
        write(ios, "KSPACING = 0.5 \n")
    end

    # for magnetic moment
    magmom = _d("magmom")
    if !isa(magmom, Missing)
        write(ios, "MAGMOM   = $magmom \n")
    end

    # for symmetry
    lsymm = _d("lsymm")
    if lsymm
        write(ios, "ISYM     = 2 \n")
    else
        write(ios, "ISYM     = 0 \n")
    end

    # for spin polarizations
    # if spin-orbit coupling is on, then spins must be polarized
    lspins = _d("lspins")
    if lspins
        write(ios, "ISPIN    = 2 \n")
    else
        write(ios, "ISPIN    = 1 \n")
    end

    # for spin-orbit coupling
    lspinorb = _d("lspinorb")
    if lspinorb
        write(ios, "LSORBIT  = .TRUE. \n")
    else
        write(ios, "LSORBIT  = .FALSE. \n")
    end

    # for optimized projectors
    window = _d("window")
    loptim = _d("loptim")
    if !isa(window, Missing) && !isa(loptim, Missing)
        if loptim
            write(ios, "LORBIT   = 14 \n")
            emin = window[1]
            write(ios, "EMIN     = $emin \n")
            emax = window[2]
            write(ios, "EMAX     = $emax \n")
        end
    end

    # for local orbitals and projectors
    lproj = _d("lproj")
    nproj = _d("nproj")
    sproj = _d("sproj")
    if !isa(lproj, Missing) && !isa(nproj, Missing) && !isa(sproj, Missing)
        if lproj
            for p = 1:nproj
                str = sproj[p]
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
        if _d("engine") === "vasp"
            cp("../POTCAR", pwd() * "/POTCAR")
            cp("../POSCAR", pwd() * "/POSCAR")
        else
            sorry()
        end
    end

    # generate essential input files
    if it.dft_dmft_iter == 0
        if _d("engine") === "vasp"
            make_incar()
        else
            sorry()
        end
    end

    # check essential input files
    if it.dft_dmft_iter >= 1
        if _d("engine") === "vasp"
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
