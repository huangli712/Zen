#
# project : pansy
# source  : vasp.jl
# author  : Li Huang (lihuang.dmft@gmail.com)
# status  : unstable
# comment :
#
# last modified: 2020/12/16
#

"""
    vasp_init(it::IterInfo)

Check the runtime environment of vasp, prepare necessary input files
"""
function vasp_init(it::IterInfo)
    # prepare essential input files
    if it.dft_dmft_iter == 0
        # copy POTCAR and POSCAR
        cp("../POTCAR", joinpath(pwd(), "POTCAR"))
        cp("../POSCAR", joinpath(pwd(), "POSCAR"))

        # generate INCAR automatically
        vasp_incar()
    end

    # check essential input files
    if it.dft_dmft_iter >= 1
        flist = ("INCAR", "POSCAR", "POTCAR")
        for i in eachindex(flist)
            filename = flist[i]
            if !isfile(filename)
                error("Please make sure the file $filename is available")
            end
        end
    end
end

"""
    vasp_run(it::IterInfo)

Execute the vasp program
"""
function vasp_run(it::IterInfo)
    # get the home directory of vasp
    dft_home = query_dft()

    # determine mpi prefix (whether the vasp is executed sequentially)
    mpi_prefix = parse_toml("MPI.toml", "dft", false)

    # select suitable vasp program
    if _d("lspinorb")
        vasp_exec = "$dft_home/vasp_ncl"
    else
        vasp_exec = "$dft_home/vasp_std"
    end

    # assemble command
    if isnothing(mpi_prefix)
        vasp_cmd = vasp_exec
    else
        vasp_cmd = split("$mpi_prefix $vasp_exec", " ")
    end

    # invoke it, the terminal output is redirected to vasp.out
    run(pipeline(`$vasp_cmd`, stdout = "vasp.out"))
end

"""
    vasp_save(it::IterInfo)

Backup the output files of vasp if necessary
"""
function vasp_save(it::IterInfo)
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
end

"""
    vasp_incar()

Generate an INCAR file. it will be used only when the dft engine is vasp
"""
function vasp_incar()
    # open the iostream
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
    @cswitch smear begin
        @case "m-p"
            write(ios, "ISMEAR   = 2 \n")
            break

        @case "gauss"
            write(ios, "ISMEAR   = 0 \n")
            break

        @case "tetra"
            write(ios, "ISMEAR   =-5 \n")
            break

        @default
            write(ios, "ISMEAR   = 2 \n")
            break
    end

    # for k-mesh density
    kmesh = _d("kmesh")
    @cswitch kmesh begin
        @case "accurate"
            write(ios, "KSPACING = 0.1 \n")
            break

        @case "medium"
            write(ios, "KSPACING = 0.2 \n")
            break

        @case "coarse"
            write(ios, "KSPACING = 0.4 \n")
            break

        @default # very coarse 
            write(ios, "KSPACING = 0.5 \n")
            break
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
    # if spin-orbit coupling is on, spin orientations must be polarized
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

    # close the iostream
    close(ios)
end

"""
    vaspio_lattice(f::AbstractString)

Reading vasp's POSCAR file, return crystallography information. Here `f`
means only the directory that contains POSCAR 
"""
function vaspio_lattice(f::AbstractString)
    # open the iostream
    fin = open(joinpath(f, "POSCAR"), "r")

    # get the case
    case = strip(readline(fin))

    # get the scaling factor
    scale = parse(F64, readline(fin))

    # get the basis vector
    bvec = zeros(F64, 3, 3)
    bvec[1, :] = parse.(F64, line_to_array(fin))
    bvec[2, :] = parse.(F64, line_to_array(fin))
    bvec[3, :] = parse.(F64, line_to_array(fin))

    # get the symbol list
    symbols = line_to_array(fin)

    # get the number of sorts of atoms
    nsorts = length(symbols)

    # get the number list 
    numbers = parse.(I64, line_to_array(fin))

    # get the total number of atoms
    natoms = sum(numbers)

    # create atom list, whose value is related to the sorts of atoms
    atom_list = zeros(I64, natoms)
    curr_index = 0
    k = 0
    for i = 1:length(numbers)
        curr_index = curr_index + 1
        for j = 1:numbers[i]
            k = k + 1
            atom_list[k] = curr_index
        end
    end

    # get the coordinates of atoms
    posi_list = zeros(F64, natoms, 3)
    readline(fin)
    for i = 1:natoms
        posi_list[i, :] = parse.(F64, line_to_array(fin)[1:3])
    end

    # close the iostream
    close(fin)

    # return the desired arrays
    return nsorts, natoms, symbols, atom_list, posi_list
end

"""
    vaspio_kmesh(f::AbstractString)

Reading vasp's IBZKPT file, return kmesh and weight. Here `f` means
only the directory that contains IBZKPT
"""
function vaspio_kmesh(f::AbstractString)
    # open the iostream
    fin = open(joinpath(f, "IBZKPT"), "r")

    # extract number of k-points
    readline(fin)
    nkpt = parse(I64, readline(fin))
    readline(fin)

    # create arrays 
    kmesh = zeros(F64, nkpt, 3)
    weight = zeros(F64, nkpt)

    # read in the k-points and their weights
    for i = 1:nkpt
        arr = parse.(F64, line_to_array(fin))
        kmesh[i, 1:3] = arr[1:3]
        weight[i] = arr[4]
    end

    # close the iostream
    close(fin)

    # return the desired arrays
    return kmesh, weight
end

"""
    vaspio_tetra(f::AbstractString)

Reading vasp's IBZKPT file, return tetrahedra information. Here `f` means
only the directory that contains IBZKPT
"""
function vaspio_tetra(f::AbstractString)
    # open the iostream
    fin = open(joinpath(f, "IBZKPT"), "r")

    # extract number of k-points
    readline(fin)
    nkpt = parse(I64, readline(fin))
    readline(fin)

    # read in the k-points and their weights
    # skip nkpt lines
    for i = 1:nkpt
        readline(fin)
    end

    # read in the tetrahedra information
    # skip one empty line
    readline(fin)

    # extract total number of tetrahedra and volume of a tetrahedron 
    arr = line_to_array(fin)
    ntet = parse(I64, arr[1])
    volt = parse(F64, arr[2])

    # create arrays
    itet = zeros(I64, ntet, 5)

    # parse the input tetrahedra information
    for t = 1:ntet
        itet[t, :] = parse.(I64, line_to_array(fin))
    end

    # close the iostream
    close(fin)

    # return the desired arrays
    return volt, itet
end

"""
    vaspio_eigen(f::AbstractString)

Reading vasp's EIGENVAL file, return energy band information. Here `f`
means only the directory that contains EIGENVAL
"""
function vaspio_eigen(f::AbstractString)
    # open the iostream
    fin = open(joinpath(f, "EIGENVAL"), "r")

    # determine number of spins
    nspin = parse(I64, line_to_array(fin)[end])

    # skip for lines
    for i = 1:4
        readline(fin)
    end

    # read in some key parameters: nelect, nkpt, nbands 
    nelect, nkpt, nband = parse.(I64, line_to_array(fin))

    # create arrays
    enk = zeros(F64, nkpt, nband, nspin)
    occupy = zeros(F64, nkpt, nband, nspin)

    # read in the energy bands and the corresponding occupations
    for i = 1:nkpt
        readline(fin)
        readline(fin)
        for j = 1:nband
            arr = line_to_array(fin)
            # for spin unpolarized case
            if nspin === 1
                enk[i, j, 1] = parse(F64, arr[2])
                occupy[i, j, 1] = parse(F64, arr[3])
                # for spin polarized case
            else
                enk[i, j, 1] = parse(F64, arr[2])
                enk[i, j, 2] = parse(F64, arr[3])
                occupy[i, j, 1] = parse(F64, arr[4])
                occupy[i, j, 2] = parse(F64, arr[5])
            end
        end
    end

    # close the iostream
    close(fin)

    # return the desired arrays
    return enk, occupy
end

"""
    vaspio_projcar(f::AbstractString)

Reading vasp's PROJCAR file, return raw projector matrix. Here `f` means
only the directory that contains PROJCAR
"""
function vaspio_projcar(f::AbstractString)
    # get key parameters from the LOCPROJ file
    nspin, nkpt, nband, nproj, nsite, sites, projs, groups = vaspio_locproj(f, true)

    # open the iostream
    fin = open(f * "/PROJCAR", "r")

    # create arrays
    chipsi = zeros(C64, nproj, nband, nkpt, nspin)

    # read in raw projector data
    for site = 1:nsite
        # extract site information
        _site = parse(I64, line_to_array(fin)[2])
        @assert _site === site

        # skip one empty line
        readline(fin)

        for spin = 1:nspin
            for kpt = 1:nkpt
                # extract k-point and spin information
                arr = line_to_array(fin)
                _kpt = parse(I64, arr[2])
                _spin = parse(I64, arr[4])
                @assert _kpt === kpt
                @assert _spin === spin

                # skip two empty lines
                readline(fin)
                readline(fin)

                # parse the input data
                for band = 1:nband
                    arr = parse.(F64, line_to_array(fin))
                    for proj = 1:groups[site]
                        cmplx = arr[2*proj] + arr[2*proj+1]im
                        chipsi[proj, band, kpt, spin] = cmplx
                    end
                end

                # skip one empty line
                readline(fin)
            end
        end
    end

    # close the iostream
    close(fin)

    # return the desired arrays
    return chipsi
end

"""
    vaspio_locproj(f::AbstractString, read_param_only::Bool)

Reading vasp's LOCPROJ file, return raw projector matrix. Here `f` means
only the directory that contains LOCPROJ
"""
function vaspio_locproj(f::AbstractString, read_param_only::Bool = false)
    # open the iostream
    fin = open(f * "/LOCPROJ", "r")

    # extract number of spins (nspin), number of k-points (nkpt),
    # number of bands (nband), and number of projectors (nproj)
    nspin, nkpt, nband, nproj = parse.(I64, line_to_array(fin)[1:4])

    # find out how many sites are there (nsite)
    # projs contains the specifications of the projectors
    sites = zeros(I64, nproj)
    projs = Array{String}(undef, nproj)
    for i = 1:nproj
        arr = line_to_array(fin)
        sites[i] = parse(I64, arr[2])
        # get rid of the ":" char
        projs[i] = replace(arr[end], ":" => "")
    end
    usites = union(sites)
    nsite = length(usites)

    # find out how many projectors are there for a given site
    groups = zeros(I64, nsite)
    for site = 1:nsite
        groups[site] = length(findall(x -> x === usites[site], sites))
    end

    # additional check, make sure nproj is equal to the sum of groups
    @assert nproj === sum(groups)

    if read_param_only
        # return only parameters
        return nspin, nkpt, nband, nproj, nsite, sites, projs, groups
    else
        # create arrays
        chipsi = zeros(C64, nproj, nband, nkpt, nspin)

        # read in raw projector data
        readline(fin)
        for spin = 1:nspin
            for kpt = 1:nkpt
                for band = 1:nband
                    # extract some indices information
                    arr = line_to_array(fin)
                    _spin = parse(I64, arr[2])
                    _kpt = parse(I64, arr[3])
                    _band = parse(I64, arr[4])

                    # check consistency of parameters
                    @assert _spin === spin
                    @assert _kpt === kpt
                    @assert _band === band

                    # parse the input data
                    _proj = 0
                    for site = 1:nsite
                        for proj = 1:groups[site]
                            _proj = _proj + 1
                            _re, _im = parse.(F64, line_to_array(fin)[2:3])
                            chipsi[_proj, band, kpt, spin] = _re + _im * im
                        end
                    end

                    # skip one empty line
                    readline(fin)
                end
            end
        end
    end

    # close the iostream
    close(fin)

    # return the desired arrays
    return chipsi
end

"""
    vaspio_chgcar(f::AbstractString)
"""
function vaspio_chgcar(f::AbstractString) end
