#
# project : pansy
# source  : vasp.jl
# author  : Li Huang (lihuang.dmft@gmail.com)
# status  : unstable
# comment :
#
# last modified: 2020/12/31
#

"""
    vasp_init(it::IterInfo)

Check the runtime environment of vasp, prepare necessary input files
"""
function vasp_init(it::IterInfo)
    # prepare essential input files
    if it.dmft_cycle == 0
        # copy POTCAR and POSCAR
        cp("../POTCAR", joinpath(pwd(), "POTCAR"), force = true)
        cp("../POSCAR", joinpath(pwd(), "POSCAR"), force = true)

        # generate INCAR automatically
        vasp_incar(it._dft_fermi)
    end

    # check essential input files
    if it.dmft_cycle >= 1
        flist = ("INCAR", "POSCAR", "POTCAR")
        for i in eachindex(flist)
            filename = flist[i]
            if !isfile(filename)
                error("Please make sure the file $filename is available")
            end
        end

        # maybe we need to update INCAR file here
        vasp_incar(it.dmft_fermi)
    end

    # well, maybe we need to generate the KPOINTS file by ourselves
    if __d("kmesh")
        vasp_kpoints()
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
    mpi_prefix = parse_toml("../MPI.toml", "dft", false)

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

    # launch it, the terminal output is redirected to vasp.out
    run(pipeline(`$vasp_cmd`, stdout = "vasp.out"))
end

"""
    vasp_save(it::IterInfo)

Backup the output files of vasp if necessary. Furthermore, the fermi level
in IterInfo struct is also updated (IterInfo._dft_fermi)
"""
function vasp_save(it::IterInfo)
    # store the data files
    if it.dmft_cycle == 0
        cp("INCAR", "INCAR.$(it.dmft_cycle)", force = true)
        cp("CHGCAR", "CHGCAR.$(it.dmft_cycle)", force = true)
        cp("OUTCAR", "OUTCAR.$(it.dmft_cycle)", force = true)
        cp("PROJCAR", "PROJCAR.$(it.dmft_cycle)", force = true)
        cp("LOCPROJ", "LOCPROJ.$(it.dmft_cycle)", force = true)
        cp("EIGENVAL", "EIGENVAL.$(it.dmft_cycle)", force = true)
        cp("vasp.out", "vasp.out.$(it.dmft_cycle)", force = true)
        cp("vasprun.xml", "vasprun.xml.$(it.dmft_cycle)", force = true)
    else
        sorry()
    end

    # anyway, the fermi level is extracted from DOSCAR, and its value will
    # be saved at IterInfo._dft_fermi
    it._dft_fermi = vaspio_fermi(pwd())
end

"""
    vasp_incar(fermi::F64)

Generate an INCAR file. it will be used only when the dft engine is vasp
"""
function vasp_incar(fermi::F64)
    # open the iostream
    ios = open("INCAR", "w")

    # standard part
    case = _c("case")
    write(ios, "System   = $case \n")
    write(ios, "PREF     = Accurate \n")
    write(ios, "EDIFF    = 1E-8 \n")
    write(ios, "ALGO     = Normal \n")
    write(ios, "LASPH    = .TRUE. \n")
    write(ios, "LMAXMIX  = 6 \n")

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

    # for kmesh density
    # if kmesh == "file", then vasp_kpoints() will be used to generate
    # the KPOINTS file
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

        @case "file"
            break

        @default # very coarse kmesh
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
    else # ignore the symmetry completely
        write(ios, "ISYM     =-1 \n")
    end

    # for spin polarizations
    # if spin-orbit coupling is on, spin orientations must be polarized
    lspins = _d("lspins") || _d("lspinorb")
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
    ewidth = _d("ewidth")
    loptim = _d("loptim")
    if !isa(ewidth, Missing) && !isa(loptim, Missing)
        if loptim
            write(ios, "LORBIT   = 14 \n")
            emin = fermi - ewidth
            write(ios, "EMIN     = $emin \n")
            emax = fermi + ewidth
            write(ios, "EMAX     = $emax \n")
        end
    end

    # for local orbitals and projectors
    lproj = _d("lproj")
    sproj = _d("sproj")
    if !isa(lproj, Missing) && !isa(sproj, Missing)
        if lproj
            for p in eachindex(sproj)
                str = sproj[p]
                write(ios, "LOCPROJ  = $str \n")
            end
        end
    end

    # close the iostream
    close(ios)
end

"""
    vasp_kpoints()

Generate a KPOINTS file for vasp
"""
function vasp_kpoints()
    # if the KPOINTS file is available, we do nothing, or else we will
    # create a new one
    if !isfile("KPOINTS")
        # open the iostream
        ios = open("KPOINTS", "w")

        # write the body
        write(ios, "Automatic mesh")
        write(ios, "0")
        write(ios, "Monkhorst-Pack")
        write(ios, "$i $i $i")
        write(ios, " 0  0  0")

        # close the iostream
        close(ios)
    end
end

"""
    vasp_files(f::String)

Check the essential output files by vasp. Here `f` means only the directory
that contains the desired files
"""
function vasp_files(f::String)
    @assert isfile(joinpath(f, "POSCAR")) &&
            isfile(joinpath(f, "IBZKPT")) &&
            isfile(joinpath(f, "DOSCAR")) &&
            isfile(joinpath(f, "CHGCAR")) &&
            isfile(joinpath(f, "LOCPROJ")) &&
            isfile(joinpath(f, "EIGENVAL"))
end

"""
    vaspio_lattice(f::String)

Reading vasp's POSCAR file, return crystallography information. Here `f`
means only the directory that contains POSCAR
"""
function vaspio_lattice(f::String)
    # open the iostream
    fin = open(joinpath(f, "POSCAR"), "r")

    # get the case
    _case = string(strip(readline(fin)))

    # get the scaling factor
    scale = parse(F64, readline(fin))

    # get the lattice vectors
    lvect = zeros(F64, 3, 3)
    lvect[1, :] = parse.(F64, line_to_array(fin))
    lvect[2, :] = parse.(F64, line_to_array(fin))
    lvect[3, :] = parse.(F64, line_to_array(fin))

    # get the symbol list
    symbols = line_to_array(fin)

    # get the number of sorts of atoms
    nsort = length(symbols)

    # get the number list
    numbers = parse.(I64, line_to_array(fin))

    # get the total number of atoms
    natom = sum(numbers)

    # now all the parameters are ready
    # we would like to create Lattice struct here
    latt = Lattice(_case, scale, nsort, natom)

    # update latt using the available data
    latt.lvect = lvect
    for i = 1:nsort
        latt.sorts[i, 1] = string(symbols[i])
        latt.sorts[i, 2] = numbers[i]
    end

    # get the atom list
    k = 0
    for i = 1:nsort
        for j = 1:numbers[i]
            k = k + 1
            latt.atoms[k] = symbols[i]
        end
    end
    # sanity check
    @assert k === natom

    # get the coordinates of atoms
    readline(fin)
    for i = 1:natom
        latt.coord[i, :] = parse.(F64, line_to_array(fin)[1:3])
    end

    # close the iostream
    close(fin)

    # return the desired struct
    return latt
end

"""
    vaspio_kmesh(f::String)

Reading vasp's IBZKPT file, return kmesh and weight. Here `f` means
only the directory that contains IBZKPT
"""
function vaspio_kmesh(f::String)
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
    vaspio_tetra(f::String)

Reading vasp's IBZKPT file, return tetrahedra information. Here `f` means
only the directory that contains IBZKPT
"""
function vaspio_tetra(f::String)
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
    vaspio_eigen(f::String)

Reading vasp's EIGENVAL file, return energy band information. Here `f`
means only the directory that contains EIGENVAL
"""
function vaspio_eigen(f::String)
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
    enk = zeros(F64, nband, nkpt, nspin)
    occupy = zeros(F64, nband, nkpt, nspin)

    # read in the energy bands and the corresponding occupations
    for i = 1:nkpt
        readline(fin)
        readline(fin)
        for j = 1:nband
            arr = line_to_array(fin)
            if nspin === 1 # for spin unpolarized case
                enk[j, i, 1] = parse(F64, arr[2])
                occupy[j, i, 1] = parse(F64, arr[3])
            else # for spin polarized case
                enk[j, i, 1] = parse(F64, arr[2])
                enk[j, i, 2] = parse(F64, arr[3])
                occupy[j, i, 1] = parse(F64, arr[4])
                occupy[j, i, 2] = parse(F64, arr[5])
            end
        end
    end

    # close the iostream
    close(fin)

    # return the desired arrays
    return enk, occupy
end

"""
    vaspio_projs(f::String)

Reading vasp's LOCPROJ file, return raw projector matrix. Here `f` means
only the directory that contains LOCPROJ
"""
function vaspio_projs(f::String)
    # open the iostream
    fin = open(joinpath(f, "LOCPROJ"), "r")

    # extract number of spins (nspin), number of k-points (nkpt),
    # number of bands (nband), and number of projectors (nproj)
    nspin, nkpt, nband, nproj = parse.(I64, line_to_array(fin)[1:4])

    # extract raw information about projectors
    sites = zeros(I64, nproj)
    descs = fill("", nproj)
    for i = 1:nproj
        arr = line_to_array(fin)
        sites[i] = parse(I64, arr[2])
        # get rid of the ":" char
        descs[i] = replace(arr[end], ":" => "")
    end

    # try to build PrTrait struct. the raw information about projectors
    # are encapsulated in it
    PT = PrTrait[]
    for i = 1:nproj
        push!(PT, PrTrait(sites[i], descs[i]))
    end

    # try to split these projectors into groups
    # at first, we collect the tuple (site, l) for all projectors
    site_l = Tuple[]
    for i in eachindex(PT)
        push!(site_l, (PT[i].site, PT[i].l))
    end
    # second, we figure out the unique (site, l)
    union!(site_l)
    # third, we create a array of PrGroup struct (except for site and
    # l, most of its member variables need to be corrected). note
    # that for a given PrGroup, the projectors indexed by PrGroup.Pr
    # should share the same site and l
    PG = PrGroup[]
    for i in eachindex(site_l)
        push!(PG, PrGroup(site_l[i]...))
    end
    # fourth, for each PrGroup, we scan all of the projectors to find
    # out those with correct site and l; record their indices; and
    # save them at PrGroup.Pr array
    for i in eachindex(PG)
        site, l = PG[i].site, PG[i].l
        PG[i].Pr = findall(x -> (x.site, x.l) === (site, l), PT)
    end
    # finally, check
    @assert nproj === sum(x -> length(x.Pr), PG)

    # create arrays
    chipsi = zeros(C64, nproj, nband, nkpt, nspin)

    # read in raw projector data
    readline(fin)
    for s = 1:nspin
        for k = 1:nkpt
            for b = 1:nband
                # extract some indices information
                arr = line_to_array(fin)
                _spin = parse(I64, arr[2])
                _kpt = parse(I64, arr[3])
                _band = parse(I64, arr[4])

                # check consistency of parameters
                @assert _spin === s
                @assert _kpt === k
                @assert _band === b

                # parse the input data
                for p = 1:nproj
                    _re, _im = parse.(F64, line_to_array(fin)[2:3])
                    chipsi[p, b, k, s] = _re + _im * im
                end

                # skip one empty line
                readline(fin)
            end
        end
    end

    # close the iostream
    close(fin)

    # return the desired arrays
    return PT, PG, chipsi
end

"""
    vaspio_fermi(f::String)

Reading vasp's DOSCAR file, return the fermi level. Here `f` means
only the directory that contains DOSCAR
"""
function vaspio_fermi(f::String)
    # open the iostream
    fin = open(joinpath(f, "DOSCAR"), "r")

    # skip five empty lines
    for i = 1:5
        readline(fin)
    end

    # extract the fermi level
    fermi = parse(F64, line_to_array(fin)[4])

    # close the iostream
    close(fin)

    # return the desired data
    return fermi
end

"""
    vaspio_charge(f::String)

Reading vasp's CHGCAR file, return the charge density. Here `f` means
only the directory that contains CHGCAR
"""
function vaspio_charge(f::String) end
