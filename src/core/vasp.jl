#
# project : pansy
# source  : vasp.jl
# author  : Li Huang (lihuang.dmft@gmail.com)
# status  : unstable
# comment :
#
# last modified: 2021/01/26
#

#
# Driver Functions
#

"""
    vasp_adaptor()

Adaptor support for vasp code. It will read the output files of vasp
code and then fulfill the DFTData dict.
"""
function vasp_adaptor()
    # S01: Print the header
    println("  < VASP Adaptor >")

    # S02: Read in lattice structure
    println("    Get Lattice")
    DFTData["latt"] = vaspio_lattice(pwd())

    # S03: Read in kmesh and the corresponding weights
    println("    Get Kmesh")
    println("    Get Weight")
    DFTData["kmesh"], DFTData["weight"] = vaspio_kmesh(pwd())

    # S04: Read in tetrahedron data if they are available
    if get_d("smear") === "tetra"
        println("    Get Tetrahedron")
        DFTData["volt"], DFTData["itet"] = vaspio_tetra(pwd())
    end

    # S05: Read in band structure and the corresponding occupancies
    println("    Get Enk")
    println("    Get Occupy")
    DFTData["enk"], DFTData["occupy"] = vaspio_eigen(pwd())

    # S06: Read in raw projectors, traits, and groups
    println("    Get Projector (Trait and Group)")
    DFTData["PT"], DFTData["PG"], DFTData["chipsi"] = vaspio_projs(pwd())

    # S07: Read in fermi level
    println("    Get Fermi Level")
    KohnShamData["fermi"] = vaspio_fermi(pwd())
end

"""
    vasp_init(it::IterInfo)

Check the runtime environment of vasp, prepare necessary input files.
"""
function vasp_init(it::IterInfo)
    # Prepare essential input files
    if it.dmft_cycle == 0
        # Copy POTCAR and POSCAR
        cp("../POTCAR", joinpath(pwd(), "POTCAR"), force = true)
        cp("../POSCAR", joinpath(pwd(), "POSCAR"), force = true)

        # Generate INCAR automatically
        vasp_incar(it._dft_fermi)
    end

    # Check essential input files
    if it.dmft_cycle >= 1
        flist = ("INCAR", "POSCAR", "POTCAR")
        for i in eachindex(flist)
            filename = flist[i]
            if !isfile(filename)
                error("Please make sure the file $filename is available")
            end
        end

        # Maybe we need to update INCAR file here
        vasp_incar(it.dmft_fermi)
    end

    # Well, perhaps we need to generate the KPOINTS file by ourselves
    if get_d("kmesh") === "file"
        vasp_kpoints()
    end
end

"""
    vasp_run(it::IterInfo)

Execute the vasp program.
"""
function vasp_run(it::IterInfo)
    # Get the home directory of vasp
    dft_home = query_dft()

    # Determine mpi prefix (whether the vasp is executed sequentially)
    mpi_prefix = inp_toml("../MPI.toml", "dft", false)

    # Select suitable vasp program
    if get_d("lspinorb")
        vasp_exec = "$dft_home/vasp_ncl"
    else
        vasp_exec = "$dft_home/vasp_std"
    end

    # Assemble command
    if isnothing(mpi_prefix)
        vasp_cmd = vasp_exec
    else
        vasp_cmd = split("$mpi_prefix $vasp_exec", " ")
    end

    # Launch it, the terminal output is redirected to vasp.out
    run(pipeline(`$vasp_cmd`, stdout = "vasp.out"))
end

"""
    vasp_save(it::IterInfo)

Backup the output files of vasp if necessary. Furthermore, the fermi level
in IterInfo struct is also updated (IterInfo._dft_fermi)
"""
function vasp_save(it::IterInfo)
    # Store the data files
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

    # Anyway, the fermi level is extracted from DOSCAR, and its value will
    # be saved at IterInfo._dft_fermi.
    it._dft_fermi = vaspio_fermi(pwd())
end

#
# Service Functions (Group A)
#

"""
    vasp_incar(fermi::F64)

Generate an INCAR file. It will be used only when the DFT engine is vasp.
"""
function vasp_incar(fermi::F64)
    # Open the iostream
    ios = open("INCAR", "w")

    # Standard part
    case = get_c("case")
    write(ios, "System   = $case \n")
    write(ios, "PREF     = Accurate \n")
    write(ios, "EDIFF    = 1E-8 \n")
    write(ios, "ALGO     = Normal \n")
    write(ios, "LASPH    = .TRUE. \n")
    write(ios, "LMAXMIX  = 6 \n")
    write(ios, "NCORE    = 4 \n")

    # Customize your INCAR according to the case.toml
    #
    # For smearing
    smear = get_d("smear")
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

    # For kmesh density
    #
    # If kmesh == "file", then vasp_kpoints() will be used to generate
    # the KPOINTS file.
    kmesh = get_d("kmesh")
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

        @default # Very coarse kmesh
            write(ios, "KSPACING = 0.5 \n")
            break
    end

    # For magnetic moment
    magmom = get_d("magmom")
    if !isa(magmom, Missing)
        write(ios, "MAGMOM   = $magmom \n")
    end

    # For symmetry
    lsymm = get_d("lsymm")
    if lsymm
        write(ios, "ISYM     = 2 \n")
    else # Ignore the symmetry completely
        write(ios, "ISYM     =-1 \n")
    end

    # For spin polarizations
    #
    # If spin-orbit coupling is on, spin orientations must be polarized.
    lspins = get_d("lspins") || get_d("lspinorb")
    if lspins
        write(ios, "ISPIN    = 2 \n")
    else
        write(ios, "ISPIN    = 1 \n")
    end

    # For spin-orbit coupling
    lspinorb = get_d("lspinorb")
    if lspinorb
        write(ios, "LSORBIT  = .TRUE. \n")
    else
        write(ios, "LSORBIT  = .FALSE. \n")
    end

    # For optimized projectors
    ewidth = get_d("ewidth")
    loptim = get_d("loptim")
    if !isa(ewidth, Missing) && !isa(loptim, Missing)
        if loptim
            write(ios, "LORBIT   = 14 \n")
            emin = fermi - ewidth
            write(ios, "EMIN     = $emin \n")
            emax = fermi + ewidth
            write(ios, "EMAX     = $emax \n")
        end
    end

    # For local orbitals and projectors
    lproj = get_d("lproj")
    sproj = get_d("sproj")
    if !isa(lproj, Missing) && !isa(sproj, Missing)
        if lproj
            for p in eachindex(sproj)
                str = sproj[p]
                write(ios, "LOCPROJ  = $str \n")
            end
        end
    end

    # Close the iostream
    close(ios)
end

"""
    vasp_kpoints(mp_scheme::Bool = true, n::I64 = 9)

Generate a valid KPOINTS file for vasp.
"""
function vasp_kpoints(mp_scheme::Bool = true, n::I64 = 9)
    # If the KPOINTS file is available, we do nothing or else we will
    # try to create a new one
    if !isfile("KPOINTS")
        # Open the iostream
        ios = open("KPOINTS", "w")

        # Write the body
        write(ios, "Automatic K-mesh Generation\n")
        write(ios, "0 \n")
        if mp_scheme
            write(ios, "Monkhorst-Pack\n")
        else
            write(ios, "Gamma\n")
        end
        write(ios, " $n  $n  $n\n")
        write(ios, " 0  0  0\n")

        # Close the iostream
        close(ios)
    end
end

"""
    vasp_files(f::String)

Check the essential output files by vasp. Here `f` means only the directory
that contains the desired files.
"""
function vasp_files(f::String)
    @assert isfile(joinpath(f, "POSCAR")) &&
            isfile(joinpath(f, "IBZKPT")) &&
            isfile(joinpath(f, "DOSCAR")) &&
            isfile(joinpath(f, "CHGCAR")) &&
            isfile(joinpath(f, "LOCPROJ")) &&
            isfile(joinpath(f, "EIGENVAL"))
end
vasp_files() = vasp_files(pwd())

#
# Service Functions (Group B)
#

"""
    vaspio_lattice(f::String)

Reading vasp's POSCAR file, return crystallography information. Here `f`
means only the directory that contains POSCAR.
"""
function vaspio_lattice(f::String)
    # Open the iostream
    fin = open(joinpath(f, "POSCAR"), "r")

    # Get the case
    _case = string(strip(readline(fin)))

    # Get the scaling factor
    scale = parse(F64, readline(fin))

    # Get the lattice vectors
    lvect = zeros(F64, 3, 3)
    lvect[1, :] = parse.(F64, line_to_array(fin))
    lvect[2, :] = parse.(F64, line_to_array(fin))
    lvect[3, :] = parse.(F64, line_to_array(fin))

    # Get the symbol list
    symbols = line_to_array(fin)

    # Get the number of sorts of atoms
    nsort = length(symbols)

    # Get the number list
    numbers = parse.(I64, line_to_array(fin))

    # Get the total number of atoms
    natom = sum(numbers)

    # Now all the parameters are ready
    # We would like to create Lattice struct here
    latt = Lattice(_case, scale, nsort, natom)

    # Update latt using the available data
    latt.lvect = lvect
    for i = 1:nsort
        latt.sorts[i, 1] = string(symbols[i])
        latt.sorts[i, 2] = numbers[i]
    end

    # Get the atom list
    k = 0
    for i = 1:nsort
        for j = 1:numbers[i]
            k = k + 1
            latt.atoms[k] = symbols[i]
        end
    end
    # Sanity check
    @assert k === natom

    # Get the coordinates of atoms
    readline(fin)
    for i = 1:natom
        latt.coord[i, :] = parse.(F64, line_to_array(fin)[1:3])
    end

    # Close the iostream
    close(fin)

    # Return the desired struct
    return latt
end

"""
    vaspio_kmesh(f::String)

Reading vasp's IBZKPT file, return kmesh and weight. Here `f` means
only the directory that contains IBZKPT.
"""
function vaspio_kmesh(f::String)
    # Open the iostream
    fin = open(joinpath(f, "IBZKPT"), "r")

    # Extract number of k-points
    readline(fin)
    nkpt = parse(I64, readline(fin))
    readline(fin)

    # Create arrays
    kmesh = zeros(F64, nkpt, 3)
    weight = zeros(F64, nkpt)

    # Read in the k-points and their weights
    for i = 1:nkpt
        arr = parse.(F64, line_to_array(fin))
        kmesh[i, 1:3] = arr[1:3]
        weight[i] = arr[4]
    end

    # Close the iostream
    close(fin)

    # Return the desired arrays
    return kmesh, weight
end

"""
    vaspio_tetra(f::String)

Reading vasp's IBZKPT file, return tetrahedra information. Here `f` means
only the directory that contains IBZKPT.
"""
function vaspio_tetra(f::String)
    # Open the iostream
    fin = open(joinpath(f, "IBZKPT"), "r")

    # Extract number of k-points
    readline(fin)
    nkpt = parse(I64, readline(fin))
    readline(fin)

    # Read in the k-points and their weights
    # Skip nkpt lines
    for i = 1:nkpt
        readline(fin)
    end

    # Read in the tetrahedra information
    # Skip one empty line
    readline(fin)

    # Extract total number of tetrahedra and volume of a tetrahedron
    arr = line_to_array(fin)
    ntet = parse(I64, arr[1])
    volt = parse(F64, arr[2])

    # Create arrays
    itet = zeros(I64, ntet, 5)

    # Parse the input tetrahedra information
    for t = 1:ntet
        itet[t, :] = parse.(I64, line_to_array(fin))
    end

    # Close the iostream
    close(fin)

    # Return the desired arrays
    return volt, itet
end

"""
    vaspio_eigen(f::String)

Reading vasp's EIGENVAL file, return energy band information. Here `f`
means only the directory that contains EIGENVAL.
"""
function vaspio_eigen(f::String)
    # Open the iostream
    fin = open(joinpath(f, "EIGENVAL"), "r")

    # Determine number of spins
    nspin = parse(I64, line_to_array(fin)[end])
    @assert nspin === 1 || nspin === 2

    # Skip for lines
    for i = 1:4
        readline(fin)
    end

    # Read in some key parameters: nelect, nkpt, nbands
    nelect, nkpt, nband = parse.(I64, line_to_array(fin))

    # Create arrays
    enk = zeros(F64, nband, nkpt, nspin)
    occupy = zeros(F64, nband, nkpt, nspin)

    # Read in the energy bands and the corresponding occupations
    for i = 1:nkpt
        readline(fin)
        readline(fin)
        for j = 1:nband
            arr = line_to_array(fin)

#
# Remarks:
#
# Here we provide two implementations. The first implementation is somewhat
# tedious, so we don't use it. It seems that the second implementation is
# more graceful.
#

#
# Implementation 1
#
            #if nspin === 1 # for spin unpolarized case
            #    enk[j, i, 1] = parse(F64, arr[2])
            #    occupy[j, i, 1] = parse(F64, arr[3])
            #end
            #if nspin === 2 # for spin polarized case
            #    enk[j, i, 1] = parse(F64, arr[2])
            #    enk[j, i, 2] = parse(F64, arr[3])
            #    occupy[j, i, 1] = parse(F64, arr[4])
            #    occupy[j, i, 2] = parse(F64, arr[5])
            #end

#
# Implementation 2
#
            for s = 1:nspin
                enk[j, i, s] = parse(F64, arr[s+1])
                occupy[j, i, s] = parse(F64, arr[s+1+nspin])
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
only the directory that contains LOCPROJ.
"""
function vaspio_projs(f::String)
    # Open the iostream
    fin = open(joinpath(f, "LOCPROJ"), "r")

    # Extract number of spins (nspin), number of k-points (nkpt),
    # number of bands (nband), and number of projectors (nproj).
    nspin, nkpt, nband, nproj = parse.(I64, line_to_array(fin)[1:4])
    @assert nspin === 1 || nspin === 2

    # Extract raw information about projectors
    sites = zeros(I64, nproj)
    descs = fill("", nproj)
    for i = 1:nproj
        arr = line_to_array(fin)
        sites[i] = parse(I64, arr[2])
        # Get rid of the ":" char
        descs[i] = replace(arr[end], ":" => "")
    end

    # Try to build PrTrait struct. The raw information about projectors
    # should be encapsulated in it.
    PT = PrTrait[]
    for i = 1:nproj
        push!(PT, PrTrait(sites[i], descs[i]))
    end

    # Try to split these projectors into groups.
    #
    # At first, we collect the tuple (site, l) for all projectors.
    site_l = Tuple[]
    for i in eachindex(PT)
        push!(site_l, (PT[i].site, PT[i].l))
    end
    #
    # Second, we figure out the unique (site, l) pairs.
    union!(site_l)
    #
    # Third, we create a array of PrGroup struct (except for site and
    # l, most of its member variables need to be corrected). Note
    # that for a given PrGroup, the projectors indexed by PrGroup.Pr
    # should share the same site and l.
    PG = PrGroup[]
    for i in eachindex(site_l)
        push!(PG, PrGroup(site_l[i]...))
    end
    #
    # Fourth, for each PrGroup, we scan all of the projectors to find
    # out those with correct site and l; record their indices; and
    # save them at PrGroup.Pr array.
    for i in eachindex(PG)
        site, l = PG[i].site, PG[i].l
        PG[i].Pr = findall(x -> (x.site, x.l) === (site, l), PT)
    end
    #
    # Finally, check correctness
    @assert nproj === sum(x -> length(x.Pr), PG)

    # Create arrays
    chipsi = zeros(C64, nproj, nband, nkpt, nspin)

    # Read in raw projector data
    readline(fin)
    for s = 1:nspin
        for k = 1:nkpt
            for b = 1:nband
                # Extract some indices information
                arr = line_to_array(fin)
                _spin = parse(I64, arr[2])
                _kpt = parse(I64, arr[3])
                _band = parse(I64, arr[4])

                # Check consistency of parameters
                @assert _spin === s
                @assert _kpt === k
                @assert _band === b

                # Parse the input data
                for p = 1:nproj
                    _re, _im = parse.(F64, line_to_array(fin)[2:3])
                    chipsi[p, b, k, s] = _re + _im * im
                end

                # Skip one empty line
                readline(fin)
            end
        end
    end

    # Close the iostream
    close(fin)

    # Return the desired arrays
    # Note: PG should be further setup at plo_group() function.
    return PT, PG, chipsi
end

"""
    vaspio_fermi(f::String)

Reading vasp's DOSCAR file, return the fermi level. Here `f` means
only the directory that contains DOSCAR.
"""
function vaspio_fermi(f::String)
    # Open the iostream
    fin = open(joinpath(f, "DOSCAR"), "r")

    # Skip five empty lines
    for i = 1:5
        readline(fin)
    end

    # Extract the fermi level
    fermi = parse(F64, line_to_array(fin)[4])

    # Close the iostream
    close(fin)

    # Return the desired data
    return fermi
end

"""
    vaspio_charge(f::String)

Reading vasp's CHGCAR file, return the charge density. Here `f` means
only the directory that contains CHGCAR.
"""
function vaspio_charge(f::String) end
