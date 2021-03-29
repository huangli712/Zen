#
# Project : Pansy
# Source  : vasp.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/03/29
#

#
# Driver Functions
#

"""
    vasp_adaptor(D::Dict{Symbol,Any})

Adaptor support for vasp code. It will read the output files of vasp
code and then fulfill the `DFTData` dict (i.e `D`).

The following vasp's files are needed: `POSCAR`, `IBZKPT`, `EIGENVAL`,
`LOCPROJ`, and `DOSCAR`.

See also: [`plo_adaptor`](@ref), [`ir_adaptor`](@ref), [`adaptor_exec`](@ref).
"""
function vasp_adaptor(D::Dict{Symbol,Any})
    # V01: Print the header
    println("  Adaptor : VASP")

    # V02: Read in lattice structure
    println("  Get lattice")
    D[:latt] = vaspio_lattice(pwd())

    # V03: Read in kmesh and the corresponding weights
    println("  Get kmesh")
    println("  Get weight")
    D[:kmesh], D[:weight] = vaspio_kmesh(pwd())

    # V04: Read in band structure and the corresponding occupancies
    println("  Get enk")
    println("  Get occupy")
    D[:enk], D[:occupy] = vaspio_eigen(pwd())

    # V05: Read in raw projectors, traits, and groups
    println("  Get projector")
    D[:PT], D[:PG], D[:chipsi] = vaspio_projs(pwd())

    # V06: Read in fermi level
    println("  Get fermi level")
    D[:fermi] = vaspio_fermi(pwd())

    # V07: Read in tetrahedron data if they are available
    if get_d("smear") === "tetra"
        println("  Get tetrahedron")
        D[:volt], D[:itet] = vaspio_tetra(pwd())
    end

    # V08: Print the footer
    println()
end

"""
    vasp_init(it::IterInfo)

Check the runtime environment of vasp, prepare necessary input files.

See also: [`vasp_exec`](@ref), [`vasp_save`](@ref).
"""
function vasp_init(it::IterInfo)
    # Prepare essential input files
    #
    # Copy POTCAR and POSCAR
    cp("../POTCAR", joinpath(pwd(), "POTCAR"), force = true)
    cp("../POSCAR", joinpath(pwd(), "POSCAR"), force = true)
    #
    # How about INCAR
    if it.dmft_cycle == 0
        # Generate INCAR automatically
        vasp_incar(it.dft_fermi)
    else
        # Maybe we need to update INCAR file here
        vasp_incar(it.dmft_fermi)
    end
    #
    # Well, perhaps we need to generate the KPOINTS file by ourselves.
    if get_d("kmesh") === "file"
        vasp_kpoints()
    end

    # Check essential input files
    flist = ("INCAR", "POSCAR", "POTCAR")
    for i in eachindex(flist)
        filename = flist[i]
        if !isfile(filename)
            error("Please make sure the file $filename is available")
        end
    end
end

"""
    vasp_exec(it::IterInfo)

Execute the vasp program.

See also: [`vasp_init`](@ref), [`vasp_save`](@ref).
"""
function vasp_exec(it::IterInfo)
    # Print the header
    println("  Engine : VASP")

    # Get the home directory of vasp
    dft_home = query_dft("vasp")

    # Determine mpi prefix (whether the vasp is executed sequentially)
    mpi_prefix = inp_toml("../MPI.toml", "dft", false)
    numproc = parse(I64, line_to_array(mpi_prefix)[3])
    println("  Para : Using $numproc processors")

    # Select suitable vasp program
    if get_d("lspinorb")
        vasp_exec = "$dft_home/vasp_ncl"
    else
        vasp_exec = "$dft_home/vasp_std"
    end
    @assert isfile(vasp_exec)
    println("  Exec : $vasp_exec")

    # Assemble command
    if isnothing(mpi_prefix)
        vasp_cmd = vasp_exec
    else
        vasp_cmd = split("$mpi_prefix $vasp_exec", " ")
    end

    # Launch it, the terminal output is redirected to vasp.out
    run(pipeline(`$vasp_cmd`, stdout = "vasp.out"))

    # Print the footer
    println()
end

"""
    vasp_save(it::IterInfo)

Backup the output files of vasp if necessary. Furthermore, the fermi level
in `IterInfo` struct is also updated (`IterInfo.dft_fermi`).

See also: [`vasp_init`](@ref), [`vasp_exec`](@ref).
"""
function vasp_save(it::IterInfo)
    # Store the data files
    #
    # Create list of files
    fl = ["INCAR", "OUTCAR", "CHGCAR", "vasp.out", "vasprun.xml"]
    #
    # Go through the file list, backup the files one by one.
    for i in eachindex(fl)
        f = fl[i]
        cp(f, "$f.$(it.dmft_cycle).$(it.dmft2_iter)", force = true)
    end

    # Anyway, the fermi level is extracted from DOSCAR, and its value
    # will be saved at IterInfo.dft_fermi.
    it.dft_fermi = vaspio_fermi(pwd())
end

#
# Service Functions (Group A)
#

"""
    vasp_incar(fermi::F64)

Generate an `INCAR` file. It will be used only when the DFT engine
is vasp.

See also: [`vasp_kpoints`](@ref).
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
        @case "mp2"
            write(ios, "ISMEAR   = 2 \n")
            break

        @case "mp1"
            write(ios, "ISMEAR   = 1 \n")
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

        # If kmesh == "file", then vasp_kpoints() will be used to
        # generate the KPOINTS file.
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
    lspins = get_d("lspins")
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
    ewidth = 4.0 # A magic number
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

#
# Remarks:
#
# The parameter `NBANDS` is quite important. Sometimes its default value
# is too small. The adaptor will fail to generate reasonable projectors.
# At this case, you will see an error thrown by the try_diag() function.
# The solution is simple, i.e., increasing `NBANDS` a bit.
#
# The current algorithm is suitable for paramagnetic systems. It has not
# been tested for magnetically ordered materials.
#
    # For number of bands
    nbands = vaspio_nband(pwd())
    write(ios, "NBANDS   = $nbands \n")

    # Close the iostream
    close(ios)
end

"""
    vasp_kpoints(mp_scheme::Bool = true, n::I64 = 9)

Generate a valid `KPOINTS` file for vasp.

See also: [`vasp_incar`](@ref).
"""
function vasp_kpoints(mp_scheme::Bool = true, n::I64 = 9)
    # If the `KPOINTS` file is available, we do nothing or else we will
    # try to create a new one.
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

Check the essential output files by vasp. Here `f` means only the
directory that contains the desired files.

See also: [`adaptor_init`](@ref).
"""
function vasp_files(f::String)
    fl = ["POSCAR", "IBZKPT", "EIGENVAL", "LOCPROJ", "DOSCAR", "CHGCAR"]
    for i in eachindex(fl)
        @assert isfile( joinpath(f, fl[i]) )
    end
end

"""
    vasp_files()

Check the essential output files by vasp in the current directory.

See also: [`adaptor_init`](@ref).
"""
vasp_files() = vasp_files(pwd())

#
# Service Functions (Group B)
#

"""
    vaspio_nband(f::String)

Reading vasp's `POSCAR` and `POTCAR` files, evaluating number of bands. It
will be used to create the `INCAR` file. Here `f` means only the directory
that contains `POSCAR` and `POTCAR`.

See also: [`vasp_incar`](@ref), [`vaspio_valence`](@ref), [`vaspio_lattice`](@ref).
"""
function vaspio_nband(f::String)
    # Extract crystallography information from `POSCAR`
    latt = vaspio_lattice(f)

    # Extract number of valence electrons from `POTCAR`
    zval = vaspio_valence(f)

    # Sanity check
    nsort, _ = size(latt.sorts)
    @assert nsort === length(zval)

    # Evaluate number of valence electrons in total
    nelect = sum(@. latt.sorts[:, 2] * zval)

#
# Remarks:
#
# In vasp, the `NBANDS` parameter is determined automatically. According
# to the wiki of vasp, it is equal to nelect / 2 + latt.natom / 2 for
# paramagnetic case by default. However, it may be insufficient for
# generating projectors. For example, for SrVO3, the default `NBANDS`
# parameter is only 20. It is too small to determine the five V-3d
# projectors. It would be better to increase it to 30.
#
# Here, we increase `NBANDS` to `1.6 * NBANDS`. The `1.6` is a magic
# number, you can adjust it by yourself.
#

    # Evaluate number of bands
    if get_d("lspins")
        nband = floor(I64, nelect * 1.2)
    else
        nband = floor(I64, (nelect / 2 + latt.natom / 2) * 1.6)
    end

    # Return the desired nband
    return nband
end

"""
    vaspio_nband()

Reading vasp's `POSCAR` and `POTCAR` files, evaluating number of bands. It
will be used to create the `INCAR` file.

See also: [`vasp_incar`](@ref), [`vaspio_valence`](@ref), [`vaspio_lattice`](@ref).
"""
vaspio_nband() = vaspio_nband(pwd())

"""
    vaspio_valence(f::String)

Reading vasp's `POTCAR` file, return `ZVAL`. Here `f` means only the
directory that contains `POTCAR`.

The information about `ZVAL` will be used to determine `NBANDS` in the
`INCAR` file.

See also: [`vasp_incar`](@ref), [`vaspio_nband`](@ref).
"""
function vaspio_valence(f::String)
    # Open the iostream
    fin = open(joinpath(f, "POTCAR"), "r")

    # Read the POTCAR
    strs = readlines(fin)

    # Extract the lines that contains the `ZVAL`
    filter!(x -> contains(x, "ZVAL"), strs)

    # Extract ZVAL, convert it into float, than save it.
    zval = map(x -> parse(F64, line_to_array(x)[6]), strs)

    # Close the iostream
    close(fin)

    # Return the desired arrays
    return zval
end

"""
    vaspio_valence()

Reading vasp's `POTCAR` file, return `ZVAL`.

The information about `ZVAL` will be used to determine `NBANDS` in the
`INCAR` file.

See also: [`vasp_incar`](@ref), [`vaspio_nband`](@ref).
"""
vaspio_valence() = vaspio_valence(pwd())

"""
    vaspio_procar(f::String)

Reading vasp's `PROCAR` file, extract orbital weight information. Here `f`
means only the directory that contains `PROCAR`.

This function is not invoked directly during the DFT + DMFT iteration. It
is designed for users merely. They can use it to judge which orbitals are
the most relevant, and then apply the obtained information to customize
their case.toml configuration file (specifically, the `window` parameter
in the `dft` block).

See also: [`tools/analyze.jl`](@ref).
"""
function vaspio_procar(f::String)
    # Open the iostream
    fin = open(joinpath(f, "PROCAR"), "r")

    # We have to determine the dimensional parameters at first
    #
    # (1) Determine key parameters: nkpt, nband, and natom.
    readline(fin)
    arr = line_to_array(fin)
    nkpt = parse(I64, arr[4])
    nband = parse(I64, arr[8])
    natom = parse(I64, arr[12])
    seekstart(fin) # Rewind the stream
    #
    # (2) Determine key parameters: norbs.
    readuntil(fin, "ion ")
    arr = line_to_array(fin)
    norbs = length(arr) - 1
    @assert norbs === 9 || norbs === 16
    seekstart(fin) # Rewind the stream
    #
    # (3) Determine key parameters: nspin.
    nspin = 0
    for line in eachline(fin)
        if contains(line, "k-points")
            nspin = nspin + 1
        end
    end
    seekstart(fin) # Rewind the stream
    #
    # (4) Determine whether spin-orbit coupling is activated.
    readuntil(fin, "ion ")
    soc = false
    if natom == 1
        lc = 0
        for line in eachline(fin)
            lc = lc + 1
            if contains(line, "ion ")
                break
            end
        end

        @cswitch lc begin
            @case 3
                soc = false
                break

            @case 6
                soc = true
                break

            @default
                error("Something wrong in PROCAR")
        end
    else # natom > 1
        lc = 0
        for line in eachline(fin)
            if contains(line, "tot ")
                lc = lc + 1
            end
            if contains(line, "ion ")
                break
            end
        end

        @cswitch lc begin
            @case 1
                soc = false
                break

            @case 4
                soc = true
                break

            @default
                error("Something wrong in PROCAR")
        end
    end
    seekstart(fin) # Rewind the stream
    #
    # (5) Additional check. `soc` must be compatible with `nspin`.
    if soc
        @assert nspin === 1
    end
    #
    # (6) Debug
    print("Parameters: ")
    @show norbs, natom, nband, nkpt, nspin, soc

    # Build `str`, which is used to tell the parser how to skip the
    # blank lines for various cases.
    #
    # (1) Init the string
    fstr = ""
    #
    # (2) Single atom vs. multiple atoms
    if natom === 1
        fstr = fstr * "1"
    else
        fstr = fstr * "2"
    end
    #
    # (3) d system vs. f system
    if norbs === 9
        fstr = fstr * "d"
    else
        fstr = fstr * "f"
    end
    #
    # (4) SOC or not
    if soc
        fstr = fstr * "s"
    else
        fstr = fstr * "n"
    end
    #
    # (5) Eight valid cases
    @assert fstr in ["1ds", "1dn", "1fs", "1fn", "2ds", "2dn", "2fs", "2fn"]
    println("PROCAR format: $fstr")

    # Create arrays
    # The `worb` is used to save the raw data, while `oab` is used to
    # save the processed data. The `enk` denotes the eigenvalues and
    # `occ` denotes the occupations.
    worb = zeros(F64, norbs, natom, nband, nkpt, nspin)
    oab = zeros(F64, norbs, natom, nband, nspin)
    enk = zeros(F64, nband, nkpt, nspin)
    occ = zeros(F64, nband, nkpt, nspin)

    # Start to parse the `PROCAR` file
    println("Parsing PROCAR...")

    # Go through each spin orientation
    for s = 1:nspin
        # Blank line
        readline(fin)

        # Go through each k-point
        for k = 1:nkpt
            # Blank lines
            readline(fin)
            readline(fin)

            # Check k-index
            arr = line_to_array(fin)
            @assert k === parse(I64, arr[2])
            println("Finishing spin $s k-point $k")

            # Go through each band for the given k-point and spin
            for b = 1:nband
                # Blank line
                readline(fin)

                # Check band index
                arr = line_to_array(fin)
                @assert b === parse(I64, arr[2])

                # Parse eigenvalues and occupations
                enk[b, k, s] = parse(F64, arr[5])
                occ[b, k, s] = parse(F64, arr[8])

                # Blank lines
                readline(fin)
                readline(fin)

                # Parse weights
                # Go through each atom and orbital
                for a = 1:natom
                    arr = line_to_array(fin)
                    @assert parse(I64, arr[1]) === a
                    @assert norbs === length(arr) - 2
                    worb[:, a, b, k, s] = parse.(F64, arr[2:end-1])
                end

                # Skip useless lines
                @cswitch fstr begin
                    @case "1ds"
                        foreach(x -> readline(fin), 1:1:3)
                        break

                    @case "1dn"
                        break

                    @case "1fs"
                        foreach(x -> readline(fin), 1:1:3)
                        break

                    @case "1fn"
                        break

                    @case "2ds"
                        foreach(x -> readline(fin), 1:1:(3*(1+natom)+1))
                        break

                    @case "2dn"
                        readline(fin)
                        break

                    @case "2fs"
                        foreach(x -> readline(fin), 1:1:(3*(2+natom)+2))
                        break

                    @case "2fn"
                        readline(fin)
                        readline(fin)
                        break

                    @default
                        sorry()
                        break
                end # @cswitch

                # Skip blank line which starts with "ion ..."
                readline(fin)

                # Parse phase factors
                # Go through each atom and orbital
                for a = 1:natom
                    readline(fin)
                end

                # Skip useless lines
                @cswitch fstr begin
                    @case "1ds"
                        break

                    @case "1dn"
                        break

                    @case "1fs"
                        break

                    @case "1fn"
                        break

                    @case "2ds"
                        readline(fin)
                        break

                    @case "2dn"
                        readline(fin)
                        break

                    @case "2fs"
                        readline(fin)
                        readline(fin)
                        break

                    @case "2fn"
                        readline(fin)
                        readline(fin)
                        break

                    @default
                        sorry()
                        break
                end # @cswitch
            end # Loop for bands
        end # Loop for k-points
    end # Loop for spins

    # Close the iostream
    close(fin)

    # Try to build `oab` by k-summation
    for s = 1:nspin
        for b = 1:nband
            for a = 1:natom
                for o = 1:norbs
                    oab[o, a, b, s] = sum(worb[o, a, b, :, s]) / float(nkpt)
                end
            end
        end
    end

    # Return the desired arrays
    return oab, enk, occ
end

"""
    vaspio_procar()

Reading vasp's `PROCAR` file, extract orbital weight information.

This function is not invoked directly during the DFT + DMFT iteration. It
is designed for users merely. They can use it to judge which orbitals are
the most relevant, and then apply the obtained information to customize
their case.toml configuration file (specifically, the `window` parameter
in the `dft` block).

See also: [`tools/analyze.jl`](@ref).
"""
vaspio_procar() = vaspio_procar(pwd())

"""
    vaspio_lattice(f::String)

Reading vasp's `POSCAR` file, return crystallography information. Here `f`
means only the directory that contains `POSCAR`.

See also: [`Lattice`](@ref), [`irio_lattice`](@ref).
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

    # Now all the parameters are ready, we would like to create
    # `Lattice` struct here.
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
    vaspio_lattice()

Reading vasp's `POSCAR` file, return crystallography information.

See also: [`Lattice`](@ref), [`irio_lattice`](@ref).
"""
vaspio_lattice() = vaspio_lattice(pwd())

"""
    vaspio_kmesh(f::String)

Reading vasp's `IBZKPT` file, return `kmesh` and `weight`. Here `f` means
only the directory that contains `IBZKPT`.

See also: [`vaspio_tetra`](@ref), [`irio_kmesh`](@ref).
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
    vaspio_kmesh()

Reading vasp's `IBZKPT` file, return `kmesh` and `weight`.

See also: [`vaspio_tetra`](@ref), [`irio_kmesh`](@ref).
"""
vaspio_kmesh() = vaspio_kmesh(pwd())

"""
    vaspio_tetra(f::String)

Reading vasp's `IBZKPT` file, return tetrahedra information. Here `f`
means only the directory that contains `IBZKPT`.

See also: [`vaspio_kmesh`](@ref), [`irio_tetra`](@ref).
"""
function vaspio_tetra(f::String)
    # Open the iostream
    fin = open(joinpath(f, "IBZKPT"), "r")

    # Extract number of k-points
    readline(fin)
    nkpt = parse(I64, readline(fin))
    readline(fin)

    # Read in the k-points and their weights
    #
    # Skip nkpt lines
    for i = 1:nkpt
        readline(fin)
    end

    # Read in the tetrahedra information
    #
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
    vaspio_tetra()

Reading vasp's `IBZKPT` file, return tetrahedra information.

See also: [`vaspio_kmesh`](@ref), [`irio_tetra`](@ref).
"""
vaspio_tetra() = vaspio_tetra(pwd())

"""
    vaspio_eigen(f::String)

Reading vasp's `EIGENVAL` file, return energy band information. Here `f`
means only the directory that contains `EIGENVAL`.

See also: [`irio_eigen`](@ref).
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
# Here we provide two implementations to read the eigenvalues. The first
# implementation is somewhat tedious, so we don't use it. It seems that
# the second implementation looks quite graceful.
#

#
# Implementation 1
#
            #if nspin === 1 # for spin unpolarized case
            #    enk[j, i, 1] = parse(F64, arr[2])
            #    occupy[j, i, 1] = parse(F64, arr[3])
            #end
            #
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
    vaspio_eigen()

Reading vasp's `EIGENVAL` file, return energy band information.

See also: [`irio_eigen`](@ref).
"""
vaspio_eigen() = vaspio_eigen(pwd())

"""
    vaspio_projs(f::String)

Reading vasp's `LOCPROJ` file, return raw projector matrix. Here `f` means
only the directory that contains `LOCPROJ`.

See also: [`irio_projs`](@ref).
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
    # Second, we figure out the unique (site, l) pairs. Here, we use the
    # union! function, which is much more efficient than the unique! one.
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
    for s = 1:nspin
        for k = 1:nkpt
            for b = 1:nband
                # Skip two empty lines
                readline(fin)
                readline(fin)

                # Parse the data line by line
                for p = 1:nproj
                    chipsi[p, b, k, s] = line_to_cmplx(fin)
                end
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
    vaspio_projs()

Reading vasp's `LOCPROJ` file, return raw projector matrix.

See also: [`irio_projs`](@ref).
"""
vaspio_projs() = vaspio_projs(pwd())

"""
    vaspio_fermi(f::String)

Reading vasp's `DOSCAR` file, return the fermi level. Here `f` means
only the directory that contains `DOSCAR`.

See also: [`irio_fermi`](@ref).
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
    vaspio_fermi()

Reading vasp's `DOSCAR` file, return the fermi level.

See also: [`irio_fermi`](@ref).
"""
vaspio_fermi() = vaspio_fermi(pwd())

"""
    vaspio_charge(f::String)

Reading vasp's `CHGCAR` file, return the charge density. Here `f` means
only the directory that contains `CHGCAR`.

See also: [`irio_charge`](@ref).
"""
function vaspio_charge(f::String) end

"""
    vaspio_charge()

Reading vasp's `CHGCAR` file, return the charge density.

See also: [`irio_charge`](@ref).
"""
vaspio_charge() = vaspio_charge(pwd())
