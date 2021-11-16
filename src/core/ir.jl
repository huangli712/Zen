#
# Project : Pansy
# Source  : ir.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/11/11
#

#=
### *Driver Functions*
=#

#=
*Remarks* :

In the `ir_adaptor()` funtion, some writting steps are optional because
the tetrahedron information might be absent.
=#

"""
    ir_adaptor(D::Dict{Symbol,Any})

Write the Kohn-Sham dataset to specified files using the IR format. Note
that the Kohn-Sham dataset are encapsulated in the `D` dict.

See also: [`plo_adaptor`](@ref), [`wannier_adaptor`](@ref).
"""
function ir_adaptor(D::Dict{Symbol,Any})
    # I01: Check the validity of the `D` dict
    key_list = [:MAP, :PG, :PW,
                :latt, :kmesh, :weight,
                :enk, :occupy, :Fchipsi, :fermi, :chipsi]
    for k in key_list
        @assert haskey(D, k)
    end

    # I02: Print the header
    println("Adaptor : IR")
    println("Try to write the processed Kohn-Sham dataset with IR format")
    println("Current directory: ", pwd())
    println("Store essential parameters")

    # I03: Write important parameters
    #
    # Key dimensional parameters
    irio_params(pwd(), D)
    #
    # Mapping struct
    irio_maps(pwd(), D[:MAP])
    #
    # PrGroup struct
    irio_groups(pwd(), D[:PG])
    #
    # PrWindow struct
    irio_windows(pwd(), D[:PW])

    # I04: Write lattice structure
    irio_lattice(pwd(), D[:latt])

    # I05: Write kmesh and the corresponding weights
    irio_kmesh(pwd(), D[:kmesh], D[:weight])

    # I06: Write band structure and the corresponding occupancies
    irio_eigen(pwd(), D[:enk], D[:occupy])

    # I07: Write normalized projectors
    irio_projs(pwd(), D[:Fchipsi])

    # I08: Write fermi level
    irio_fermi(pwd(), D[:fermi])

    # I09: Write raw projectors
    irio_rawcp(pwd(), D[:chipsi])

    # I10: Check the validity of the `D` dict further (optional)
    if get_d("smear") == "tetra" && !is_qe()
        key_list = [:volt, :itet]
        for k in key_list
            @assert haskey(D, k)
        end
    end

    # I11: Write tetrahedron data if they are available
    if get_d("smear") == "tetra" && !is_qe()
        irio_tetra(pwd(), D[:volt], D[:itet])
    end
end

"""
    ir_save(it::IterInfo)

Backup the files outputed by the adaptor.

See also: [`ir_adaptor`](@ref).
"""
function ir_save(it::IterInfo)
    # Create a list of files that need to be backup
    fir1 = ["params", "maps", "groups", "windows"]
    fir2 = ["lattice", "kmesh", "eigen", "projs", "fermi"]
    file_list = union(fir1, fir2)
    #
    # If tetrahedron data are available
    if get_d("smear") == "tetra" && !is_qe()
        push!(file_list, "tetra")
    end

    # Store the data files
    for i in eachindex(file_list)
        file_src = file_list[i] * ".ir"
        file_dst = file_list[i] * ".ir.$(it.I₃)"
        cp(file_src, file_dst, force = true)
    end
end

"""
    ir_read(f::String)

Read and parse `maps.ir`, `groups.ir`, `windows.ir`, `lattice.ir`,
`kmesh.ir`, `eigen.ir`, `projs.ir`, `fermi.ir`, and `rawcp.ir`. The
data are encapsulated in a dictionary. Here `f` means the directory
where the files as mentioned above are available.

See also: [`ir_adaptor`](@ref).
"""
function ir_read(f::String)
    # Create a dict which is used to store the Kohn-Sham data.
    D = Dict{Symbol,Any}()

    # Get mapping between groups, windows and quantum impurity problems.
    D[:MAP] = irio_maps(f)

    # Get groups of projectors
    D[:PG] = irio_groups(f)

    # Get windows for projectors
    D[:PW] = irio_windows(f)

    # Get crystal structure
    D[:latt] = irio_lattice(f)

    # Get 𝑘-mesh and weights
    D[:kmesh], D[:weight] = irio_kmesh(f)

    # Get eigenvalues and occupations
    D[:enk], D[:occupy] = irio_eigen(f)

    # Get normalized projectors
    D[:Fchipsi] = irio_projs(f)

    # Get fermi level
    D[:fermi] = irio_fermi(f)

    # Get raw projectors
    D[:chipsi] = irio_rawcp(f)

    # Return the desired dict
    return D
end

#=
### *Service Functions* : *Files I/O*
=#

"""
    irio_params(f::String, D::Dict{Symbol,Any})

Write the key parameters extracted from the Kohn-Sham dataset. Here `f`
means only the directory that we want to use.

See also: [`PrGroup`](@ref), [`PrWindow`](@ref).
"""
function irio_params(f::String, D::Dict{Symbol,Any})
    # Extract crystallography information
    _case = D[:latt]._case
    scale = D[:latt].scale
    nsort = D[:latt].nsort
    natom = D[:latt].natom

    # Extract the fermi level
    fermi = D[:fermi]

    # Extract some key parameters
    nband, nkpt, nspin = size(D[:enk])

    # Extract `ntet`, it is optional.
    if haskey(D, :itet)
        ntet, _ = size(D[:itet])
    else
        ntet = 1 # To avoid null array
    end

    # Extract `volt`, it is optional.
    if haskey(D, :volt)
        volt = D[:volt]
    else
        volt = 1.0 # To avoid wrong renormalization
    end

    # Extract `ngrp`
    ngrp, = size(D[:PG])
    #
    # Extract `qdim`, maximum number of projectors among the groups.
    # `size(D[:PG][g].Tr,1)` is actually ndim
    qdim = maximum( [ size(D[:PG][g].Tr,1) for g = 1:ngrp ] )

    # Extract `nwnd`
    nwnd, = size(D[:PW])
    #
    # Extract qbnd, maximum number of bands among the windows.
    qbnd = maximum( [ D[:PW][w].nbnd for w = 1:nwnd ] )
    #
    # Extract xbnd, maximum number of bands included in all windows.
    _bmax = maximum( [ D[:PW][w].bmax for w = 1:nwnd ] )
    _bmin = minimum( [ D[:PW][w].bmin for w = 1:nwnd ] )
    xbnd = _bmax - _bmin + 1
    @assert xbnd ≥ qbnd

    # D[:PW] and D[:PG] should have the same size
    @assert ngrp == nwnd

    # Extract `nsite` and `nmesh`
    nmesh = get_m("nmesh")
    nsite = get_i("nsite")

    # To make sure the validaty of nsite
    @assert nsite ≤ ngrp

    # Output the data
    open(joinpath(f, "params.ir"), "w") do fout
        # Write the header
        println(fout, "# File: params.ir")
        println(fout, "# Data: some necessary parameters")
        println(fout)

        # Write basic parameters
        println(fout, "# Lattice:")
        println(fout, "model -> $_case")
        println(fout, "scale -> $scale")
        println(fout, "nsort -> $nsort")
        println(fout, "natom -> $natom")
        println(fout)

        println(fout, "# Eigen:")
        println(fout, "nband -> $nband")
        println(fout, "nkpt  -> $nkpt")
        println(fout, "nspin -> $nspin")
        println(fout, "fermi -> $fermi")
        println(fout)

        println(fout, "# Tetra:")
        println(fout, "ntet  -> $ntet")
        println(fout, "volt  -> $volt")
        println(fout)

        println(fout, "# Group:")
        println(fout, "ngrp  -> $ngrp")
        println(fout, "qdim  -> $qdim")
        println(fout)

        println(fout, "# Window:")
        println(fout, "nwnd  -> $nwnd")
        println(fout, "qbnd  -> $qbnd")
        println(fout, "xbnd  -> $xbnd")
        println(fout)

        println(fout, "# Sigma:")
        println(fout, "nsite -> $nsite")
        println(fout, "nmesh -> $nmesh")
        println(fout)
    end # END OF IOSTREAM

    # Print some useful information
    println("  > Open and write the file params.ir")
end

"""
    irio_maps(f::String, MAP::Mapping)

Write the information contained in `Mapping` struct. Here `f` means only
the directory that we want to use.

See also: [`Mapping`](@ref).
"""
function irio_maps(f::String, MAP::Mapping)
    # Extract key parameters
    nsite = length(MAP.i_grp)
    ngrp = length(MAP.g_imp)
    nwnd = length(MAP.w_imp)

    # Output the data
    open(joinpath(f, "maps.ir"), "w") do fout
        # Write the header
        println(fout, "# File: maps.ir")
        println(fout, "# Data: some necessary data structures")
        println(fout)
        println(fout, "nsite -> $nsite")
        println(fout, "ngrp  -> $ngrp")
        println(fout, "nwnd  -> $nwnd")
        println(fout)

        # Write the body
        # For i_grp part
        println(fout, "# i_grp")
        foreach(x -> @printf(fout, "%8i", x), MAP.i_grp)
        println(fout)
        println(fout)

        # For i_wnd part
        println(fout, "# i_wnd")
        foreach(x -> @printf(fout, "%8i", x), MAP.i_wnd)
        println(fout)
        println(fout)

        # For g_imp part
        println(fout, "# g_imp")
        foreach(x -> @printf(fout, "%8i", x), MAP.g_imp)
        println(fout)
        println(fout)

        # For w_imp part
        println(fout, "# w_imp")
        foreach(x -> @printf(fout, "%8i", x), MAP.w_imp)
        println(fout)
        println(fout)
    end # END OF IOSTREAM

    # Print some useful information
    println("  > Open and write the file maps.ir")
end

"""
    irio_maps(f::String)

Extract the `Mapping` struct from `maps.ir`. Here `f` means the
directory that this file exists.

See also: [`Mapping`](@ref).
"""
function irio_maps(f::String)
    # Print the header
    println("Open and parse the file maps.ir")

    # Check file's status
    fn = joinpath(f, "maps.ir")
    @assert isfile(fn)

    # Declare a Mapping struct
    MAP = nothing

    # Input the data
    open(fn, "r") do fin
        # Skip the header
        readline(fin)
        readline(fin)
        readline(fin)

        # Extract some dimensional parameters
        nsite = parse(I64, line_to_array(fin)[3])
        ngrp  = parse(I64, line_to_array(fin)[3])
        nwnd  = parse(I64, line_to_array(fin)[3])
        @assert ngrp == nwnd
        readline(fin)

        # For i_grp part
        readline(fin)
        i_grp = parse.(I64, line_to_array(fin))
        @assert length(i_grp) == nsite
        readline(fin)

        # For i_wnd part
        readline(fin)
        i_wnd = parse.(I64, line_to_array(fin))
        @assert length(i_wnd) == nsite
        readline(fin)

        # For g_imp part
        readline(fin)
        g_imp = parse.(I64, line_to_array(fin))
        @assert length(g_imp) == ngrp
        readline(fin)

        # For w_imp part
        readline(fin)
        w_imp = parse.(I64, line_to_array(fin))
        @assert length(w_imp) == nwnd
        readline(fin)

        # Build the MAP struct
        # Call the default constructor
        MAP = Mapping(i_grp, i_wnd, g_imp, w_imp)
    end

    # Print some useful information
    println("  > Number of impurity sites: ", length(MAP.i_grp))
    println("  > Number of groups: ", length(MAP.g_imp))
    println("  > Number of windows: ", length(MAP.w_imp))

    # Return the desired value
    return MAP
end

"""
    irio_groups(f::String, PG::Array{PrGroup,1})

Write the information contained in `PrGroup` struct. Here `f` means only
the directory that we want to use.

See also: [`PrGroup`](@ref).
"""
function irio_groups(f::String, PG::Array{PrGroup,1})
    # Output the data
    open(joinpath(f, "groups.ir"), "w") do fout
        # Write the header
        println(fout, "# File: groups.ir")
        println(fout, "# Data: some necessary data structures")
        println(fout)

        # Write each PrGroup
        ngrp = length(PG)
        #
        println(fout, "ngrp  -> $ngrp")
        println(fout)
        #
        for p in eachindex(PG)
            println(fout, "# PrGroup : $p")
            println(fout, "site  -> $(PG[p].site)")
            println(fout, "l     -> $(PG[p].l)")
            println(fout, "corr  -> $(PG[p].corr)")
            println(fout, "shell -> $(PG[p].shell)")
            println(fout, "ndim  -> $(size(PG[p].Tr,1))")
            println(fout)
        end
    end # END OF IOSTREAM

    # Print some useful information
    println("  > Open and write the file groups.ir")
end

"""
    irio_groups(f::String)

Extract the `PrGroup` struct from `groups.ir`. Here `f` means the
directory that this file exists. Be careful, the returned `PrGroup`
is not completely the same with the true one. For example, its `Pr`
and `Tr` fields are not correct. But it doesn't matter.

See also: [`PrGroup`](@ref).
"""
function irio_groups(f::String)
    # Print the header
    println("Open and parse the file groups.ir")

    # Check file's status
    fn = joinpath(f, "groups.ir")
    @assert isfile(fn)

    # Define array of PrGroup struct
    PG = PrGroup[]

    # Input the data
    open(fn, "r") do fin
        # Skip the header
        readline(fin)
        readline(fin)
        readline(fin)

        # Get number of groups
        ngrp = parse(I64, line_to_array(fin)[3])
        @assert ngrp ≥ 1
        readline(fin)

        # Go through each group
        for p = 1:ngrp
            # Check index of the current group
            _p = parse(I64, line_to_array(fin)[4])
            @assert _p == p

            # Get some key parameters for this group
            site = parse(I64, line_to_array(fin)[3])
            l    = parse(I64, line_to_array(fin)[3])
            corr = parse(Bool, line_to_array(fin)[3])
            shell = line_to_array(fin)[3]
            ndim = parse(I64, line_to_array(fin)[3])
            readline(fin)

            # Create a PrGroup struct and further setup it
            _pg = PrGroup(site, l)
            _pg.corr = corr
            _pg.shell = shell

            # Store this group in PG
            push!(PG, _pg)
            println("  > Extract group [$p]")
        end # END OF P LOOP
    end # END OF IOSTREAM

    # Print the summary
    println("  > Summary of groups:")
    for i in eachindex(PG)
        print("    [ Group $i ]")
        print("  site -> ", PG[i].site)
        print("  l -> ", PG[i].l)
        print("  corr -> ", PG[i].corr)
        print("  shell -> ", PG[i].shell)
        println("  Pr -> ", PG[i].Pr)
    end

    # Return the desired array
    return PG
end

"""
    irio_windows(f::String, PW::Array{PrWindow,1})

Write the information contained in `PrWindow` struct. Here `f` means only
the directory that we want to use.

See also: [`PrWindow`](@ref).
"""
function irio_windows(f::String, PW::Array{PrWindow,1})
    # Output the data
    open(joinpath(f, "windows.ir"), "w") do fout
        # Write the header
        println(fout, "# File: windows.ir")
        println(fout, "# Data: some necessary data structures")
        println(fout)

        # Write each PrWindow
        nwnd = length(PW)
        #
        println(fout, "nwnd  -> $nwnd")
        println(fout)
        #
        for p in eachindex(PW)
            println(fout, "# PrWindow: $p")
            println(fout, "bmin  -> $(PW[p].bmin)")
            println(fout, "bmax  -> $(PW[p].bmax)")
            println(fout, "nbnd  -> $(PW[p].nbnd)")
            println(fout, "kwin  ->")
            #
            nkpt, nspin, ndir = size(PW[p].kwin)
            @assert ndir == 2 # For lower and upper boundaries
            #
            for s = 1:nspin
                for k = 1:nkpt
                    @printf(fout, "%8i %8i %8i %8i\n", k, s, PW[p].kwin[k, s, :]...)
                end
            end
            #
            println(fout)
        end # END OF P LOOP
    end # END OF IOSTREAM

    # Print some useful information
    println("  > Open and write the file windows.ir")
end

"""
    irio_windows(f::String)

Extract the `PrWindow` struct from `windows.ir`. Here `f` means the
directory that this file exists.

See also: [`PrWindow`](@ref).
"""
function irio_windows(f::String)
    # Print the header
    println("Open and parse the file windows.ir")

    # Check file's status
    fn = joinpath(f, "windows.ir")
    @assert isfile(fn)

    # Define array of PrWindow struct
    PW = PrWindow[]

    # Input the data
    open(fn, "r") do fin
        # Skip the header
        readline(fin)
        readline(fin)
        readline(fin)

        # Get number of windows
        nwnd = parse(I64, line_to_array(fin)[3])
        @assert nwnd ≥ 1
        readline(fin)

        # Go through each window
        for p = 1:nwnd
            # Check index of the current window
            _p = parse(I64, line_to_array(fin)[3])
            @assert _p == p

            # Get some key parameters for this window
            bmin = parse(I64, line_to_array(fin)[3])
            bmax = parse(I64, line_to_array(fin)[3])
            nbnd = parse(I64, line_to_array(fin)[3])
            bwin = (bmin, bmax)
            readline(fin)

            # Read data for kwin. All the data will be stored in lines.
            lines = []
            c = 0 # Counter
            str = readline(fin)
            while length(str) > 0
                c  = c + 1 # Increase counter
                push!(lines, str)
                str = readline(fin)
            end

            # Extract dimensional parameters for this PrWindow (nkpt and nspin)
            nkpt, nspin, _, _ = parse.(I64, line_to_array(lines[end]))

            # Create array for kwin
            kwin = zeros(I64, nkpt, nspin, 2)

            # Now we parse the data in lines and fill in kwin
            c = 0 # Counter
            for s = 1:nspin
                for k = 1:nkpt
                    c = c + 1 # Increase counter
                    arr = line_to_array(lines[c])
                    _k, _s = parse.(I64, arr[1:2])
                    @assert _k == k && _s == s
                    kwin[k,s,:] = parse.(I64, arr[3:4])
                end
            end

            # Initialize a PrWindow struct and store it in PW
            push!(PW, PrWindow(kwin,bwin))
            println("  > Extract window [$p]")
        end # END OF P LOOP
    end # END OF IOSTREAM

    # Print the summary
    println("  > Summary of windows:")
    for i in eachindex(PW)
        print("    [ Window $i ]")
        print("  bmin -> ", PW[i].bmin)
        print("  bmax -> ", PW[i].bmax)
        print("  nbnd -> ", PW[i].nbnd)
        println("  bwin -> ", PW[i].bwin)
    end

    # Return the desired array
    return PW
end

"""
    irio_lattice(f::String, latt::Lattice)

Write the lattice information to `lattice.ir` using the IR format. Here
`f` means only the directory that we want to use.

See also: [`vaspio_lattice`](@ref), [`qeio_lattice`](@ref).
"""
function irio_lattice(f::String, latt::Lattice)
    # Print the header
    println("Store essential Kohn-Sham dataset")

    # Extract some key parameters
    _case, scale, nsort, natom = latt._case, latt.scale, latt.nsort, latt.natom

    # Output the data
    open(joinpath(f, "lattice.ir"), "w") do fout
        # Write the header
        println(fout, "# File: lattice.ir")
        println(fout, "# Data: Lattice struct")
        println(fout)
        println(fout, "_case -> $_case")
        println(fout, "scale -> $scale")
        println(fout, "nsort -> $nsort")
        println(fout, "natom -> $natom")
        println(fout)

        # Write the body
        # For sorts part
        println(fout, "[sorts]")
        for i = 1:nsort # Symbols
            @printf(fout, "%8s", latt.sorts[i, 1])
        end
        println(fout)
        for i = 1:nsort # Numbers
            @printf(fout, "%8i", latt.sorts[i, 2])
        end
        println(fout)
        println(fout)

        # For atoms part
        println(fout, "[atoms]")
        for i = 1:natom
            @printf(fout, "%8s", latt.atoms[i])
        end
        println(fout)
        println(fout)

        # For lvect part
        println(fout, "[lvect]")
        for i = 1:3
            @printf(fout, "%16.12f %16.12f %16.12f\n", latt.lvect[i, 1:3]...)
        end
        println(fout)

        # For coord part
        println(fout, "[coord]")
        for i = 1:natom
            @printf(fout, "%16.12f %16.12f %16.12f\n", latt.coord[i, 1:3]...)
        end
    end # END OF IOSTREAM

    # Print some useful information
    println("  > Open and write the file lattice.ir (lattice)")
end

"""
    irio_lattice(f::String)

Extract the lattice information from `lattice.ir`. Here `f` means the
directory that this file exists.

See also: [`Lattice`](@ref).
"""
function irio_lattice(f::String)
    # Print the header
    println("Open and parse the file lattice.ir (lattice)")

    # Check file's status
    fn = joinpath(f, "lattice.ir")
    @assert isfile(fn)

    # Define lattice struct
    latt = nothing

    # Input the data
    open(fn, "r") do fin
        # Skip the header
        readline(fin)
        readline(fin)
        readline(fin)

        # Extract some key parameters
        _case = String(line_to_array(fin)[3])
        scale = parse(F64, line_to_array(fin)[3])
        nsort = parse(I64, line_to_array(fin)[3])
        natom = parse(I64, line_to_array(fin)[3])
        readline(fin)

        # Now all the parameters are ready, we would like to create
        # a `Lattice` struct here.
        latt = Lattice(_case, scale, nsort, natom)

        # For sorts part
        readline(fin)
        latt.sorts[:,1] = String.(line_to_array(fin))
        latt.sorts[:,2] = parse.(I64, line_to_array(fin))
        readline(fin)

        # For atoms part
        readline(fin)
        latt.atoms = String.(line_to_array(fin))
        readline(fin)

        # For lvect part
        readline(fin)
        for i = 1:3
            latt.lvect[i,:] = parse.(F64, line_to_array(fin))
        end
        readline(fin)

        # For coord part
        readline(fin)
        for i = 1:natom
            latt.coord[i,:] = parse.(F64, line_to_array(fin))
        end
    end # END OF IOSTREAM

    # Print some useful information
    println("  > System: ", latt._case)
    println("  > Atoms: ", latt.atoms)

    # Return the desired struct
    return latt
end

"""
    irio_kmesh(f::String, kmesh::Array{F64,2}, weight::Array{F64,1})

Write the kmesh and weight information to `kmesh.ir` using the IR format.
Here `f` means only the directory that we want to use.

See also: [`vaspio_kmesh`](@ref).
"""
function irio_kmesh(f::String, kmesh::Array{F64,2}, weight::Array{F64,1})
    # Extract some key parameters
    nkpt, ndir = size(kmesh)

    # Extract some key parameters
    _nkpt, = size(weight)

    # Sanity check
    @assert nkpt == _nkpt

    # Output the data
    open(joinpath(f, "kmesh.ir"), "w") do fout
        # Write the header
        println(fout, "# File: kmesh.ir")
        println(fout, "# Data: kmesh[nkpt,ndir] and weight[nkpt]")
        println(fout)
        println(fout, "nkpt -> $nkpt")
        println(fout, "ndir -> $ndir")
        println(fout)

        # Write the body
        for k = 1:nkpt
            @printf(fout, "%16.12f %16.12f %16.12f", kmesh[k, 1:3]...)
            @printf(fout, "%8.2f\n", weight[k])
        end
    end # END OF IOSTREAM

    # Print some useful information
    println("  > Open and write the file kmesh.ir (kmesh and weight)")
end

"""
    irio_kmesh(f::String)

Extract the kmesh and weight information from `kmesh.ir`. Here `f` means
the directory that this file exists.
"""
function irio_kmesh(f::String)
    # Print the header
    println("Open and parse the file kmesh.ir (kmesh and weight)")

    # Check file's status
    fn = joinpath(f, "kmesh.ir")
    @assert isfile(fn)

    # Define kmesh and weight
    kmesh = nothing
    weight = nothing

    # Input the data
    open(fn, "r") do fin
        # Skip the header
        readline(fin)
        readline(fin)
        readline(fin)

        # Extract some dimensional parameters
        nkpt = parse(I64, line_to_array(fin)[3])
        ndir = parse(I64, line_to_array(fin)[3])
        @assert ndir == 3
        readline(fin)

        # Allocate memories
        kmesh = zeros(F64, nkpt, ndir)
        weight = zeros(F64, nkpt)

        # Read the body
        for k = 1:nkpt
            line = line_to_array(fin)
            kmesh[k,:] = parse.(F64, line[1:3])
            weight[k] = parse(F64, line[4])
        end
    end # END OF IOSTREAM

    # Print some useful information
    println("  > Number of k-points: ", length(weight))
    println("  > Total sum of weights: ", sum(weight))
    println("  > Shape of Array kmesh: ", size(kmesh))
    println("  > Shape of Array weight: ", size(weight))

    # Return the desired arrays
    return kmesh, weight
end

"""
    irio_tetra(f::String, volt::F64, itet::Array{I64,2})

Write the tetrahedra information to `tetra.ir` using the IR format. Here
`f` means only the directory that we want to use.

See also: [`vaspio_tetra`](@ref).
"""
function irio_tetra(f::String, volt::F64, itet::Array{I64,2})
    # Extract some key parameters
    ntet, ndim = size(itet)

    # Sanity check
    @assert ndim == 5

    # Output the data
    open(joinpath(f, "tetra.ir"), "w") do fout
        # Write the header
        println(fout, "# File: tetra.ir")
        println(fout, "# Data: itet[ntet,5]")
        println(fout)
        println(fout, "ntet -> $ntet")
        println(fout, "volt -> $volt")
        println(fout)

        # Write the body
        for t = 1:ntet
            @printf(fout, "%8i %8i %8i %8i %8i\n", itet[t, :]...)
        end
    end # END OF IOSTREAM

    # Print some useful information
    println("  > Open and write the file tetra.ir (itet and volt)")
end

"""
    irio_tetra(f::String)

Extract the tetrahedra information from `tetra.ir`. Here `f` means the
directory that this file exists.
"""
function irio_tetra(f::String)
    # Print the header
    println("Open and parse the file tetra.ir (itet and volt)")

    # Check file's status
    fn = joinpath(f, "tetra.ir")
    @assert isfile(fn)

    # Define tetrahedra
    volt = 0.0
    itet = nothing

    # Input the data
    open(fn, "r") do fin
        # Skip the header
        readline(fin)
        readline(fin)
        readline(fin)

        # Extract some dimensional parameters
        ntet = parse(I64, line_to_array(fin)[3])
        volt = parse(F64, line_to_array(fin)[3])
        readline(fin)

        # Allocate memory
        itet = zeros(I64, ntet, 5)

        # Read the body
        for t = 1:ntet
            itet[t, :] = parse.(I64, line_to_array(fin))
        end
    end # END OF IOSTREAM

    # Print some useful information
    println("  > Number of tetrahedra: ", size(itet)[1])
    println("  > Volume of one tetrahedron: ", volt)
    println("  > Shape of Array itet: ", size(itet))

    # Return the desired arrays
    return volt, itet
end

"""
    irio_eigen(f::String, enk::Array{F64,3}, occupy::Array{F64,3})

Write the eigenvalues to `eigen.ir` using the IR format. Here `f` means
only the directory that we want to use.

See also: [`vaspio_eigen`](@ref), [`qeio_eigen`](@ref).
"""
function irio_eigen(f::String, enk::Array{F64,3}, occupy::Array{F64,3})
    # Extract some key parameters
    nband, nkpt, nspin = size(enk)

    # Extract some key parameters
    _nband, _nkpt, _nspin = size(occupy)

    # Sanity check
    @assert nband == _nband && nkpt == _nkpt && nspin == _nspin

    # Output the data
    open(joinpath(f, "eigen.ir"), "w") do fout
        # Write the header
        println(fout, "# File: eigen.ir")
        println(fout, "# Data: enk[nband,nkpt,nspin] and occupy[nband,nkpt,nspin]")
        println(fout)
        println(fout, "nband -> $nband")
        println(fout, "nkpt  -> $nkpt ")
        println(fout, "nspin -> $nspin")
        println(fout)

        # Write the body
        for s = 1:nspin
            for k = 1:nkpt
                for b = 1:nband
                    @printf(fout, "%16.12f %16.12f\n", enk[b, k, s], occupy[b, k, s])
                end # END OF B LOOP
            end # END OF K LOOP
        end # END OF S LOOP
    end # END OF IOSTREAM

    # Print some useful information
    println("  > Open and write the file eigen.ir (enk and occupy)")
end

"""
    irio_eigen(f::String)

Extract the eigenvalues from `eigen.ir`. Here `f` means the directory
that this file exists.
"""
function irio_eigen(f::String)
    # Print the header
    println("Open and parse the file eigen.ir (enk and occupy)")

    # Check file's status
    fn = joinpath(f, "eigen.ir")
    @assert isfile(fn)

    # Define eigenvalues and occupations
    enk = nothing
    occupy = nothing

    # Input the data
    open(fn, "r") do fin
        # Skip the header
        readline(fin)
        readline(fin)
        readline(fin)

        # Extract some dimensional parameters
        nband = parse(I64, line_to_array(fin)[3])
        nkpt  = parse(I64, line_to_array(fin)[3])
        nspin = parse(I64, line_to_array(fin)[3])
        readline(fin)

        # Allocate memories
        enk = zeros(F64, nband, nkpt, nspin)
        occupy = zeros(F64, nband, nkpt, nspin)

        # Read the body
        for s = 1:nspin
            for k = 1:nkpt
                for b = 1:nband
                    line = line_to_array(fin)
                    enk[b, k, s] = parse(F64, line[1])
                    occupy[b, k, s] = parse(F64, line[2])
                end # END OF B LOOP
            end # END OF K LOOP
        end # END OF S LOOP
    end # END OF IOSTREAM

    # Print some useful information
    println("  > Number of DFT bands: ", size(enk)[1])
    println("  > Number of k-points: ", size(enk)[2])
    println("  > Number of spins: ", size(enk)[3])
    println("  > Shape of Array enk: ", size(enk))
    println("  > Shape of Array occupy: ", size(occupy))

    # Return the desired arrays
    return enk, occupy
end

"""
    irio_projs(f::String, chipsi::Array{Array{C64,4},1})

Write the projectors to `projs.ir` using the IR format. Here `f` means
only the directory that we want to use.

The projectors have been processed to fulfill the requirement of the
DMFT engine.

See also: [`vaspio_projs`](@ref).
"""
function irio_projs(f::String, chipsi::Array{Array{C64,4},1})
    # Output the data
    open(joinpath(f, "projs.ir"), "w") do fout
        # Write the header
        println(fout, "# File: projs.ir")
        println(fout, "# Data: chipsi[nproj,nband,nkpt,nspin]")
        println(fout)

        # Go through each PrGroup / PrWindow
        for p in eachindex(chipsi)
            # Extract some key parameters
            ndim, nbnd, nkpt, nspin = size(chipsi[p])

            # Write the header
            println(fout, "group -> $p")
            println(fout, "nproj -> $ndim")
            println(fout, "nband -> $nbnd")
            println(fout, "nkpt  -> $nkpt ")
            println(fout, "nspin -> $nspin")
            println(fout)

            # Write the body
            for s = 1:nspin
                for k = 1:nkpt
                    for b = 1:nbnd
                        for d = 1:ndim
                            z = chipsi[p][d, b, k, s]
                            @printf(fout, "%16.12f %16.12f\n", real(z), imag(z))
                        end # END OF D LOOP
                    end # END OF B LOOP
                end # END OF K LOOP
            end # END OF S LOOP
            println(fout)
        end # END OF P LOOP
    end # END OF IOSTREAM

    # Print some useful information
    println("  > Open and write the file projs.ir (chipsi)")
end

"""
    irio_projs(f::String)

Extract the normalized projectors from `projs.ir`. Here `f` means the
directory that this file exists.
"""
function irio_projs(f::String)
    # Print the header
    println("Open and parse the file projs.ir (chipsi)")

    # Check file's status
    fn = joinpath(f, "projs.ir")
    @assert isfile(fn)

    # Count how many groups there are
    lines = readlines(fn)
    filter!(x -> contains(x, "group"), lines)
    ngroup = length(lines)
    @assert ngroup ≥ 1

    # Define the projectors. They will be filled later.
    chipsi = Array{C64,4}[]

    # Input the data
    open(fn, "r") do fin
        # Skip the header
        readline(fin)
        readline(fin)
        readline(fin)

        # Go through each PrGroup / PrWindow
        for g = 1:ngroup
            # Extract some dimensional parameters
            _g    = parse(I64, line_to_array(fin)[3])
            ndim  = parse(I64, line_to_array(fin)[3])
            nbnd  = parse(I64, line_to_array(fin)[3])
            nkpt  = parse(I64, line_to_array(fin)[3])
            nspin = parse(I64, line_to_array(fin)[3])
            readline(fin)

            # Allocate memory
            P = zeros(C64, ndim, nbnd, nkpt, nspin)

            # Read the body
            for s = 1:nspin
                for k = 1:nkpt
                    for b = 1:nbnd
                        for d = 1:ndim
                            line = line_to_array(fin)
                            _re = parse(F64, line[1])
                            _im = parse(F64, line[2])
                            P[d, b, k, s] = _re + _im * im
                        end # END OF D LOOP
                    end # END OF B LOOP
                end # END OF K LOOP
            end # END OF S LOOP
            readline(fin)

            # Store P in chipsi
            push!(chipsi, P)

            # Print some useful information
            println("  > Shape of Array chipsi: ", size(P))
        end # END OF G LOOP
    end # END OF IOSTREAM

    # Return the desired arrays
    return chipsi
end

"""
    irio_fermi(f::String, fermi::F64)

Write the fermi level to `fermi.ir` using the IR format. Here `f` means
only the directory that we want to use.

See also: [`vaspio_fermi`](@ref), [`qeio_fermi`](@ref).
"""
function irio_fermi(f::String, fermi::F64)
    # Output the data
    open(joinpath(f, "fermi.ir"), "w") do fout
        # Write the header
        println(fout, "# File: fermi.ir")
        println(fout, "# Data: fermi")
        println(fout)
        println(fout, "fermi -> $fermi")
        println(fout)

        # Write the body
        # N/A
    end # END OF IOSTREAM

    # Print some useful information
    println("  > Open and write the file fermi.ir (fermi)")
end

"""
    irio_fermi(f::String)

Extract the fermi level from `fermi.ir`. Here `f` means the directory
that this file exists.
"""
function irio_fermi(f::String)
    # Print the header
    println("Open and parse the file fermi.ir (fermi)")

    # Check file's status
    fn = joinpath(f, "fermi.ir")
    @assert isfile(fn)

    # Define the fermi level
    fermi = 0.0

    # Input the data
    open(fn, "r") do fin
        # Skip the header
        readline(fin)
        readline(fin)
        readline(fin)

        # Extract the fermi level
        fermi = parse(F64, line_to_array(fin)[3])
    end # END OF IOSTREAM

    # Print some useful information
    println("  > Fermi level: $fermi eV")

    # Return the desired value
    return fermi
end

"""
    irio_rawcp(f::String, chipsi::Array{C64,4})

Write the projectors to `rawcp.ir` using the IR format. Here `f` means
only the directory that we want to use.

The projectors are original data. They have not been modified.

See also: [`irio_projs`](@ref).
"""
function irio_rawcp(f::String, chipsi::Array{C64,4})
    # Extract some key parameters
    nproj, nband, nkpt, nspin = size(chipsi)

    # Output the data
    open(joinpath(f, "rawcp.ir"), "w") do fout
        # Write the header
        println(fout, "# File: rawcp.ir")
        println(fout, "# Data: chipsi[nproj,nband,nkpt,nspin]")
        println(fout)
        println(fout, "nproj -> $nproj")
        println(fout, "nband -> $nband")
        println(fout, "nkpt  -> $nkpt ")
        println(fout, "nspin -> $nspin")
        println(fout)

        # Write the body
        for s = 1:nspin
            for k = 1:nkpt
                for b = 1:nband
                    for p = 1:nproj
                        z = chipsi[p, b, k, s]
                        @printf(fout, "%16.12f %16.12f\n", real(z), imag(z))
                    end # END OF P LOOP
                end # END OF B LOOP
            end # END OF K LOOP
        end # END OF S LOOP
    end # END OF IOSTREAM

    # Print some useful information
    println("  > Open and write the file rawcp.ir (chipsi)")
end

"""
    irio_rawcp(f::String)

Extract the raw projectors from `rawcp.ir`. Here `f` means the directory
that this file exists.
"""
function irio_rawcp(f::String)
    # Print the header
    println("Open and parse the file rawcp.ir (chipsi)")

    # Check file's status
    fn = joinpath(f, "rawcp.ir")
    @assert isfile(fn)

    # Define the projectors
    chipsi = nothing

    # Input the data
    open(fn, "r") do fin
        # Skip the header
        readline(fin)
        readline(fin)
        readline(fin)

        # Extract some dimensional parameters
        nproj = parse(I64, line_to_array(fin)[3])
        nband = parse(I64, line_to_array(fin)[3])
        nkpt  = parse(I64, line_to_array(fin)[3])
        nspin = parse(I64, line_to_array(fin)[3])
        readline(fin)

        # Allocate memory
        chipsi = zeros(C64, nproj, nband, nkpt, nspin)

        # Read the body
        for s = 1:nspin
            for k = 1:nkpt
                for b = 1:nband
                    for p = 1:nproj
                        _re, _im = parse.(F64, line_to_array(fin))
                        chipsi[p, b, k, s] = _re + _im * im
                    end # END OF P LOOP
                end # END OF B LOOP
            end # END OF K LOOP
        end # END OF S LOOP
    end # END OF IOSTREAM

    # Print some useful information
    println("  > Shape of Array chipsi: ", size(chipsi))

    # Return the desired arrays
    return chipsi
end

"""
    irio_charge(f::String)

Write the charge density to `charge.ir` using the IR format. Here `f`
means only the directory that we want to use.

See also: [`vaspio_charge`](@ref).
"""
function irio_charge(f::String) end
