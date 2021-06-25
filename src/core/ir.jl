#
# Project : Pansy
# Source  : ir.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/06/25
#

#=
### *Driver Functions*
=#

#=
*Remarks*:

In the `ir_adaptor()` funtion, some writting steps are optional because
the tetrahedron information might be absent.
=#

"""
    ir_adaptor(D::Dict{Symbol,Any})

Write the Kohn-Sham dataset to specified files using the IR format. Note
that the Kohn-Sham dataset are encapsulated in the `D` dict.

See also: [`vasp_adaptor`](@ref), [`plo_adaptor`](@ref).
"""
function ir_adaptor(D::Dict{Symbol,Any})
    # I01: Check the validity of the `D` dict
    key_list = [:MAP, :PG, :PW, :latt, :kmesh, :weight, :enk, :occupy, :Fchipsi, :fermi]
    for k in key_list
        @assert haskey(D, k)
    end

    # I02: Print the header
    println("Adaptor : IR")
    println("Try to write the processed Kohn-Sham dataset with IR format")
    println("Current directory: ", pwd())

    # I03: Write important parameters
    irio_params(pwd(), D)
    #
    irio_maps(pwd(), D[:MAP])
    #
    irio_groups(pwd(), D[:PG])
    #
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

    # I09: Check the validity of the `D` dict further (optional)
    if get_d("smear") === "tetra"
        key_list = [:volt, :itet]
        for k in key_list
            @assert haskey(D, k)
        end
    end

    # I10: Write tetrahedron data if they are available
    if get_d("smear") === "tetra"
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
    if get_d("smear") === "tetra"
        push!(file_list, "tetra")
    end

    # Store the data files
    for i in eachindex(file_list)
        file_src = file_list[i] * ".ir"
        file_dst = file_list[i] * ".ir.$(it.I₃)"
        cp(file_src, file_dst, force = true)
    end
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
    # Print the header
    println("Store essential parameters")

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
    # Extract `qdim`, maximum number of projectors in all groups.
    # `size(D[:PG][g].Tr,1)` is actually ndim
    qdim = maximum( [ size(D[:PG][g].Tr,1) for g = 1:ngrp ] )

    # Extract `nwnd`
    nwnd, = size(D[:PW])
    #
    # Extract qbnd, maximum number of bands in all windows.
    qbnd = maximum( [ D[:PW][w].nbnd for w = 1:nwnd ] )

    # D[:PW] and D[:PG] should have the same size
    @assert ngrp === nwnd

    # Extract `nsite` and `nmesh`
    nmesh = get_m("nmesh")
    nsite = get_i("nsite")

    # To make sure the validaty of nsite
    @assert nsite <= ngrp

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

Write the information contained in `Mapping`. Here `f` means only the
directory that we want to use.

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
    irio_groups(f::String, PG::Array{PrGroup,1})

Write the information contained in `PrGroup`. Here `f` means only the
directory that we want to use.

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
        println(fout, "ngrp  -> $ngrp")
        println(fout)
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
    irio_windows(f::String, PW::Array{PrWindow,1})

Write the information contained in `PrWindow`. Here `f` means only the
directory that we want to use.

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
        println(fout, "nwnd  -> $nwnd")
        println(fout)
        for p in eachindex(PW)
            println(fout, "# PrWindow: $p")
            println(fout, "bmin  -> $(PW[p].bmin)")
            println(fout, "bmax  -> $(PW[p].bmax)")
            println(fout, "nbnd  -> $(PW[p].nbnd)")
            println(fout, "kwin  ->")
            nkpt, nspin, ndir = size(PW[p].kwin)
            @assert ndir === 2 # For lower and upper boundaries
            for s = 1:nspin
                for k = 1:nkpt
                    @printf(fout, "%8i %8i %8i %8i\n", k, s, PW[p].kwin[k, s, :]...)
                end
            end
            println(fout)
        end # END OF P LOOP
    end # END OF IOSTREAM

    # Print some useful information
    println("  > Open and write the file windows.ir")
end

"""
    irio_lattice(f::String, latt::Lattice)

Write the lattice information to lattice.ir using the IR format. Here `f`
means only the directory that we want to use.

See also: [`vaspio_lattice`](@ref).
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
    irio_kmesh(f::String, kmesh::Array{F64,2}, weight::Array{F64,1})

Write the kmesh and weight information to kmesh.ir using the IR format.
Here `f` means only the directory that we want to use.

See also: [`vaspio_kmesh`](@ref).
"""
function irio_kmesh(f::String, kmesh::Array{F64,2}, weight::Array{F64,1})
    # Extract some key parameters
    nkpt, ndir = size(kmesh)

    # Extract some key parameters
    _nkpt, = size(weight)

    # Sanity check
    @assert nkpt === _nkpt

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
    irio_tetra(f::String, volt::F64, itet::Array{I64,2})

Write the tetrahedra information to tetra.ir using the IR format. Here `f`
means only the directory that we want to use.

See also: [`vaspio_tetra`](@ref).
"""
function irio_tetra(f::String, volt::F64, itet::Array{I64,2})
    # Extract some key parameters
    ntet, ndim = size(itet)

    # Sanity check
    @assert ndim === 5

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
    irio_eigen(f::String, enk::Array{F64,3}, occupy::Array{F64,3})

Write the eigenvalues to eigen.ir using the IR format. Here `f` means
only the directory that we want to use.

See also: [`vaspio_eigen`](@ref).
"""
function irio_eigen(f::String, enk::Array{F64,3}, occupy::Array{F64,3})
    # Extract some key parameters
    nband, nkpt, nspin = size(enk)

    # Extract some key parameters
    _nband, _nkpt, _nspin = size(enk)

    # Sanity check
    @assert nband === _nband && nkpt === _nkpt && nspin === _nspin

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
    irio_projs(f::String, chipsi::Array{C64,4})

Write the projectors to projs.ir using the IR format. Here `f` means
only the directory that we want to use.

The projectors are original data. They have not been modified.

See also: [`vaspio_projs`](@ref).
"""
function irio_projs(f::String, chipsi::Array{C64,4})
    # Extract some key parameters
    nproj, nband, nkpt, nspin = size(chipsi)

    # Output the data
    open(joinpath(f, "projs.ir"), "w") do fout
        # Write the header
        println(fout, "# File: projs.ir")
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
    println("  > Open and write the file projs.ir (chipsi)")
end

"""
    irio_projs(f::String, chipsi::Array{Array{C64,4},1})

Write the projectors to projs.ir using the IR format. Here `f` means
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
    irio_fermi(f::String, fermi::F64)

Write the fermi level to fermi.ir using the IR format. Here `f` means
only the directory that we want to use.

See also: [`vaspio_fermi`](@ref).
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
    irio_charge(f::String)

Write the charge density to charge.ir using the IR format. Here `f` means
only the directory that we want to use.

See also: [`vaspio_charge`](@ref).
"""
function irio_charge(f::String) end
