#
# project : pansy
# source  : ir.jl
# author  : Li Huang (lihuang.dmft@gmail.com)
# status  : unstable
# comment :
#
# last modified: 2021/02/10
#

#
# Driver Functions
#

"""
    ir_adaptor(D::Dict{Symbol,Any})

Write the Kohn-Sham data to specified files using the IR format. Note
that the Kohn-Sham data are encapsulated in the `D` dict.

See also: [`vasp_adaptor`](@ref), [`plo_adaptor`](@ref).
"""
function ir_adaptor(D::Dict{Symbol,Any})
    # I01: Print the header
    println("  < IR Adaptor >")

    # I02: Check the validity of the `D` dict
    key_list = [:latt, :kmesh, :weight, :enk, :occupy, :chipsi_f, :fermi]
    for k in key_list
        @assert haskey(D, k)
    end

    # I03: Write important parameters
    println("    Put Params")
    irio_params(pwd(), D)

    # I04: Write lattice structure
    println("    Put Lattice")
    irio_lattice(pwd(), D[:latt])

    # I05: Write kmesh and the corresponding weights
    println("    Put Kmesh")
    println("    Put Weight")
    irio_kmesh(pwd(), D[:kmesh], D[:weight])

    # I06: Write band structure and the corresponding occupancies
    println("    Put Enk")
    println("    Put Occupy")
    irio_eigen(pwd(), D[:enk], D[:occupy])

    # I07: Write projectors, traits, and groups
    println("    Put Projector")
    irio_projs(pwd(), D[:chipsi_f])

    # I08: Write fermi level
    println("    Put Fermi Level")
    irio_fermi(pwd(), D[:fermi])

#
# Remarks:
#
# The following steps are optional, because the tetrahedron information
# might be absent.
#

    # I09: Check the validity of the `D` dict further (optional)
    if get_d("smear") === "tetra"
        key_list = [:volt, :itet]
        for k in key_list
            @assert haskey(D, k)
        end
    end

    # I10: Write tetrahedron data if they are available
    if get_d("smear") === "tetra"
        println("    Put Tetrahedron")
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
    file_list = ["params", "lattice", "kmesh", "eigen", "projs", "fermi"]
    if get_d("smear") === "tetra"
        push!(file_list, "tetra")
    end

    # Store the data files
    for i in eachindex(file_list)
        file_src = file_list[i] * ".ir"
        file_dst = file_list[i] * ".ir.$(it.dmft_cycle)"
        cp(file_src, file_dst, force = true)
    end
end

#
# Service Functions
#

"""
    irio_params(f::String, D::Dict{Symbol,Any})

Write the key parameters extracted from the Kohn-Sham data. Here `f`
means only the directory that we want to use.

See also: [`PrGroup`](@ref), [`PrWindow`](@ref).
"""
function irio_params(f::String, D::Dict{Symbol,Any})
    # Extract some key parameters
    nband, nkpt, nspin = size(D[:enk])

    # Extract `ntet`, it is optional.
    if haskey(D, :itet)
        ntet, _ = size(D[:itet])
    else
        ntet = 0
    end

    # Extract `volt`, it is optional.
    if haskey(D, :volt)
        volt = D[:volt]
    else
        volt = 0.0
    end

    # Extract `ngrp`
    ngrp, = size(D[:PG])

    # Extract `nwnd`
    nwnd, = size(D[:PW])
    @assert ngrp === nwnd

    # Output the data
    open(joinpath(f, "params.ir"), "w") do fout
        # Write the header
        println(fout, "# file: traits.ir")
        println(fout, "# data: some necessary parameters")
        println(fout)

        # Write basic parameters
        println(fout, "# Common  :")
        println(fout, "nband -> $nband")
        println(fout, "nkpt  -> $nkpt")
        println(fout, "nspin -> $nspin")
        println(fout, "ntet  -> $ntet")
        println(fout, "volt  -> $volt")
        println(fout, "ngrp  -> $ngrp")
        println(fout, "nwnd  -> $nwnd")
        println(fout)

        # Write PrGroup[]
        for p in eachindex(D[:PG])
            PG = D[:PG]
            println(fout, "# PrGroup : $p")
            println(fout, "site  -> $(PG[p].site)")
            println(fout, "l     -> $(PG[p].l)")
            println(fout, "corr  -> $(PG[p].corr)")
            println(fout, "shell -> $(PG[p].shell)")
            println(fout, "ndim  -> $(size(PG[p].Tr,1))")
        end
        println(fout)

        # Write PrWindow[]
        for p in eachindex(D[:PW])
            PW = D[:PW]
            println(fout, "# PrWindow: $p")
            println(fout, "bmin  -> $(PW[p].bmin)")
            println(fout, "bmax  -> $(PW[p].bmax)")
            println(fout, "nbnd  -> $(PW[p].nbnd)")
        end
        println(fout)
    end
end

"""
    irio_lattice(f::String, latt::Lattice)

Write the lattice information to lattice.ir using the IR format. Here `f`
means only the directory that we want to use.

See also: [`vaspio_lattice`](@ref).
"""
function irio_lattice(f::String, latt::Lattice)
    # Extract some key parameters
    _case, scale, nsort, natom = latt._case, latt.scale, latt.nsort, latt.natom

    # Output the data
    open(joinpath(f, "lattice.ir"), "w") do fout
        # Write the header
        println(fout, "# file: lattice.ir")
        println(fout, "# data: Lattice struct")
        println(fout)
        println(fout, "scale -> $_case")
        println(fout, "scale -> $scale")
        println(fout, "nsort -> $nsort")
        println(fout, "natom -> $natom")
        println(fout)

        # Write the body
        # For sorts part
        println(fout, "[sorts]")
        for i = 1:nsort # Symbols
            @printf(fout, "%6s", latt.sorts[i, 1])
        end
        println(fout)
        for i = 1:nsort # Numbers
            @printf(fout, "%6i", latt.sorts[i, 2])
        end
        println(fout)
        println(fout)

        # For atoms part
        println(fout, "[atoms]")
        for i = 1:natom
            @printf(fout, "%6s", latt.atoms[i])
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
    end
end

"""
    irio_kmesh(f::String, kmesh::Array{F64,2}, weight::Array{F64,1})

Write the kmesh and weight information to kmesh.ir using the IR format. Here
`f` means only the directory that we want to use.

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
        println(fout, "# file: kmesh.ir")
        println(fout, "# data: kmesh[nkpt,ndir] and weight[nkpt]")
        println(fout)
        println(fout, "nkpt -> $nkpt")
        println(fout, "ndir -> $ndir")
        println(fout)

        # Write the body
        for k = 1:nkpt
            @printf(fout, "%16.12f %16.12f %16.12f", kmesh[k, 1:3]...)
            @printf(fout, "%8.2f\n", weight[k])
        end
    end
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
        println(fout, "# file: tetra.ir")
        println(fout, "# data: itet[ntet,5]")
        println(fout)
        println(fout, "ntet -> $ntet")
        println(fout, "volt -> $volt")
        println(fout)

        # Write the body
        for t = 1:ntet
            @printf(fout, "%8i %8i %8i %8i %8i\n", itet[t, :]...)
        end
    end
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
        println(fout, "# file: eigen.ir")
        println(fout, "# data: enk[nband,nkpt,nspin] and occupy[nband,nkpt,nspin]")
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
                end
            end
        end
    end
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
        println(fout, "# file: projs.ir")
        println(fout, "# data: chipsi[nproj,nband,nkpt,nspin]")
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
                    end
                end
            end
        end
    end
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
        println(fout, "# file: projs.ir")
        println(fout, "# data: chipsi[nproj,nband,nkpt,nspin]")
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
                        end
                    end
                end
            end
            println(fout)
        end
    end
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
        println(fout, "# file: fermi.ir")
        println(fout, "# data: fermi")
        println(fout)
        println(fout, "fermi -> $fermi")
        println(fout)

        # Write the body
        # N/A
    end
end

"""
    irio_charge(f::String)

Write the charge density to charge.ir using the IR format. Here `f` means
only the directory that we want to use.

See also: [`vaspio_charge`](@ref).
"""
function irio_charge(f::String) end
