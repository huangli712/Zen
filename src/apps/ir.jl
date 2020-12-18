#
# project : pansy
# source  : ir.jl
# author  : Li Huang (lihuang.dmft@gmail.com)
# status  : unstable
# comment :
#
# last modified: 2020/12/18
#

"""
    irio_lattice()

Write the lattice information using the IR format
"""
function irio_lattice() end

"""
    irio_kmesh(f::AbstractString, kmesh::Array{F64,2}, weight::Array{F64,1})

Write the kmesh information using the IR format
"""
function irio_kmesh(f::AbstractString, kmesh::Array{F64,2}, weight::Array{F64,1})
    # extract some key parameters
    nkpt, ndir = size(kmesh)

    # output the data
    open(joinpath(f, "kmesh.ir"), "w") do fout
        println(fout, "# file: kmesh.ir")
        println(fout, "# data: kmesh[nkpt,ndir] and weight[nkpt]")
        println(fout)
        println(fout, "nkpt -> $nkpt")
        println(fout, "ndir -> $ndir")
        println(fout)
        for k = 1:nkpt
            @printf(fout, "%16.12f %16.12f %16.12f %8.2f\n", kmesh[k, 1:3]..., weight[k])
        end
    end
end

"""
    irio_tetra(f::AbstractString, volt::F64, itet::Array{I64,2})

Write the tetrahedra information using the IR format
"""
function irio_tetra(f::AbstractString, volt::F64, itet::Array{I64,2})
    # extract some key parameters
    ntet, = size(itet)

    # output the data
    open(joinpath(f, "tetra.ir"), "w") do fout
        println(fout, "# file: tetra.ir")
        println(fout, "# data: itet[ntet,5]")
        println(fout)
        println(fout, "ntet -> $ntet")
        println(fout, "volt -> $volt")
        println(fout)
        for t = 1:ntet
            @printf(fout, "%8i %8i %8i %8i %8i\n", itet[t, :]...)
        end
    end
end

"""
    irio_eigen(f::AbstractString, enk::Array{F64,3}, occupy::Array{F64,3})

Write the eigenvalues using the IR format
"""
function irio_eigen(f::AbstractString, enk::Array{F64,3}, occupy::Array{F64,3})
    # extract some key parameters
    nband, nkpt, nspin = size(enk)

    # output the data
    open(joinpath(f, "eigen.ir"), "w") do fout
        println(fout, "# file: eigen.ir")
        println(fout, "# data: enk[nband,nkpt,nspin] and occupy[nband,nkpt,nspin]")
        println(fout)
        println(fout, "nband -> $nband")
        println(fout, "nkpt  -> $nkpt ")
        println(fout, "nspin -> $nspin")
        println(fout)
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
    irio_projs(f::AbstractString, chipsi::Array{C64,4})

Write the projectors using the IR format
"""
function irio_projs(f::AbstractString, chipsi::Array{C64,4})
    # extract some key parameters
    nproj, nband, nkpt, nspin = size(chipsi)

    # output the data
    open(joinpath(f, "projs.ir"), "w") do fout
        println(fout, "# file: projs.ir")
        println(fout, "# data: chipsi[nproj,nband,nkpt,nspin]")
        println(fout)
        println(fout, "nproj -> $nproj")
        println(fout, "nband -> $nband")
        println(fout, "nkpt  -> $nkpt ")
        println(fout, "nspin -> $nspin")
        println(fout)
        for s = 1:nspin
            for k = 1:nkpt
                for b = 1:nband
                    for p = 1:nproj
                        _re = real(chipsi[p,b,k,s])
                        _im = imag(chipsi[p,b,k,s])
                        @printf(fout, "%16.12f %16.12f\n", _re, _im)
                    end
                end
            end
        end
    end
end

"""
    irio_fermi(f::AbstractString, fermi::F64)

Write the fermi level using the IR format
"""
function irio_fermi(f::AbstractString, fermi::F64)
    # output the data
    open(joinpath(f, "fermi.ir"), "w") do fout
        println(fout, "# file: fermi.ir")
        println(fout, "# data: fermi")
        println(fout)
        println(fout, "fermi -> $fermi")
        println(fout)
    end
end

"""
    irio_charge()
"""
function irio_charge() end
