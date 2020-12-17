#
# project : pansy
# source  : ir.jl
# author  : Li Huang (lihuang.dmft@gmail.com)
# status  : unstable
# comment :
#
# last modified: 2020/12/17
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
    nkpt, nband, nspin = size(enk)

    # output the data
    open(joinpath(f, "eigen.ir"), "w") do fout
        println(fout, "# file: eigen.ir")
        println(fout, "# data: enk[nkpt,nband,nspin] and occupy[nkpt,nband,nspin]")
        println(fout)
        println(fout, "nkpt  -> $nkpt ")
        println(fout, "nband -> $nband")
        println(fout, "nspin -> $nspin")
        println(fout)
        for s = 1:nspin
            for b = 1:nband
                for k = 1:nkpt
                    @printf(fout, "%16.12f %16.12f\n", enk[k, b, s], occupy[k, b, s])
                end
            end
        end
    end
end

"""
    irio_projs()

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
    irio_fermi()
"""
function irio_fermi() end

"""
    irio_charge()
"""
function irio_charge() end
