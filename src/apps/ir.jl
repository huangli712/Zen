function ortho() end

function density_matrix() end

"""
    irio_lattice()

Write the lattice information using the IR format
"""
function irio_lattice() end

"""
    irio_kmesh()

Write the k-mesh information using the IR format 
"""
function irio_kmesh(f::AbstractString, kmesh::Array{F64,2}, weight::Array{F64,1})
    # extract some key parameters
    nkpt, ndir = size(kmesh)

    # output the data
    open(f * "/kmesh.ir", "w") do fout
        for k = 1:nkpt
            println(fout, kmesh[k, :], weight[k])
        end
    end
end

"""
    irio_tetra(f::AbstractString, ntet::I64, volt::F64, itet::Array{I64,2})

Write the tetrahedra information using the IR format
"""
function irio_tetra(f::AbstractString, ntet::I64, volt::F64, itet::Array{I64,2})
    open(f * "/tetra.ir", "w") do fout
        @assert (ntet, 5) === size(itet)
        println(fout, "ntet: $ntet")
        println(fout, "volt: $volt")
        for t = 1:ntet
            println(fout, itet[t, :])
        end
    end
end

"""
    irio_projs()

Write the projectors using the IR format
"""
function irio_projs() end

"""
    irio_eigen(f::AbstractString, enk::Array{F64,3}, occupy::Array{F64,3})

Write the eigenvalues using the IR format
"""
function irio_eigen(f::AbstractString, enk::Array{F64,3}, occupy::Array{F64,3})
    # extract some key parameters
    nkpt, nband, nspin = size(enk)

    # output the data
    open(f * "/eigen.ir", "w") do fout
        println(fout, "nkpt : $nkpt ")
        println(fout, "nband: $nband")
        println(fout, "nspin: $nspin")
        for s = 1:nspin
            for b = 1:nband
                for k = 1:nkpt
                    println(fout, enk[k, b, s], " ", occupy[k, b, s])
                end
            end
        end
    end
end
