module Adaptor

include("vasp.jl")
export from_poscar
export from_ibzkpt
export from_projcar
export from_locproj
export from_eigenval
export from_chgcar
export to_chgcar

end
function from_poscar(f::AbstractString)
end

function from_projcar(f::AbstractString)
end

function from_locproj(f::AbstractString)
end

"""
    from_ibzkpt

Reading IBZKPT, return k-mesh and k-weight
"""
function from_ibzkpt(f::AbstractString)
    fin = open(f * "/IBZKPT", "r")

    readline(fin)
    str = readline(fin)
    nkpt = parse(Int64, str)
    readline(fin)

    kmesh = zeros(Float64, nkpt, 3)
    weight = zeros(Float64, nkpt) 

    i = 0
    for line in eachline(fin)
        substr = split(line, " ", keepempty = false)
        subarr = map(x -> parse(Float64, x),  substr[1:3]) 
        i = i + 1
        kmesh[i,1:3] = subarr
        weight[i] = parse(Int64, substr[4])
        println(kmesh[i,:], " ", weight[i])
    end
    @assert( i === nkpt )
    println(sum(weight))

    close(fin)
end

function from_eigenval(f::AbstractString)
end

function from_chgcar(f::AbstractString)
end

function to_chgcar()
end
