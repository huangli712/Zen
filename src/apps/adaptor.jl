function from_poscar(f::AbstractString)
end

function from_projcar(f::AbstractString)
end

function from_locproj(f::AbstractString)
end

"""
    from_ibzkpt(f::AbstractString)

Reading vasp's IBZKPT file, return k-mesh and k-weight. Here `f` means
only the directory that contains IBZKPT
"""
function from_ibzkpt(f::AbstractString)
    # open the iostream
    fin = open(f * "/IBZKPT", "r")

    # get number of k-points
    readline(fin)
    str = readline(fin)
    nkpt = parse(I64, str)
    readline(fin)

    # create arrays 
    kmesh = zeros(F64, nkpt, 3)
    weight = zeros(F64, nkpt) 

    # read in the k-points and their weights
    i = 0
    for line in eachline(fin)
        substr = split(line, " ", keepempty = false)
        subarr = map(x -> parse(F64, x),  substr) 
        i = i + 1
        kmesh[i,1:3] = subarr[1:3]
        weight[i] = substr[4]
        println(kmesh[i,:], " ", weight[i])
    end

    # make sure the number of k-points is correct
    @assert( i === nkpt )
    println(sum(weight))

    # close the iostream
    close(fin)
end

function from_eigenval(f::AbstractString)
end

function from_chgcar(f::AbstractString)
end

function to_chgcar()
end
