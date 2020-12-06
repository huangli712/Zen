"""
    from_poscar(f::AbstractString)

Reading vasp's POSCAR file, return crystallography information. Here `f`
means only the directory that contains POSCAR 
"""
function from_poscar(f::AbstractString)
    # open the iostream
    fin = open(f * "/POSCAR", "r")

    # get the case
    case = strip(readline(fin))

    # get the scaling factor
    scale = parse(F64, readline(fin))

    # get the basis vector
    bvec = zeros(F64, 3, 3)
    bvec[1,:] = parse.(F64, split(readline(fin), " ", keepempty = false))
    bvec[2,:] = parse.(F64, split(readline(fin), " ", keepempty = false))
    bvec[3,:] = parse.(F64, split(readline(fin), " ", keepempty = false))

    # get the symbol list
    symbols = split(readline(fin), " ", keepempty = false)
    @show symbols

    # close the iostream
    close(fin)
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
        weight[i] = subarr[4]
    end

    # make sure the number of k-points is correct
    @assert( i === nkpt )

    # close the iostream
    close(fin)

    # return the desired arrays
    return kmesh, weight
end

function from_eigenval(f::AbstractString)
end

function from_chgcar(f::AbstractString)
end

function to_chgcar()
end
