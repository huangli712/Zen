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

    # get the number list 
    numbers = parse.(I64, split(readline(fin), " ", keepempty = false))

    # get the total number of atoms
    total_atoms = sum(numbers)

    # create atom lists
    atom_list = zeros(I64, total_atoms)
    curr_index = 0
    k = 0
    for i in 1:length(numbers)
        curr_index = curr_index + 1
        for j in 1:numbers[i]
            k = k + 1
            atom_list[k] = curr_index
        end
    end

    # get the coordinates of atoms
    posi_list = zeros(F64, total_atoms, 3)
    readline(fin)
    for i in 1:total_atoms
        posi_list[i,:] = parse.(F64, split(readline(fin), " ", keepempty = false)[1:3])
    end

    # close the iostream
    close(fin)

    return symbols, atom_list, posi_list
end

"""
    from_projcar(f::AbstractString)

Reading vasp's PROJCAR file, return raw projector matrix. Here `f` means
only the directory that contains PROJCAR
"""
function from_projcar(f::AbstractString)
    # get key parameters from the LOCPROJ file
    nspin, nkpt, nband, nproj, nsite = from_locproj(f, true)

    # open the iostream
    fin = open(f * "/PROJCAR", "r")

    # create arrays
    chipsi = zeros(C64, nproj, nband, nkpt, nspin)

    # read in raw projector data
    for site in 1:nsite
        readline(fin)
        readline(fin)
        for spin in 1:nspin
            for kpt in 1:nkpt
                arr = split(readline(fin), " ", keepempty = false)
                curr_kpt = parse(I64, arr[2])
                curr_spin = parse(I64, arr[4])
                @assert curr_kpt === kpt
                @assert curr_spin === spin
                @show arr, curr_kpt, curr_spin
                exit(-1)
            end
        end
    end
    

    # close the iostream
    close(fin)
end

"""
    from_locproj(f::AbstractString, read_param_only::Bool)

Reading vasp's LOCPROJ file, return raw projector matrix. Here `f` means
only the directory that contains LOCPROJ
"""
function from_locproj(f::AbstractString, read_param_only::Bool)
    # open the iostream
    fin = open(f * "/LOCPROJ", "r")
    
    if read_param_only
        arr = split(readline(fin), " ", keepempty = false)
        nspin, nkpt, nband, nproj = tuple(map(x -> parse(I64,x), arr[1:4])...)
        sites = zeros(I64, nproj)
        for i in 1:nproj
            sites[i] = parse(I64, split(readline(fin), " ", keepempty = false)[2])
        end
        nsite = length(union(sites))
        return nspin, nkpt, nband, nproj, nsite
    else
        error("Sorry, this feature has not been implemented")
    end

    # close the iostream
    close(fin)
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

"""
    from_eigenval(f::AbstractString)

Reading vasp's EIGENVAL file, return energy band information. Here `f`
means only the directory that contains EIGENVAL
"""
function from_eigenval(f::AbstractString)
    # open the iostream
    fin = open(f * "/EIGENVAL", "r")

    # determine number of spins
    nspin = parse(I64, split(readline(fin), " ", keepempty = false)[end])

    # skip for lines
    for i in 1:4
        readline(fin)
    end

    # read in some key parameters: nelect, nkpt, nbands 
    nelect, nkpt, nband = tuple(parse.(I64, split(readline(fin), " ", keepempty = false))...)

    # create arrays
    enk = zeros(F64, nkpt, nband, nspin)
    occupy = zeros(F64, nkpt, nband, nspin)

    # read in the energy bands and the corresponding occupations
    for i in 1:nkpt
        readline(fin)
        readline(fin)
        for j in 1:nband
            line = split(readline(fin), " ", keepempty = false)
            if nspin === 1
                enk[i,j,1] = parse(F64, line[2])
                occupy[i,j,1] = parse(F64, line[3])
            else # nspin == 2
                enk[i,j,1] = parse(F64, line[2])
                enk[i,j,2] = parse(F64, line[3])
                occupy[i,j,1] = parse(F64, line[4])
                occupy[i,j,2] = parse(F64, line[5])
            end
            @show i, j, enk[i,j,:], occupy[i,j,:]
        end
    end

    # close the iostream
    close(fin)

    # return the desired arrays
    return enk, occupy
end

function from_chgcar(f::AbstractString)
end

function to_chgcar()
end
