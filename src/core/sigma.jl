#
# Project : Pansy
# Source  : sigma.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/03/26
#

"""
    sigma_reset(lr::Logger)

Create initial self-energy functions and write them to `sigma.bare`.

See also: [`sigma_dcount`](@ref).
"""
function sigma_reset(lr::Logger)
    # Print the log
    prompt(lr.log, "sigma::reset")
    println("< Sigma : Reset >")

    # The sdim creates a mapping from shell (string) to ndim (integer).
    # It is used to parse get_i("shell") to extract the `ndim` parameter.
    sdim = Dict{String,I64}(
               "s"     => 1,
               "p"     => 3,
               "d"     => 5,
               "f"     => 7,
               "d_t2g" => 3, # Only a subset of d orbitals
               "d_eg"  => 2, # Only a subset of d orbitals
           )

    # Extract some necessary parameters
    axis = get_m("axis")
    nmesh = get_m("nmesh")
    beta = get_m("beta")
    nsite = get_i("nsite")
    nspin = 2

    # Create frequency mesh
    fmesh = zeros(F64, nmesh)
    if axis === 1 # Imaginary axis
        for i = 1:nmesh
            fmesh[i] = (2 * i - 1) * pi / beta
        end
    else # Real axis
        sorry()
    end

    # Create self-energy functions
    #
    # Initialize an array for self-energy functions
    SA = Array{C64,3}[]
    D = I64[]
    #
    # Go through the impurity problems
    for i = 1:nsite
        # Retrieve specification for impurity problem
        str = get_i("shell")[i]

        # Get the dimension of impurity problem
        ndim = get(sdim, str, 1)
        push!(D, ndim)

        # Create a temporary array for self-energy function
        S = zeros(C64, nmesh, ndim, nspin)

        # Push S into SA to save it
        push!(SA, S)
    end

    # Write self-energy functions to sigma.bare
    open("sigma.bare", "w") do fout
        # Write the header
        println(fout, "# File: sigma.bare")
        println(fout, "# Data: bare self-energy functions")
        println(fout)
        println(fout, "axis  -> $axis")
        println(fout, "beta  -> $beta")
        println(fout, "nsite -> $nsite") 
        println(fout, "nmesh -> $nmesh")
        println(fout, "nspin -> $nspin")
        for i = 1:nsite
            println(fout, "ndim$i -> $(D[i])")
        end
        println(fout)

        # Write the body
        # Go through each impurity problem
        for i = 1:nsite
            for s = 1:nspin
                println(fout, "# site: $i spin: $s")
                for m = 1:nmesh
                    @printf(fout, "%16.12f", fmesh[m])
                    foreach(x -> @printf(fout, "%16.12f %16.12f", real(x), imag(x)), SA[i][m, :, s])
                    println(fout)
                end
                println(fout)
                println(fout)
            end
        end
    end
end

"""
    sigma_dcount(lr::Logger)

Calculate double counting terms for self-energy functions and write
them to `sigma.dc`.

See also: [`sigma_reset`](@ref).
"""
function sigma_dcount(lr::Logger)
    # Print the log
    prompt(lr.log, "sigma::dcount")
    prompt("Sigma : Dcount")
end

"""
    sigma_split(lr::Logger)

Split the hybridization functions (or similar local functions) and then
distribute them into the `impurity.i` folder.

See also: [`sigma_gather`](@ref).
"""
function sigma_split(lr::Logger)
    # Print the log
    prompt(lr.log, "sigma::split")
    prompt("Sigma : Split")
end

"""
    sigma_gather(lr::Logger)

Gather the self-energy functions (or similar local functions) from the
`impurity.i` folder and then combine them into a single `sigma.bare` file.

See also: [`sigma_split`](@ref).
"""
function sigma_gather(lr::Logger)
    # Print the log
    prompt(lr.log, "sigma::gather")
    prompt("Sigma : Gather")
end
