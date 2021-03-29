#
# Project : Pansy
# Source  : sigma.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/03/29
#

"""
    sigma_reset()

Create initial self-energy functions and write them to `sigma.bare`.

See also: [`sigma_dcount`](@ref).
"""
function sigma_reset()
    # Print the log
    println("Sigma : Reset")

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
    println("  Create frequency mesh")

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

        # Get the dimension of impurity problem (save it)
        ndim = get(sdim, str, 1)
        push!(D, ndim)

        # Create a temporary array for self-energy function
        S = zeros(C64, nmesh, ndim, nspin)

        # Push S into SA to save it
        push!(SA, S)
    end
    println("  Create self-energy functions")

    # Write self-energy functions to sigma.bare
    open("dmft1/sigma.bare", "w") do fout
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

    # Print blank line for better visualization
    println()
end

"""
    sigma_dcount()

Calculate double counting terms for self-energy functions and write
them to `sigma.dc`.

See also: [`sigma_reset`](@ref).
"""
function sigma_dcount()
    # Print the log
    println("[ Sigma : Dcount ]")

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
    nsite = get_i("nsite")
    nspin = 2

    # Create double counting terms for self-energy functions
    #
    # Initialize an array for dc
    sigdc = zeros(F64, nsite)
    #
    # Go through the impurity problems
    for i = 1:nsite
        # Get interaction parameters
        U = get_i("upara")[i]
        J = get_i("jpara")[i]

        # Get occupation numbers
        N = get_i("occup")[i]

        # Choose suitable double counting scheme
        @cswitch get_m("dcount") begin
            # Fully localized limit scheme (occupation number is fixed)
            @case "fll1"
                sigdc[i] = cal_dc_fll(U, J, N)
                break

            # Fully localized limit scheme (occupation number is adjusted)
            @case "fll2"
                sorry()
                break

            # Around mean-field scheme
            @case "amf"
                sorry()
                break

            # Exact double counting scheme
            @case "exact"
                sorry()
                break
        end
    end

end

"""
    sigma_split()

Split the hybridization functions (or similar local functions) and then
distribute them into the `impurity.i` folder.

See also: [`sigma_gather`](@ref).
"""
function sigma_split()
    # Print the log
    println("[ Sigma : Split ]")
end

"""
    sigma_gather()

Gather the self-energy functions (or similar local functions) from the
`impurity.i` folder and then combine them into a single `sigma.bare` file.

See also: [`sigma_split`](@ref).
"""
function sigma_gather()
    # Print the log
    println("[ Sigma : Gather ]")
end

"""
    cal_dc_fll(U::F64, J::F64, N::F64)

Evaluate the double counting term by the fully localized limit scheme.

See also: [`cal_dc_amf`](@ref), [`cal_dc_exact`](@ref).
"""
function cal_dc_fll(U::F64, J::F64, N::F64)
    U * ( N - 0.5 ) - J / 2.0 * ( N - 1.0 )
end

"""
    cal_dc_amf(U::F64, J::F64, N::F64)

Evaluate the double counting term by the around mean-field scheme.

See also: [`cal_dc_fll`](@ref), [`cal_dc_exact`](@ref).
"""
function cal_dc_amf(U::F64, J::F64, N::F64)
    sorry()
end

"""
    cal_dc_exact(U::F64, J::F64, N::F64)

Evaluate the double counting term by the exact scheme.

See also: [`cal_dc_fll`](@ref), [`cal_dc_amf`](@ref).
"""
function cal_dc_exact(U::F64, J::F64, N::F64)
    sorry()
end
