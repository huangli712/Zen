#
# Project : Pansy
# Source  : sigma.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/05/25
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
    nspin = ( get_d("lspins") ? 2 : 1 )

    # Create frequency mesh
    fmesh = zeros(F64, nmesh)
    if axis === 1 # Imaginary axis
        for i = 1:nmesh
            fmesh[i] = (2 * i - 1) * pi / beta
        end
        println("  Create Matsubara frequency mesh")
    else # Real axis
        sorry()
    end

    # Create self-energy functions
    #
    # Initialize an array for self-energy functions
    SA = Array{C64,4}[]
    #
    # `D` is used to record the dimensional parameter of the self-energy
    D = I64[]
    #
    # Go through the impurity problems
    for i = 1:nsite
        # Retrieve specification for impurity problem
        str = get_i("shell")[i]

        # Get the dimension of impurity problem
        ndim = get(sdim, str, 1)
        #
        # Save `ndim` in `D`
        push!(D, ndim)

        # Create a temporary array for self-energy function
        S = zeros(C64, ndim, ndim, nmesh, nspin)

        # Push S into SA to save it
        push!(SA, S)
    end
    println("  Create local self-energy functions")

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
                    @printf(fout, "%4s %6i %16.12f\n", "w:", m, fmesh[m])
                    # There are 2 columns and ndim * ndim rows
                    for a = 1:D[i]
                        for b = 1:D[i]
                            x = SA[i][b, a, m, s]
                            @printf(fout, "%16.12f %16.12f\n", real(x), imag(x))
                        end
                    end
                end
                println(fout)
            end
        end
    end

    # Print blank line for better visualization
    println("The local self-energy functions are written to sigma.bare")
    println()
end

"""
    sigma_dcount(it::IterInfo)

Calculate double counting terms for local self-energy functions and
write them to `sigma.dc`.

See also: [`sigma_reset`](@ref).
"""
function sigma_dcount(it::IterInfo)
    # Print the log
    println("Sigma : Dcount")

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
    nspin = ( get_d("lspins") ? 2 : 1 )

    # Create double counting terms for self-energy functions
    #
    # Initialize an array for dc
    DCA = Array{F64,3}[]
    #
    # `D` is used to record the dimensional parameter of the self-energy
    D = I64[]
    #
    # Go through the impurity problems and calculate dc
    for i = 1:nsite
        # Get interaction parameters
        U = get_i("upara")[i]
        J = get_i("jpara")[i]

        # Get occupation numbers
        N = get_i("occup")[i]

        # Get number of orbitals: `ndim`
        #
        # Retrieve specification for impurity problem
        str = get_i("shell")[i]
        #
        # Get the dimension of impurity problem
        ndim = get(sdim, str, 1)
        #
        # Save `ndim` in `D`
        push!(D, ndim)

        # Create a temporary array for the double counting terms
        DC = zeros(F64, ndim, ndim, nspin)

        # Choose suitable double counting scheme
        @cswitch get_m("dcount") begin
            # Fully localized limit scheme with fixed occupation number
            @case "fll1"
                sigdc = cal_dc_fll(U, J, N)
                fill!(DC, sigdc)
                break

            # Fully localized limit scheme with dynamic occupation number
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

        # Special treatment for the first iteration
        @show it
        if it.dmft_cycle <= 1 && it.dmft1_iter <= 1
            fill!(DC, 0.0)
        end

        # Push DC into DCA to save it
        push!(DCA, DC)
    end
    println("  Create double counting terms: $(get_m("dcount"))")

    # Write double counting terms to sigma.dc
    open("dmft1/sigma.dc", "w") do fout
        # Write the header
        println(fout, "# File: sigma.dc")
        println(fout, "# Data: double counting terms")
        println(fout)
        println(fout, "nsite -> $nsite")
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
                # There are 2 columns and ndim * ndim rows
                # The double counting terms are assumed to be complex
                # numbers with zero imaginary parts.
                for a = 1:D[i]
                    for b = 1:D[i]
                        if a == b
                            @printf(fout, "%16.12f %16.12f\n", DCA[i][b, a, s], 0.0)
                        else
                            @printf(fout, "%16.12f %16.12f\n", 0.0, 0.0)
                        end
                    end
                end
                println(fout)
            end
        end
    end

    # Print blank line for better visualization
    println("The double counting terms are written to sigma.dc")
    println()
end

"""
    sigma_split()

Split the hybridization functions (or similar local functions) and then
distribute them into the `impurity.i` folder.

See also: [`sigma_gather`](@ref).
"""
function sigma_split()
    # Print the log
    println("Sigma : Split")

    # Declare the frequency mesh and hybridization function
    fmesh = []
    Delta = []
    ndim  = []

    # Filename for hybridization functions
    fhyb = "dmft1/dmft.hyb_l"

    # Make sure the existence of hybridization functions
    @assert isfile(fhyb)

    # Parse `fhyb`, extract the hybridization functions 
    open(fhyb, "r") do fin

        # Get the dimensional parameters
        nsite = parse(I64, line_to_array(fin)[3])
        nspin = parse(I64, line_to_array(fin)[3])
        nmesh = parse(I64, line_to_array(fin)[3])
        qdim = parse(I64, line_to_array(fin)[4])

        # Skip two lines
        readline(fin)
        readline(fin)

        # Create an array for frequency mesh
        fmesh = zeros(F64, nmesh)

        # Create an array for hybridization functions
        Delta = zeros(C64, qdim, qdim, nmesh, nspin, nsite)
        ndim = zeros(I64, nsite)

        # Read the data
        for t = 1:nsite
            for s = 1:nspin
                # Parse indices and dimensional parameter
                strs = readline(fin)
                _t = parse(I64, line_to_array(strs)[3])
                _s = parse(I64, line_to_array(strs)[5])
                cdim = parse(I64, line_to_array(strs)[7])
                ndim[t] = cdim
                @assert _t == t && _s == s
                for m = 1:nmesh
                    # Parse frequency mesh
                    fmesh[m] = parse(F64, line_to_array(fin)[3])
                    # Parse hybridization functions
                    for q = 1:cdim
                        for p = 1:cdim
                            _re, _im = parse.(F64, line_to_array(fin)[3:4])
                            Delta[p,q,m,s,t] = _re + _im * im
                        end
                    end
                end
                # Skip two lines
                readline(fin)
                readline(fin)
            end
        end
    end

    # Next, we are going to split the hybridization functions according
    # to the quantum impurity problems

    # Extract the dimensional parameters
    _, qdim, nmesh, nspin, nsite = size(Delta)

    # Go through each quantum impurity problems
    for t = 1:nsite

        fhyb = "impurity.$t/dmft.hyb_l"

        open(fhyb, "w") do fout
            @printf(fout, "%9s %4i\n", "# nsite: ", nsite)
            @printf(fout, "%9s %4i\n", "# nspin: ", nspin)
            @printf(fout, "%9s %4i\n", "# nmesh: ", nmesh)
            @printf(fout, "%9s %4i\n", "# qdim : ", qdim)

            println(fout)
            println(fout)

            for s = 1:nspin
                @printf(fout, "# site:%4i  spins:%4i  dims:%4i\n", t, s, ndim[t])
                for m = 1:nmesh
                    @printf(fout, "w:%6i %16.8f\n", m, fmesh[m])
                    for q = 1:ndim[t]
                        for p = 1:ndim[t]
                            z = Delta[p,q,m,s,t]
                            @printf(fout, "%4i %4i %16.8f %16.8f\n", p, q, real(z), imag(z))
                        end
                    end
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
    sigma_gather()

Gather the self-energy functions (or similar local functions) from the
`impurity.i` folder and then combine them into a single `sigma.bare` file.

See also: [`sigma_split`](@ref).
"""
function sigma_gather()
    # Print the log
    println("Sigma : Gather")

    # Print blank line for better visualization
    println()
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
