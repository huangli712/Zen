#
# Project : Pansy
# Source  : sigma.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/06/08
#

#
# Driver Functions
#

"""
    sigma_reset(ai::Array{Impurity,1})

Create initial self-energy functions and write them to `sigma.bare`. The
`sigma.bare` file is key input for the dynamical mean-field theory engine.
Now this function only supports Matsubara self-energy functions.

See also: [`sigma_dcount`](@ref).
"""
function sigma_reset(ai::Array{Impurity,1})
    # Print the log
    println("Sigma : Reset")

    # Extract some necessary parameters
    axis = get_m("axis")
    nmesh = get_m("nmesh")
    beta = get_m("beta")
    nsite = get_i("nsite")
    nspin = ( get_d("lspins") ? 2 : 1 )
    @assert nsite == length(ai)

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

    # Create default self-energy functions
    #
    # Initialize an array for self-energy functions
    SA = Array{C64,4}[]
    #
    # Go through the quantum impurity problems
    for i = 1:nsite
        # Get the dimension of impurity problem
        nband = ai[i].nband

        # Create a temporary array for self-energy function
        S = zeros(C64, nband, nband, nmesh, nspin)

        # Push S into SA to save it
        push!(SA, S)
    end # END OF I LOOP
    println("  Create local self-energy functions")

    # Write self-energy functions and the corresponding frequency mesh
    write_sigma(fmesh, SA, ai)

    # Print blank line for better visualization
    println()
end

"""
    sigma_dcount(it::IterInfo, ai::Array{Impurity,1})

Calculate double counting terms for local self-energy functions and
write them to `sigma.dc`, which is key input for the dynamical mean-
field theory engine.

The field `it.dc` will be updated in this function as well.

See also: [`sigma_reset`](@ref).
"""
function sigma_dcount(it::IterInfo, ai::Array{Impurity,1})
    # Print the log
    println("Sigma : Dcount")

    # Extract some necessary parameters
    nsite = get_i("nsite")
    nspin = ( get_d("lspins") ? 2 : 1 )
    @assert nsite == length(ai)

    # Create double counting terms for self-energy functions
    #
    # Initialize an array for dc
    DCA = Array{F64,3}[]
    #
    # Go through the quantum impurity problems and calculate dc
    for i = 1:nsite
        # Get interaction parameters
        U = ai[i].upara
        J = ai[i].jpara

        # Get nominal occupation number
        N = get_i("occup")[i]

        # Get realistic occupation number
        GetNimpx(ai[i])
        occup = ai[i].occup

        # Get number of orbitals
        nband = ai[i].nband

        # Create a temporary array for the double counting terms
        DC = zeros(F64, nband, nband, nspin)

        # Choose suitable double counting scheme
        @cswitch get_m("dcount") begin
            # Fully localized limit scheme with fixed occupation number
            @case "fll1"
                sigdc = cal_dc_fll(U, J, N)
                fill!(DC, sigdc)
                break

            # Fully localized limit scheme with dynamic occupation number
            @case "fll2"
                sigdc = cal_dc_fll(U, J, occup)
                fill!(DC, sigdc)
                break

            # Around mean-field scheme
            @case "amf"
                sorry()
                break

            # K. Held scheme
            @case "held"
                sorry()
                break
            
            # Exact double counting scheme
            @case "exact"
                sorry()
                break
        end

        # Special treatment for the first iteration
        if it.I₃ <= 1 && it.I₁ <= 1
            sigdc = 0.0
            fill!(DC, sigdc)
        end

        # Use `sigdc` to update the IterInfo struct
        it.dc[i] = DC[1,1,1]

        # Push DC into DCA to save it
        push!(DCA, DC)
    end # END OF I LOOP
    println("  Create double counting terms: $(get_m("dcount"))")

    # Write double counting terms
    write_sigdc(DCA, ai)

    # Print blank line for better visualization
    println()
end

"""
    sigma_split(ai::Array{Impurity,1})

Split the hybridization functions (and local impurity levels) and then
distribute them into the `impurity.i` folder.

See also: [`sigma_gather`](@ref).
"""
function sigma_split(ai::Array{Impurity,1})
    # Print the log
    println("Sigma : Split")

    # Split the hybridization functions Δ
    fmesh, Delta = read_delta(ai)
    write_delta(fmesh, Delta, ai)

    # Split the local impurity levels εᵢ
    Eimpx = read_eimpx(ai)
    write_eimpx(Eimpx, ai)

    # Print blank line for better visualization
    println()
end

"""
    sigma_gather(it::IterInfo, ai::Array{Impurity,1})

Gather the self-energy functions Σ (or similar local functions) from
the all the `impurity.i` folders and then combine them into a single
`sigma.bare` file.

See also: [`sigma_split`](@ref).
"""
function sigma_gather(it::IterInfo, ai::Array{Impurity,1})
    # Print the log
    println("Sigma : Gather")

    # Extract some necessary parameters
    nmesh = get_m("nmesh")
    nsite = get_i("nsite")
    @assert nsite == length(ai)

    # Create empty array for self-energy functions
    SA = Array{C64,4}[]

    # Declare frequency mesh
    fmesh = nothing

    # Go through each quantum impurity problems
    for t = 1:nsite
        # Extract the frequency mesh and self-energy function
        fmesh, sig_l = GetSigma(ai[t])
        println("  Read self-energy functions for impurity: $t")

        # Extract and verify the dimensional parameters
        _, _b, _m, _ = size(sig_l)
        @assert _b == ai[t].nband
        @assert _m == nmesh

        # Store sig_l in SA
        push!(SA, sig_l)
    end # END OF T LOOP

    # Now the self-energy functions for all quantum impurity problems and
    # the corresponding frequency mesh are ready. We are going to write
    # them into the `dmft1/sigma.bare` file.

    # Write self-energy functions to sigma.bare
    write_sigma(fmesh, SA, ai)

    # Backup the self-energy functions
    fsig = "dmft1/sigma.bare"
    cp(fsig, "$fsig.$(it.I₃).$(it.I₁)", force = true)

    # Print blank line for better visualization
    println()
end

#
# Service Functions: For Double Counting Terms
#

#=
*Theory*:

```math
\Sigma^{\text{FLL}}_{\text{dc}}
    =
    U\left(N - \frac{1}{2}\right) - \frac{J}{2} (N - 1).
```
=#

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
    cal_dc_held()

Evaluate the double counting term by the K. Held scheme.

See also: [`cal_dc_fll`](@ref), [`cal_dc_exact`](@ref).
"""
function cal_dc_held()
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

#
# Service Functions: For I/O Operations
#

"""
    read_sigma(ai::Array{Impurity,1}, fsig::String = "dmft1/sigma.bare")

Read the self-energy functions from the `dmft1/sigma.bare` file. The
working directory of this function must be the root folder.

See also: [`read_sigdc`](@ref).
"""
function read_sigma(ai::Array{Impurity,1}, fsig::String = "dmft1/sigma.bare")
    # Declare the frequency mesh and self-energy functions
    fmesh = nothing
    SA = nothing

    # Make sure the existence of self-energy functions
    @assert isfile(fsig)

    # Parse `fsig`, extract the self-energy functions
    open(fsig, "r") do fin
        # Skip some lines
        readline(fin)
        readline(fin)
        readline(fin)

        # Get the dimensional parameters
        #
        # Skip two lines
        readline(fin)
        readline(fin)
        #
        nsite = parse(I64, line_to_array(fin)[3])
        nmesh = parse(I64, line_to_array(fin)[3])
        nspin = parse(I64, line_to_array(fin)[3])
        #
        for t = 1:nsite
            qdim  = parse(I64, line_to_array(fin)[3])
            @assert qdim == ai[t].nband
        end
        #
        readline(fin)

        # Create an array for frequency mesh
        fmesh = zeros(F64, nmesh)

        # Create an array for self-energy functions
        SA = Array{C64,4}[]

        # Read the body
        # Go through each quantum impurity problem
        for t = 1:nsite
            # Create an array for the site-dependent self-energy functions
            Sigma = zeros(C64, ai[t].nband, ai[t].nband, nmesh, nspin)

            # Go through each spin
            for s = 1:nspin
                # Parse indices and dimensional parameter
                strs = readline(fin)
                _t = parse(I64, line_to_array(strs)[3])
                _s = parse(I64, line_to_array(strs)[5])
                @assert _t == t
                @assert _s == s

                # Go through each frequency point
                for m = 1:nmesh
                    # Parse frequency mesh
                    fmesh[m] = parse(F64, line_to_array(fin)[3])
                    # Parse self-energy functions
                    for q = 1:ai[t].nband
                        for p = 1:ai[t].nband
                            _re, _im = parse.(F64, line_to_array(fin)[1:2])
                            Sigma[p,q,m,s] = _re + _im * im
                        end
                    end
                end # END OF M LOOP

                # Skip separator
                readline(fin)
            end # END OF S LOOP

            # Store Sigma in SA
            push!(SA, Sigma)
        end # END OF T LOOP
    end # END OF IOSTREAM
    println("  Read self-energy functions from: $fsig")

    # Return the desired arrays
    return fmesh, SA
end

"""
    read_sigdc(ai::Array{Impurity,1}, fsig::String = "dmft1/sigma.dc")

Read the double counting terms from the `dmft1/sigma.dc` file. The
working directory of this function must be the root folder.

See also: [`read_sigma`](@ref).
"""
function read_sigdc(ai::Array{Impurity,1}, fsig::String = "dmft1/sigma.dc")
    # Declare the double counting terms
    DCA = nothing

    # Make sure the existence of double counting terms
    @assert isfile(fsig)

    # Parse `fsig`, extract the double counting terms
    open(fsig, "r") do fin
        # Get the dimensional parameters
        #
        # Skip three lines
        readline(fin)
        readline(fin)
        readline(fin)
        #
        nsite = parse(I64, line_to_array(fin)[3])
        nspin = parse(I64, line_to_array(fin)[3])
        #
        for t = 1:nsite
            qdim  = parse(I64, line_to_array(fin)[3])
            @assert qdim == ai[t].nband
        end
        #
        readline(fin)

        # Create an array for double counting terms
        DCA = Array{F64,3}[]

        # Read the body
        # Go through each quantum impurity problem
        for t = 1:nsite
            # Create an array for the site-dependent double counting terms
            DC = zeros(F64, ai[t].nband, ai[t].nband, nspin)

            # Go through each spin
            for s = 1:nspin
                # Parse indices and dimensional parameter
                strs = readline(fin)
                _t = parse(I64, line_to_array(strs)[3])
                _s = parse(I64, line_to_array(strs)[5])
                @assert _t == t
                @assert _s == s

                # Parse self-energy functions
                for q = 1:ai[t].nband
                    for p = 1:ai[t].nband
                        DC[p,q,s] = parse.(F64, line_to_array(fin)[1])
                    end
                end

                # Skip separator
                readline(fin)
            end # END OF S LOOP

            # Store DC in DCA
            push!(DCA, DC)
        end # END OF T LOOP
    end # END OF IOSTREAM
    println("  Read double counting terms from: $fsig")

    # Return the desire array
    return DCA
end

#
# Service Functions: For I/O Operations
#

"""
    write_sigma(fmesh::Array{F64,1}, SA::Array{Array{C64,4},1}, ai::Array{Impurity,1})

Write the self-energy functions and the corresponding frequency mesh into
the `dmft1/sigma.bare` file, which is key input for the dynamical mean-
field theory engine. The working directory of this function must be the
root folder.

See also: [`write_sigdc`](@ref).
"""
function write_sigma(fmesh::Array{F64,1}, SA::Array{Array{C64,4},1}, ai::Array{Impurity,1})
    # Extract some necessary parameters
    axis = get_m("axis")
    nmesh = get_m("nmesh")
    beta = get_m("beta")
    nsite = get_i("nsite")
    nspin = ( get_d("lspins") ? 2 : 1 )

    # Write self-energy functions to dmft1/sigma.bare
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
            println(fout, "ndim$i -> $(ai[i].nband)")
        end
        println(fout)

        # Write the body
        # Go through each quantum impurity problem and spin
        for i = 1:nsite
            for s = 1:nspin
                println(fout, "# site: $i spin: $s")
                for m = 1:nmesh
                    # Write frequency grid
                    @printf(fout, "%4s %6i %16.12f\n", "w:", m, fmesh[m])
                    # Write self-energy function data
                    # There are 2 columns and nband * nband rows
                    for p = 1:ai[i].nband
                        for q = 1:ai[i].nband
                            x = SA[i][q, p, m, s]
                            @printf(fout, "%16.12f %16.12f\n", real(x), imag(x))
                        end
                    end
                end # END OF M LOOP
                # Write separator for each frequency point
                println(fout)
            end # END OF S LOOP
        end # END OF I LOOP
    end # END OF IOSTREAM

    # Print message to the screen
    println("  Write self-energy functions into: dmft1/sigma.bare")
end

"""
    write_sigdc(DCA::Array{Array{F64,3},1}, ai::Array{Impurity,1})

Write the double counting terms into the `dmft1/sigma.dc` file, which is
the key input for the dynamical mean-field theory engine. The working
directory of this function must be the root folder.

See also: [`write_sigma`](@ref).
"""
function write_sigdc(DCA::Array{Array{F64,3},1}, ai::Array{Impurity,1})
    # Extract some necessary parameters
    nsite = get_i("nsite")
    nspin = ( get_d("lspins") ? 2 : 1 )
    @assert nsite == length(ai)

   # Write double counting terms to sigma.dc
    open("dmft1/sigma.dc", "w") do fout
        # Write the header
        println(fout, "# File: sigma.dc")
        println(fout, "# Data: double counting terms")
        println(fout)
        println(fout, "nsite -> $nsite")
        println(fout, "nspin -> $nspin")
        for i = 1:nsite
            println(fout, "ndim$i -> $(ai[i].nband)")
        end
        println(fout)

        # Write the body
        # Go through each impurity problem
        for i = 1:nsite
            for s = 1:nspin
                println(fout, "# site: $i spin: $s")
                # There are 2 columns and nband * nband rows
                # The double counting terms are assumed to be complex
                # numbers but with zero imaginary parts.
                for p = 1:ai[i].nband
                    for q = 1:ai[i].nband
                        # Only the diagonal elements are useful
                        if p == q
                            @printf(fout, "%16.12f %16.12f\n", DCA[i][q, p, s], 0.0)
                        else
                            @printf(fout, "%16.12f %16.12f\n", 0.0, 0.0)
                        end
                    end # END OF Q LOOP
                end # END OF P LOOP
                # Write separator for each site and spin block
                println(fout)
            end # END OF S LOOP
        end # END OF I LOOP
    end # END OF IOSTREAM

    # Print message to the screen
    println("  Write double counting terms into: dmft1/sigma.dc")
end
