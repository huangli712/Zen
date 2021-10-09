#
# Project : Pansy
# Source  : sigma.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/10/09
#

#=
### *Driver Functions*
=#

"""
    sigma_reset(ai::Array{Impurity,1}, with_init_dc::Bool = true)

Create initial self-energy functions and write them to `sigma.bare`. The
`sigma.bare` file is key input for the dynamical mean-field theory engine.
The word `bare` means that the double counting term has not been removed
from the self-energy functions. Now this function only supports Matsubara
self-energy functions Σ(𝑖ωₙ).

If `with_init_dc = true`, then the real parts of self-energy functions
are initialized by the double counting terms within the fully localized
limited scheme. If `with_init_dc = false`, then the self-energy functions
are set to be complex zero.

See also: [`sigma_dcount`](@ref).
"""
function sigma_reset(ai::Array{Impurity,1}, with_init_dc::Bool = true)
    # Print the header
    println("Sigma : Reset")
    println("Try to create bare self-energy functions")
    println("Current directory: ", pwd())

    # Extract some necessary parameters
    axis = get_m("axis")
    nmesh = get_m("nmesh")
    beta = get_m("beta")
    nsite = get_i("nsite")
    nspin = ( get_d("lspins") ? 2 : 1 )
    @assert nsite == length(ai)

    # Create frequency mesh
    println("Create frequency mesh")
    fmesh = zeros(F64, nmesh)
    if axis == 1 # Imaginary axis
        for i = 1:nmesh
            fmesh[i] = (2 * i - 1) * pi / beta
        end
        #
        println("  > Create Matsubara frequency mesh: $nmesh points")
    else # Real axis
        println("  > Create real frequency mesh: $nmesh points")
        sorry()
    end

    # Create default self-energy functions
    println("Create self-energy functions")
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

        # Setup initial values for the self-energy functions if needed
        if with_init_dc
            # Try to calculate the double counting terms within the fully
            # localized limited scheme at first.
            @assert ai[i].occup == get_i("occup")[i]
            sigdc, _ = cal_dc_fll(ai[i].upara, ai[i].jpara, ai[i].occup)

            # Go through spins, meshes, and bands, setup elements of S.
            # Only the diagonal elements in orbital space is filled.
            for s = 1:nspin
                for m = 1:nmesh
                    for b = 1:nband
                        S[b,b,m,s] = sigdc + 0.0im
                    end # END OF B LOOP
                end # END OF M LOOP
            end # END OF S LOOP
            #
            println("  > Initial value: [$i] -> ", sigdc)
        end

        # Push S into SA to save it
        push!(SA, S)
        #
        println("  > Shape of Array Σ: [$i] -> ", size(S))
    end # END OF I LOOP

    # Write self-energy functions and the corresponding frequency mesh
    println("Write self-energy functions")
    write_sigma(fmesh, SA, ai)
end

"""
    sigma_dcount(it::IterInfo,
                 ai::Array{Impurity,1},
                 reset_dc::Bool = false)

Calculate double counting terms for local self-energy functions and
write them to `sigma.dc`, which is an essential input for the dynamical
mean-field theory engine.

If `reset_dc = true`, it will reset the double counting terms to zero.
This is particularly useful for the first DFT + DMFT iteration. However,
if `reset_dc = false`, it will retain the double counting terms.

The field `it.dc` will be updated in this function as well.

See also: [`sigma_reset`](@ref).
"""
function sigma_dcount(it::IterInfo,
                      ai::Array{Impurity,1},
                      reset_dc::Bool = false)
    # Print the header
    println("Sigma : Dcount")
    println("Try to build double counting terms for self-energy functions")
    println("Current directory: ", pwd())

    # Extract some necessary parameters
    nsite = get_i("nsite")
    nspin = ( get_d("lspins") ? 2 : 1 )
    @assert nsite == length(ai)

    # Create double counting terms for self-energy functions
    println("Create double counting terms")
    #
    # Initialize an array for dc
    DCA = Array{F64,3}[]
    Edc = F64[]
    #
    # Go through the quantum impurity problems and calculate dc
    for i = 1:nsite
        # Get interaction parameters
        U = ai[i].upara
        J = ai[i].jpara

        # Determine impurity occupancy
        #
        # Get nominal occupation number
        N = get_i("occup")[i]
        #
        # Get realistic occupation number
        GetNimpx(ai[i])
        occup = ai[i].occup
        #
        # How about spin-resolved impurity occupancy
        nup = ai[i].nup
        ndown = ai[i].ndown

        # Get number of orbitals
        nband = ai[i].nband

        # Create a temporary array for the double counting terms
        DC = zeros(F64, nband, nband, nspin)

        # Choose suitable double counting scheme
        scheme = get_m("dcount")
        @cswitch scheme begin
            # Fully localized limit scheme (nominal occupation number)
            @case "fll1"
                sigup, sigdn, edc = cal_dc_fll(U, J, N / 2.0, N / 2.0)
                break

            # Fully localized limit scheme (dynamic occupation number)
            @case "fll2"
                sigup, sigdn, edc = cal_dc_fll(U, J, nup, ndown)
                break

            # Around mean-field scheme
            @case "amf"
                sigup, sigdn, edc = cal_dc_amf(U, J, nup, ndown, nband)
                break

            # K. Held scheme
            @case "held"
                sigup, edc = cal_dc_held(U, J, occup, nband)
                sigdn = sigup
                break

            # Exact double counting scheme
            @case "exact"
                sorry()
                break
        end
        #
        # Setup the DC arrays
        if nspin == 1
            fill!(DC, sigup)
        else
            fill!(view(DC,:,:,1), sigup)
            fill!(view(DC,:,:,2), sigdn)
        end
        #
        # Print some useful information
        println("  > Using the $scheme scheme: Σdc = $sigup (spin up)")
        println("  > Using the $scheme scheme: Σdc = $sigdn (spin down)")
        println("  > Using the $scheme scheme: Edc = $edc eV")
        println("  > Shape of Array Σdc: [$i] -> ", size(DC))

        # Special treatment for the first iteration
        if reset_dc && ( it.I₃ ≤ 1 && it.I₁ ≤ 1 )
            fill!(DC, 0.0)
            edc = 0.0
            println("  > Reset Σdc to: ", 0.0)
            println("  > Reset Edc to: ", 0.0)
        end

        # Use `sigup` to update the IterInfo struct
        it.dc[i] = sigup

        # Push DC into DCA to save it
        push!(DCA, DC)
        push!(Edc, edc)
    end # END OF I LOOP

    # Update the energy due to double counting
    it.et.dc = sum(Edc)

    # Write double counting terms
    println("Write double counting terms")
    write_sigdc(DCA, ai)
end

"""
    sigma_split(ai::Array{Impurity,1})

Split the hybridization functions (and local impurity levels) and then
distribute them into the `impurity.i` folder.

See also: [`sigma_gather`](@ref).
"""
function sigma_split(ai::Array{Impurity,1})
    # Print the header
    println("Sigma : Split")
    println("Try to split the hybridization functions and local impurity levels")
    println("Current directory: ", pwd())

    # Split the hybridization functions Δ
    println("Treat hybridization functions")
    #
    # Read it from dmft1/dmft.delta
    fmesh, Delta = read_delta(ai)
    #
    # Write it into impurity.i/dmft.delta
    write_delta(fmesh, Delta, ai)

    # Split the local impurity levels εᵢ
    println("Treat local impurity levels")
    #
    # Read it from dmft1/dmft.eimpx
    Eimpx = read_eimpx(ai)
    #
    # Write it into impurity.i/dmft.eimpx
    write_eimpx(Eimpx, ai)
end

"""
    sigma_gather(it::IterInfo, ai::Array{Impurity,1})

Gather the self-energy functions Σ (or similar local functions) from
the all the `impurity.i` folders and then combine them into a single
`sigma.bare` file.

See also: [`sigma_split`](@ref).
"""
function sigma_gather(it::IterInfo, ai::Array{Impurity,1})
    # Print the header
    println("Sigma : Gather")
    println("Try to combine self-energy functions from quantum impurities")
    println("Current directory: ", pwd())

    # Extract some necessary parameters
    nmesh = get_m("nmesh")
    nsite = get_i("nsite")
    @assert nsite == length(ai)

    # Create empty array for self-energy functions
    SA = Array{C64,4}[]

    # Declare frequency mesh
    fmesh = nothing

    # Go through each quantum impurity problems
    println("Gather self-energy functions")
    for t = 1:nsite
        # Extract the frequency mesh and self-energy function
        fmesh, sigma = GetSigma(ai[t])

        # Print some useful information
        println("  > Read self-energy functions for impurity: $t")
        println("  > Shape of Array ω: ", size(fmesh))
        println("  > Shape of Array Σ: ", size(sigma))

        # Extract and verify the dimensional parameters
        _, _b, _m, _ = size(sigma)
        @assert _b == ai[t].nband
        @assert _m == nmesh

        # Store sigma in SA
        push!(SA, sigma)
    end # END OF T LOOP

    # Now the self-energy functions for all quantum impurity problems and
    # the corresponding frequency mesh are ready. We are going to write
    # them into the `dmft1/sigma.bare` file.

    # Write self-energy functions to sigma.bare
    println("Write self-energy functions")
    write_sigma(fmesh, SA, ai)

    # Backup the self-energy functions. This operation is necessary
    # because we have to mix the self-energy functions for the two
    # successive runs.
    fsig = "dmft1/sigma.bare"
    cp(fsig, "$fsig.$(it.I₃).$(it.I₁)", force = true)
end

#=
### *Service Functions* : *For Double Counting Terms*
=#

#=
*Fully Localized Limit Scheme*:

For spin-unpolarized case,
```math
\begin{equation}
\Sigma^{\text{FLL}}_{\text{dc}}
    =
    U\left(N - \frac{1}{2}\right)
    -
    J \left(\frac{N}{2} - \frac{1}{2}\right).
\end{equation}
```

For spin-polarized case,
```math
\begin{equation}
\Sigma^{\text{FLL}}_{\text{dc},\uparrow}
    =
    U\left(N_{\uparrow} + N_{\downarrow} - \frac{1}{2}\right)
    -
    J \left(N_{\uparrow} - \frac{1}{2}\right),
\end{equation}
```

and

```math
\begin{equation}
\Sigma^{\text{FLL}}_{\text{dc},\downarrow}
    =
    U\left(N_{\uparrow} + N_{\downarrow} - \frac{1}{2}\right)
    -
    J \left(N_{\downarrow} - \frac{1}{2}\right).
\end{equation}
```

The energy due to double counting reads

```math
\begin{equation}
E_{\text{dc}} = \frac{U}{2}N(N-1) - \frac{J}{4}N(N-2).
\end{equation}
```
=#

"""
    cal_dc_fll(U::F64, J::F64, N::F64)

Evaluate the double counting term by the fully localized limit scheme.
This function is for the spin-unpolarized case.

See also: [`cal_dc_amf`](@ref), [`cal_dc_exact`](@ref).
"""
function cal_dc_fll(U::F64, J::F64, N::F64)
    Vdc = U * ( N - 0.5 ) - J / 2.0 * ( N - 1.0 )
    Edc = U / 2.0 * N * ( N - 1.0 ) - J / 4.0 * N * ( N - 2.0 )
    return Vdc, Edc
end

"""
    cal_dc_fll(U::F64, J::F64, Nup::F64, Ndn::F64)

Evaluate the double counting term by the fully localized limit scheme.
This function is for the spin-polarized case.

See also: [`cal_dc_amf`](@ref), [`cal_dc_exact`](@ref).
"""
function cal_dc_fll(U::F64, J::F64, Nup::F64, Ndn::F64)
    N = Nup + Ndn
    Vup = U * ( Nup + Ndn - 0.5 ) - J * ( Nup - 0.5 )
    Vdn = U * ( Nup + Ndn - 0.5 ) - J * ( Ndn - 0.5 )
    Edc = U / 2.0 * N * ( N - 1.0 ) - J / 4.0 * N * ( N - 2.0 )
    return Vup, Vdn, Edc
end

#=
*Around Mean-Field Scheme*:

The original formula are:

```math
\begin{equation}
\Sigma^{\text{AMF}}_{\text{dc},\uparrow}
    =
    UN_{\downarrow} + (U - J) N_{\uparrow} \frac{2l}{2l+1}.
\end{equation}
```

and

```math
\begin{equation}
\Sigma^{\text{AMF}}_{\text{dc},\downarrow}
    =
    UN_{\uparrow} + (U - J) N_{\downarrow} \frac{2l}{2l+1}.
\end{equation}
```

In the following codes, we implement a slightly different version:


```math
\begin{equation}
\Sigma^{\text{AMF}}_{\text{dc},\uparrow}
    =
    U\left(N_{\uparrow} + N_{\downarrow} - \frac{N_{\uparrow}}{M}\right)
    -
    J\left(N_{\uparrow} - \frac{N_{\uparrow}}{M}\right).
\end{equation}
```

and

```math
\begin{equation}
\Sigma^{\text{AMF}}_{\text{dc},\downarrow}
    =
    U\left(N_{\uparrow} + N_{\downarrow} - \frac{N_{\downarrow}}{M}\right)
    -
    J\left(N_{\downarrow} - \frac{N_{\downarrow}}{M}\right).
\end{equation}
```

In these equations, ``l`` means the quantum number of angular momentum,
while ``M`` is the number of correlated orbitals.

In addition, the energy due to double counting reads

```math
\begin{equation}
E_{\text{dc}}
    =
    \frac{UN^2}{2}\left(1 - \frac{1}{2M}\right)
    -
    \frac{JN^2}{2}\left(\frac{1}{2} - \frac{1}{2M}\right).
\end{equation}
```
=#

"""
    cal_dc_amf(U::F64, J::F64, N::F64, M::I64)

Evaluate the double counting term by the around mean-field scheme.
This function is for the spin-unpolarized case.

See also: [`cal_dc_fll`](@ref), [`cal_dc_exact`](@ref).
"""
function cal_dc_amf(U::F64, J::F64, N::F64, M::I64)
    Vdc = U * ( N - N / ( 2.0 * M ) ) - J * ( N / 2.0 - N / ( 2.0 * M ) )
    Edc = U * N^2 / 2.0 * ( 1.0 - 1.0 / ( 2.0 * M ) )
          - J * N^2 / 2.0 * ( 1.0 / 2.0 - 1.0 / ( 2.0 * M ) )
    return Vdc, Edc
end

"""
    cal_dc_amf(U::F64, J::F64, Nup::F64, Ndn::F64, M::I64)

Evaluate the double counting term by the around mean-field scheme.
This function is for the spin-polarized case.

See also: [`cal_dc_fll`](@ref), [`cal_dc_exact`](@ref).
"""
function cal_dc_amf(U::F64, J::F64, Nup::F64, Ndn::F64, M::I64)
    N = Nup + Ndn
    Vup = U * ( Nup + Ndn - Nup / M ) - J * ( Nup - Nup / M )
    Vdn = U * ( Nup + Ndn - Ndn / M ) - J * ( Ndn - Ndn / M )
    Edc = U * N^2 / 2.0 * ( 1.0 - 1.0 / ( 2.0 * M ) )
          - J * N^2 / 2.0 * ( 1.0 / 2.0 - 1.0 / ( 2.0 * M ) )
    return Vup, Vdn, Edc
end

#=
*K. Held's Double Counting Scheme*:

```math
\begin{equation}
\Sigma_{\text{dc}} = \bar{U} \left(N - \frac{1}{2}\right),
\end{equation}
```

where the averaged interaction ``\bar{U}`` reads

```math
\begin{equation}
\bar{U} = \frac{U + (M - 1)(2U - 5J)}{2M - 1}.
\end{equation}
```

Here ``M`` is the number of correlated orbitals.

The energy due to double counting reads

```math
\begin{equation}
E_{\text{dc}} = \frac{\bar{U}}{2} N (N - 1).
\end{equation}
```
=#

"""
    cal_dc_held(U::F64, J::F64, N::F64, M::I64)

Evaluate the double counting term by the K. Held scheme.

See also: [`cal_dc_fll`](@ref), [`cal_dc_amf`](@ref), [`cal_dc_exact`](@ref).
"""
function cal_dc_held(U::F64, J::F64, N::F64, M::I64)
    Uav = ( U + ( M - 1.0 ) * ( 2.0 * U - 5.0 * J ) ) / ( 2.0 * M - 1.0 )
    Vdc = Uav * ( N - 0.5 )
    Edc = Uav / 2.0 * N * ( N - 1.0 )
    return Vdc, Edc
end

"""
    cal_dc_exact(U::F64, J::F64, N::F64)

Evaluate the double counting term by the exact scheme.

See also: [`cal_dc_fll`](@ref), [`cal_dc_amf`](@ref).
"""
function cal_dc_exact(U::F64, J::F64, N::F64)
    sorry()
end

#=
### *Service Functions* : *Files I/O Operations*
=#

"""
    read_sigma(ai::Array{Impurity,1}, fsig::String = "dmft1/sigma.bare")

Read the self-energy functions from the `dmft1/sigma.bare` file. The
working directory of this function must be the root folder.

This function is usually called by `mixer_sigma()` function.

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

    # Print some useful information
    println("  > Read self-energy functions from: $fsig")
    println("  > Shape of Array ω: ", size(fmesh))
    for t in eachindex(SA)
        println("  > Shape of Array Σ: [$t] -> ", size(SA[t]))
    end

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

    # Print some useful information
    println("  > Read double counting terms from: $fsig")
    for t in eachindex(DCA)
        println("  > Shape of Array Σdc: [$t] -> ", size(DCA[t]))
    end

    # Return the desire array
    return DCA
end

#=
### *Service Functions* : *Files I/O Operations*
=#

"""
    write_sigma(fmesh::Array{F64,1},
                SA::Array{Array{C64,4},1},
                ai::Array{Impurity,1})

Write the self-energy functions and the corresponding frequency mesh into
the `dmft1/sigma.bare` file, which is key input for the dynamical mean-
field theory engine. The working directory of this function must be the
root folder.

See also: [`write_sigdc`](@ref).
"""
function write_sigma(fmesh::Array{F64,1},
                     SA::Array{Array{C64,4},1},
                     ai::Array{Impurity,1})
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

    # Print some useful information
    println("  > Write self-energy functions into: dmft1/sigma.bare")
    println("  > Shape of Array ω: ", size(fmesh))
    for t in eachindex(SA)
        println("  > Shape of Array Σ: [$t] -> ", size(SA[t]))
    end

    # Copy sigma.bare to the dmft2 directory
    cp("dmft1/sigma.bare", "dmft2/sigma.bare", force = true)

    # Print some useful information
    println("  > Write self-energy functions into: dmft2/sigma.bare")
    println("  > Shape of Array ω: ", size(fmesh))
    for t in eachindex(SA)
        println("  > Shape of Array Σ: [$t] -> ", size(SA[t]))
    end
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

    # Print some useful information
    println("  > Write double counting terms into: dmft1/sigma.dc")
    for t in eachindex(DCA)
        println("  > Shape of Array Σdc: [$t] -> ", size(DCA[t]))
    end

    # Copy sigma.dc to the dmft2 directory
    cp("dmft1/sigma.dc", "dmft2/sigma.dc", force = true)

    # Print some useful information
    println("  > Write double counting terms into: dmft2/sigma.dc")
    for t in eachindex(DCA)
        println("  > Shape of Array Σdc: [$t] -> ", size(DCA[t]))
    end
end
