#
# Project : Pansy
# Source  : mixer.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2025/03/26
#

#=
### *Multiple Dispatchers*
=#

"""
    mixer_call(::NULLMixer, it::IterInfo, ai::Array{Impurity,1})
    mixer_call(::ΣMixer, it::IterInfo, ai::Array{Impurity,1})
    mixer_call(::ΔMixer, it::IterInfo, ai::Array{Impurity,1})
    mixer_call(::EMixer, it::IterInfo, ai::Array{Impurity,1})
    mixer_call(::ΓMixer, it::IterInfo, ai::Array{Impurity,1})

Try to mix various functions during the iterations, including self-energy
functions Σ, hybridization functions Δ, effective impurity levels E, and
correction for density matrix Γ. It acts as a dispatcher.

See also: [`_mixer_`](@ref).
"""
function mixer_call(::NULLMixer, it::IterInfo, ai::Array{Impurity,1})
    sorry()
end
#
function mixer_call(::ΣMixer, it::IterInfo, ai::Array{Impurity,1})
    # Try to mix the self-energy functions Σ.
    mixer_sigma(it, ai)
end
#
function mixer_call(::ΔMixer, it::IterInfo, ai::Array{Impurity,1})
    # Try to mix the hybridization functions Δ.
    mixer_delta(it, ai)
end
#
function mixer_call(::EMixer, it::IterInfo, ai::Array{Impurity,1})
    # Try to mix the local impurity levels E.
    mixer_eimpx(it, ai)
end
#
function mixer_call(::ΓMixer, it::IterInfo, ai::Array{Impurity,1})
    # Try to mix the correction for density matrix Γ.
    mixer_gcorr(it)
end

#=
### *Driver Functions*
=#

"""
    mixer_sigma(it::IterInfo, ai::Array{Impurity,1})

Try to mix the self-energy functions Σ and then use the mixed values
to update the `dmft1/sigma.bare` file.

See also: [`mixer_core`](@ref), [`amix`](@ref), [`distance`](@ref).
"""
function mixer_sigma(it::IterInfo, ai::Array{Impurity,1})
    # Print the header
    println("Mixer : Sigma")
    println("Try to mix self-energy functions")
    println("Current directory: ", pwd())

    # Get current dmft loop
    cycle = it.I₃

    # Get current iteration
    curr = it.I₁

    # Get previous iteration
    if it.sc == 1
        _cycle, _prev = prev_it(it)
        @assert _cycle == cycle
        @assert _prev == curr - 1
        @assert _prev ≥ 1
    else
        _cycle, _prev = prev_it(it, 1)
        @assert cycle ≥ _cycle ≥ 1
        @assert _prev ≥ 1
    end
    println("Determine previous and current objects")
    println("  > Curr: (I₃, I₁) -> ($cycle, $curr)")
    println("  > Prev: (I₃, I₁) -> ($_cycle, $_prev)")

    # Determine filenames for self-energy functions
    fcurr = "dmft1/sigma.bare.$cycle.$curr"
    fprev = "dmft1/sigma.bare.$_cycle.$_prev"

    # Check whether these files are available
    @assert isfile(fcurr) && isfile(fprev)

    # Read in the self-energy functions (previous and current)
    println("Read self-energy functions")
    Mcurr, Scurr = read_sigma(ai, fcurr)
    Mprev, Sprev = read_sigma(ai, fprev)
    @assert size(Scurr) == size(Sprev) && size(Mcurr) == size(Mprev)

    # Mix the self-energy functions using linear mixing algorithm
    println("Mix self-energy functions from two successive iterations")
    α = amix(it)
    Snew = Scurr * α + Sprev * (1.0 - α)
    println("  > Mixing parameter α = $α")

    # Write the new self-energy functions into `dmft1/sigma.bare`
    println("Write self-energy functions")
    write_sigma(Mcurr, Snew, ai)

    # Check the convergence condition
    println("Evaluate the convergence condition for self-energy functions")
    dist = distance(Scurr, Sprev)
    it.cs = ( dist < get_m("sc") )
    println("  > Averaged ΔΣ = $dist ( convergence is $(it.cs) )")
end

"""
    mixer_delta(it::IterInfo, ai::Array{Impurity,1})

Try to mix the hybridization functions Δ and then use the mixed values
to update the `dmft1/dmft.delta` file.

See also: [`mixer_core`](@ref), [`amix`](@ref).
"""
function mixer_delta(it::IterInfo, ai::Array{Impurity,1})
    # Print the header
    println("Mixer : Delta")
    println("Try to mix hybridization functions")
    println("Current directory: ", pwd())

    # Get current dmft loop
    cycle = it.I₃

    # Get current iteration
    curr = it.I₁

    # Get previous iteration
    if it.sc == 1
        _cycle, _prev = prev_it(it)
        @assert _cycle == cycle
        @assert _prev == curr - 1
        @assert _prev ≥ 1
    else
        _cycle, _prev = prev_it(it, 1)
        @assert cycle ≥ _cycle ≥ 1
        @assert _prev ≥ 1
    end
    println("Determine previous and current objects")
    println("  > Curr: (I₃, I₁) -> ($cycle, $curr)")
    println("  > Prev: (I₃, I₁) -> ($_cycle, $_prev)")

    # Determine filenames for hybridization functions
    fcurr = "dmft1/dmft.delta.$cycle.$curr"
    fprev = "dmft1/dmft.delta.$_cycle.$_prev"

    # Check whether these files are available
    @assert isfile(fcurr) && isfile(fprev)

    # Read in the hybridization functions (previous and current)
    println("Read hybridization functions")
    Mcurr, Dcurr = read_delta(ai, fcurr)
    Mprev, Dprev = read_delta(ai, fprev)
    @assert size(Dcurr) == size(Dprev) && size(Mcurr) == size(Mprev)

    # Mix the hybridization functions using linear mixing algorithm
    println("Mix hybridization functions from two successive iterations")
    α = amix(it)
    Dnew = Dcurr * α + Dprev * (1.0 - α)
    println("  > Mixing parameter α = $α")

    # Write the new hybridization functions into `dmft1/dmft.delta`
    println("Write hybridization functions")
    write_delta(Mcurr, Dnew, ai, "dmft1/dmft.delta")
end

"""
    mixer_eimpx(it::IterInfo, ai::Array{Impurity,1})

Try to mix the local impurity levels εᵢ and then use the mixed value
to update the `dmft1/dmft.eimpx` file.

See also: [`mixer_core`](@ref), [`amix`](@ref).
"""
function mixer_eimpx(it::IterInfo, ai::Array{Impurity,1})
    # Print the header
    println("Mixer : Eimpx")
    println("Try to mix local impurity levels")
    println("Current directory: ", pwd())

    # Get current dmft loop
    cycle = it.I₃

    # Get current iteration
    curr = it.I₁

    # Get previous iteration
    if it.sc == 1
        _cycle, _prev = prev_it(it)
        @assert _cycle == cycle
        @assert _prev == curr - 1
        @assert _prev ≥ 1
    else
        _cycle, _prev = prev_it(it, 1)
        @assert cycle ≥ _cycle ≥ 1
        @assert _prev ≥ 1
    end
    println("Determine previous and current objects")
    println("  > Curr: (I₃, I₁) -> ($cycle, $curr)")
    println("  > Prev: (I₃, I₁) -> ($_cycle, $_prev)")

    # Determine filenames for local impurity levels
    fcurr = "dmft1/dmft.eimpx.$cycle.$curr"
    fprev = "dmft1/dmft.eimpx.$_cycle.$_prev"

    # Check whether these files are available
    @assert isfile(fcurr) && isfile(fprev)

    # Read in the local impurity levels (previous and current)
    println("Read local impurity levels")
    Ecurr = read_eimpx(ai, fcurr)
    Eprev = read_eimpx(ai, fprev)
    @assert size(Ecurr) == size(Eprev)

    # Mix the local impurity levels using linear mixing algorithm
    println("Mix local impurity levels from two successive iterations")
    α = amix(it)
    Enew = Ecurr * α + Eprev * (1.0 - α)
    println("  > Mixing parameter α = $α")

    # Write the new local impurity levels into `dmft1/dmft.eimpx`
    println("Write local impurity levels")
    write_eimpx(Enew, ai, "dmft1/dmft.eimpx")
end

#=
*Kerker mixing algorithm* :

```math
\begin{equation}
\Gamma_{\text{mix}}(\mathbf{G})
    =
    \Gamma_{\text{in}}(\mathbf{G})
    +
    \alpha \frac{\mathbf{G}^2}{\mathbf{G}^2 + \gamma^2}
    [\Gamma_{\text{out}}(\mathbf{G}) - \Gamma_{\text{in}}(\mathbf{G})],
\end{equation}
```

where ``\mathbf{G}`` is the ``k`` vector, and ``\alpha`` and ``\gamma``
are two predefined parameters. In the present implementation, we set
``\alpha = 0.1`` and ``\gamma = 1.0``. Usually they work quite well.
=#

"""
    mixer_gcorr(it::IterInfo)

Try to mix the correction for density matrix Γ and then use the mixed value
to update the `dmft2/dmft.gcorr` file. Here we use the Kerker algorithm,
instead of the linear mixing algorithm.

See also: [`mixer_core`](@ref), [`distance`](@ref).
"""
function mixer_gcorr(it::IterInfo)
    # Print the header
    println("Mixer : Gamma")
    println("Try to mix correction for density matrix")
    println("Current directory: ", pwd())

    # Get current dmft loop
    cycle = it.I₃

    # Get current iteration
    curr = it.I₂

    # Get previous iteration
    @assert it.sc == 2
    _cycle, _prev = prev_it(it, 2)
    @assert cycle ≥ _cycle ≥ 1
    @assert _prev ≥ 1
    #
    println("Determine previous and current objects")
    println("  > Curr: (I₃, I₂) -> ($cycle, $curr)")
    println("  > Prev: (I₃, I₂) -> ($_cycle, $_prev)")

    # Determine filenames for correction for density matrix
    fcurr = "dmft2/dmft.gcorr.$cycle.$curr"
    fprev = "dmft2/dmft.gcorr.$_cycle.$_prev"

    # Check whether these files are available
    @assert isfile(fcurr) && isfile(fprev)

    # Read in the correction for density matrix (previous and current)
    println("Read correction for density matrix")
    kmesh_curr, kwin_curr, gcorr_curr = read_gcorr(fcurr)
    kmesh_prev, kwin_prev, gcorr_prev = read_gcorr(fprev)
    @assert size(kmesh_curr) == size(kmesh_prev)
    @assert size(kwin_curr) == size(kwin_prev)
    if size(gcorr_curr) != size(gcorr_prev)
        print("  > Size of density matrix does not match each other")
        println(red(" (Very dangerous)"))
        return
    end

    # Mix the correction for density matrix using Kerker algorithm
    println("Mix correction for density matrix from two successive iterations")
    #
    # Extract and setup key parameters
    _, _, nkpt, nspin = size(gcorr_curr)
    α = 0.1 # Mixing parameters
    γ = 1.0
    #
    # Apply the Kerker algorithm
    for s = 1:nspin
        for k = 1:nkpt
            # Evaluate the mixing factor
            G₂ = sum(kmesh_curr[k,:] .^ 2)
            amix = α * G₂ / (G₂ + γ^2)
            @printf("  > Mixing parameter α = %10.7f (for 𝑘 %4i and σ %1i)\n", amix, k, s)
            #
            # Create a view for the diagonal elements only
            ind = diagind(gcorr_curr[:,:,k,s])
            Γcurr = view(view(gcorr_curr,:,:,k,s), ind)
            Γprev = view(view(gcorr_prev,:,:,k,s), ind)
            #
            # Mix the diagonal elements only
            @. Γcurr = amix * Γcurr + (1.0 - amix) * Γprev
        end # END OF K LOOP
    end # END OF S LOOP

    # Write the new correction for density matrix into `dmft2/dmft.gcorr`
    println("Write correction for density matrix")
    write_gcorr(kmesh_curr, kwin_curr, gcorr_curr, "dmft2/dmft.gcorr")

    # Check the convergence condition
    println("Evaluate the convergence condition for density matrix")
    dist = distance(gcorr_curr, gcorr_prev)
    it.cc = ( dist < get_m("cc") )
    println("  > Averaged ΔΓ = $dist ( convergence is $(it.cc) )")
end

#=
### *Service Functions*
=#

"""
    amix(it::IterInfo)

Return the mixing factor for mixer component. It should depend on the
current iteration number, instead of being a constant.

See also: [`IterInfo`](@ref).
"""
function amix(it::IterInfo)
    factor = 1.0
    #
    if it.sc == 1
        factor = exp(-(it.I₁ - 1) * get_m("mixer"))
    else
        factor = exp(-(it.I₃ - 1) * get_m("mixer"))
    end
    #
    return factor
end

"""
    distance(SA::Vector{Array{C64,4}}, SB::Vector{Array{C64,4}})

Calculate the difference between two multi-dimensional arrays. Usually
We apply this function to calculate the difference between two
self-energy functions.

See also: [`mixer_sigma`](@ref).
"""
function distance(SA::Vector{Array{C64,4}}, SB::Vector{Array{C64,4}})
    # Check the dimensional parameters to make sure SA is similar to SB
    @assert length(SA) == length(SB)
    foreach((A, B) -> ( @assert size(A) == size(B) ), SA, SB)

    # Evaluate the difference
    SC = SA - SB
    diff = zero(C64)
    for i in eachindex(SC)
        # Actually, the non-diagonal elements are zero!
        num_zeros_elements = count(x -> x == zero(C64), SC[i])
        diff = diff + sum(SC[i]) / (length(SC[i]) - num_zeros_elements)
    end
    diff = diff / length(SC)

    # Return the desired value
    return abs(diff)
end

"""
    distance(GA::Array{C64,4}, GB::Array{C64,4})

Calculate the difference between two multi-dimensional arrays. Usually
We apply this function to calculate the difference between two
corrections for density matrix.

See also: [`mixer_gcorr`](@ref).
"""
function distance(GA::Array{C64,4}, GB::Array{C64,4})
    # Check the dimensional parameters to make sure GA is similar to GB
    @assert size(GA) == size(GB)

    # Evaluate the difference
    GC = GA - GB
    diff = sum(GC) / length(GC)

    # Return the desired value
    return abs(diff)
end
