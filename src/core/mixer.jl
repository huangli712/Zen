#
# Project : Pansy
# Source  : mixer.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/06/25
#

"""
    mixer_sigma(it::IterInfo, ai::Array{Impurity,1})

Try to mix the self-energy functions Σ and then use the mixed values
to update the `dmft1/sigma.bare` file.

See also: [`mixer_core`](@ref).
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
    @printf("  > Curr: (I₃, I₁) -> (%4i,%4i)\n", cycle, curr)
    @printf("  > Prev: (I₃, I₁) -> (%4i,%4i)\n", _cycle, _prev)

    # Determine filenames for self-energy functions
    fcurr = "dmft1/sigma.bare.$cycle.$curr"
    fprev = "dmft1/sigma.bare.$_cycle.$_prev"

    # Check whether these files are available
    @assert isfile(fcurr) && isfile(fprev)

    # Read in the self-energy functions (previous and current)
    println("Read self-energy functions")
    fcurr, Scurr = read_sigma(ai, fcurr)
    fprev, Sprev = read_sigma(ai, fprev)
    @assert size(Scurr) == size(Sprev) && size(fcurr) == size(fprev)

    # Mix the self-energy functions using linear mixing algorithm
    println("Mix self-energy functions for two successive iterations")
    α = amix(it)
    Snew = Scurr * α + Sprev * (1.0 - α)
    println("  > Mixing parameter α = $α")

    # Write the new self-energy functions into `dmft1/sigma.bare`
    println("Write self-energy functions")
    write_sigma(fcurr, Snew, ai)

    # Check the convergence condition
    println("Evaluate the convergence condition for self-energy functions")
    dist = distance(Scurr, Sprev)
    it.cs = ( dist < get_m("sc") )
    println("  > Averaged ΔΣ = $dist ( convergence is $(it.cs) )" )
end

"""
    mixer_delta(it::IterInfo, ai::Array{Impurity,1})

Try to mix the hybridization functions Δ and then use the mixed values
to update the `dmft1/dmft.delta` file.

See also: [`mixer_core`](@ref).
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
    @printf("  > Curr: (I₃, I₁) -> (%4i,%4i)\n", cycle, curr)
    @printf("  > Prev: (I₃, I₁) -> (%4i,%4i)\n", _cycle, _prev)

    # Determine filenames for hybridization functions
    fcurr = "dmft1/dmft.delta.$cycle.$curr"
    fprev = "dmft1/dmft.delta.$_cycle.$_prev"

    # Check whether these files are available
    @assert isfile(fcurr) && isfile(fprev)

    # Read in the hybridization functions (previous and current)
    println("Read hybridization functions")
    fcurr, Dcurr = read_delta(ai, fcurr)
    fprev, Dprev = read_delta(ai, fprev)
    @assert size(Dcurr) == size(Dprev) && size(fcurr) == size(fprev)

    # Mix the hybridization functions using linear mixing algorithm
    println("Mix hybridization functions for two successive iterations")
    α = amix(it)
    Dnew = Dcurr * α + Dprev * (1.0 - α)
    println("  > Mixing parameter α = $α")

    # Write the new hybridization functions into `dmft1/dmft.delta`
    println("Write hybridization functions")
    write_delta(fcurr, Dnew, ai, "dmft1/dmft.delta")
end

"""
    mixer_eimpx(it::IterInfo, ai::Array{Impurity,1})

Try to mix the local impurity levels εᵢ and then use the mixed value
to update the `dmft1/dmft.eimpx` file.

See also: [`mixer_core`](@ref).
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
    @printf("  > Curr: (I₃, I₁) -> (%4i,%4i)\n", cycle, curr)
    @printf("  > Prev: (I₃, I₁) -> (%4i,%4i)\n", _cycle, _prev)

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
    println("Mix local impurity levels for two successive iterations")
    α = amix(it)
    Enew = Ecurr * α + Eprev * (1.0 - α)
    println("  > Mixing parameter α = $α")

    # Write the new local impurity levels into `dmft1/dmft.eimpx`
    println("Write local impurity levels")
    write_eimpx(Enew, ai, "dmft1/dmft.eimpx")
end

"""
    mixer_gamma(it::IterInfo)

See also: [`mixer_core`](@ref).
"""
function mixer_gamma(it::IterInfo)
    # Print the header
    println("Mixer : Gamma")
end

"""
    amix(it::IterInfo)

Return the mixing factor for mixer component. It should depend on the
current iteration number.

See also: [`IterInfo`](@ref).
"""
function amix(it::IterInfo)
    factor = 1.0
    if it.sc == 1
        factor = exp(-(it.I₁ - 1) * get_m("mixer"))
    else
        factor = exp(-(it.I₃ - 1) * get_m("mixer"))
    end
    return factor
end

"""
    distance(SA::Vector{Array{C64,4}}, SB::Vector{Array{C64,4}})

Calculate the difference between two multi-dimensional arrays.
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
