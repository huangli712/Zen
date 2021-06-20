#
# Project : Pansy
# Source  : mixer.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/06/21
#

"""
    mixer_sigma(it::IterInfo, ai::Array{Impurity,1})

Try to mix the self-energy functions Σ and then use the mixed values
to update the `dmft1/sigma.bare` file.

See also: [`mixer_core`](@ref).
"""
function mixer_sigma(it::IterInfo, ai::Array{Impurity,1})
    # Print the log
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

    # Mix the self-energy functions
    Snew = Scurr * get_m("mixer") + Sprev * (1.0 - get_m("mixer"))

    # Write the new self-energy functions into `dmft1/sigma.bare`
    println("Write self-energy functions")
    write_sigma(fcurr, Snew, ai)

    # Print blank line for better visualization
    println()
end

"""
    mixer_delta(it::IterInfo, ai::Array{Impurity,1})

Try to mix the hybridization functions Δ and then use the mixed values
to update the `dmft1/dmft.delta` file.

See also: [`mixer_core`](@ref).
"""
function mixer_delta(it::IterInfo, ai::Array{Impurity,1})
    # Print the log
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

    # Mix the hybridization functions
    Dnew = Dcurr * get_m("mixer") + Dprev * (1.0 - get_m("mixer"))

    # Write the new hybridization functions into `dmft1/dmft.delta`
    println("Write hybridization functions")
    write_delta(fcurr, Dnew, ai, "dmft1/dmft.delta")

    # Print blank line for better visualization
    println()
end

"""
    mixer_eimpx(it::IterInfo, ai::Array{Impurity,1})

Try to mix the local impurity levels εᵢ and then use the mixed value
to update the `dmft1/dmft.eimpx` file.

See also: [`mixer_core`](@ref).
"""
function mixer_eimpx(it::IterInfo, ai::Array{Impurity,1})
    # Print the log
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

    # Mix the local impurity levels
    Enew = Ecurr * get_m("mixer") + Eprev * (1.0 - get_m("mixer"))

    # Write the new local impurity levels into `dmft1/dmft.eimpx`
    println("Write local impurity levels")
    write_eimpx(Enew, ai, "dmft1/dmft.eimpx")

    # Print blank line for better visualization
    println()
end

"""
    mixer_gamma(it::IterInfo)

See also: [`mixer_core`](@ref).
"""
function mixer_gamma(it::IterInfo)
    # Print the log
    println("Mixer : Gamma")

    # Print blank line for better visualization
    println()
end
