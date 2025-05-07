#
# Project : Lily
# Source  : inout.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2025/05/07
#

"""
    write_spectrum(am::AbstractMesh, Aout::Vector{F64})

Write spectrum A(ω) to `image.data`. The grid is defined in `am`, and
the spectral data are contained in `Aout`.

### Arguments
* am   -> Real frequency mesh, ω.
* Aout -> Spectral function, A(ω).

### Returns
N/A
"""
function write_spectrum(am::AbstractMesh, Aout::Vector{F64})
    @assert length(am) == length(Aout)

    open("image.data", "w") do fout
        for i in eachindex(am)
            @printf(fout, "%16.12f %16.12f\n", am[i], Aout[i])
        end
    end
end

"""
    write_spectrum(ind::I64, sf::SpectralFunction)

Write spectrum A(ω) to `image.data.i`. All information about the spectral
function is included in `sf`.

### Arguments
* ind -> Index for the spectral function.
* sf -> A SpectralFunction struct.

### Returns
N/A

See also: [`SpectralFunction`](@ref).
"""
function write_spectrum(ind::I64, sf::SpectralFunction)
    @assert ind ≥ 1

    open("image.data." * string(ind), "w") do fout
        for i in eachindex(sf.mesh)
            @printf(fout, "%16.12f %16.12f\n", sf.mesh[i], sf.image[i])
        end
    end
end

"""
    write_backward(ag::AbstractGrid, G::Vector{F64})

We can use the calculated spectrum in real axis to generate the Green's
function data in imaginary axis. This function will write the data to
`green.data`, which can be fed into the analytic continuation tools.
Here, `G` is the constructed Green's function data.

### Arguments
* ag -> Grid for input data.
* G  -> Constructed Green's function.

### Returns
N/A

See also: [`reprod`](@ref).
"""
function write_backward(ag::AbstractGrid, G::Vector{F64})
    ngrid = length(ag)
    ng = length(G)
    @assert ngrid == ng || ngrid * 2 == ng

    # The reproduced data are defined in imaginary time axis.
    if ngrid == ng
        open("green.data", "w") do fout
            for i in eachindex(ag)
                @printf(fout, "%16.12f %16.12f\n", ag[i], G[i])
            end
        end
    # The reproduced data are defined in Matsubara frequency axis.
    else
        open("green.data", "w") do fout
            for i in eachindex(ag)
                @printf(fout, "%16.12f %16.12f %16.12f\n", ag[i], G[i], G[i+ngrid])
            end
        end
    end
end

"""
    write_backward(ind::I64, gf::GreenFunction)

Write the Green's function data to `green.data.i`. All information about
the Green's function is included in `gf`.

### Arguments
* ind -> Index for the Green's function.
* gf -> A GreenFunction struct.

### Returns
N/A

See also: [`reprod`](@ref).
"""
function write_backward(ind::I64, gf::GreenFunction)
    @assert ind ≥ 1

    ag = gf.grid
    G = gf.green
    err = gf.error

    ngrid = length(ag)
    ng = length(G)
    @assert ngrid == ng || ngrid * 2 == ng

    # The reproduced data are defined in imaginary time axis.
    if ngrid == ng
        open("green.data." * string(ind), "w") do fout
            for i in eachindex(ag)
                @printf(fout, "%16.12f ", ag[i])
                @printf(fout, "%16.12f %16.12f\n", G[i], err[i])
            end
        end
    # The reproduced data are defined in Matsubara frequency axis.
    else
        open("green.data." * string(ind), "w") do fout
            for i in eachindex(ag)
                @printf(fout, "%16.12f ", ag[i])
                @printf(fout, "%16.12f %16.12f ", G[i], G[i+ngrid])
                @printf(fout, "%16.12f %16.12f\n", err[i], err[i+ngrid])
            end
        end
    end
end

"""
    write_backward(ind::I64, gf::Vector{GreenFunction})

Write the Green's function data to `green.bin.i`. All information about
the Green's function is included in `gf`. Note that `gf` contains multiple
data bins. The number of data bins is `length(gf)`. A GreenFunction struct
is related to a data bin.

### Arguments
* ind -> Index for the Green's function.
* gf -> A vector of GreenFunction struct.

### Returns
N/A

See also: [`reprod`](@ref).
"""
function write_backward(ind::I64, gf::Vector{GreenFunction})
    @assert ind ≥ 1
    nbins = length(gf)
    @assert nbins > 1

    # The reproduced data must be defined in imaginary time axis.
    open("green.bin." * string(ind), "w") do fout
        for i = 1:nbins
            # Unpack the GreenFunction struct
            ag = gf[i].grid
            G = gf[i].green
            err = gf[i].error
            #
            # Check the dimensional parameters
            ngrid = length(ag)
            ng = length(G)
            @assert ngrid == ng
            #
            # Write data bin
            @printf(fout, "# data bin : %6i\n", i)
            for i in eachindex(ag)
                @printf(fout, "%16.12f ", ag[i])
                @printf(fout, "%16.12f %16.12f\n", G[i], err[i])
            end
            println(fout)
            println(fout)
        end
    end
end

"""
    Base.show(io::IO, 𝑝::GaussianPeak)

Write a GaussianPeak struct.

See also: [`GaussianPeak`](@ref).
"""
function Base.show(io::IO, 𝑝::GaussianPeak)
    println("peak type : gaussian")
    @printf("  amplitude   : %16.12f\n", 𝑝.A)
    @printf("  broadening  : %16.12f\n", 𝑝.Γ)
    @printf("  shift       : %16.12f  ", 𝑝.ϵ)
end

"""
    Base.show(io::IO, 𝑝::LorentzianPeak)

Write a LorentzianPeak struct.

See also: [`LorentzianPeak`](@ref).
"""
function Base.show(io::IO, 𝑝::LorentzianPeak)
    println("peak type : lorentzian")
    @printf("  amplitude   : %16.12f\n", 𝑝.A)
    @printf("  broadening  : %16.12f\n", 𝑝.Γ)
    @printf("  shift       : %16.12f  ", 𝑝.ϵ)
end

"""
    Base.show(io::IO, 𝑝::DeltaPeak)

Write a DeltaPeak struct.

See also: [`DeltaPeak`](@ref).
"""
function Base.show(io::IO, 𝑝::DeltaPeak)
    println("peak type : delta")
    @printf("  amplitude   : %16.12f\n", 𝑝.A)
    @printf("  broadening  : %16.12f\n", 𝑝.Γ)
    @printf("  shift       : %16.12f  ", 𝑝.ϵ)
end

"""
    Base.show(io::IO, 𝑝::RectanglePeak)

Write a RectanglePeak struct.

See also: [`RectanglePeak`](@ref).
"""
function Base.show(io::IO, 𝑝::RectanglePeak)
    println("peak type : rectangle")
    @printf("  center : %16.12f\n", 𝑝.c)
    @printf("  width  : %16.12f\n", 𝑝.w)
    @printf("  height : %16.12f  ", 𝑝.h)
end

"""
    Base.show(io::IO, 𝑝::RiseDecayPeak)

Write a RiseDecayPeak struct.

See also: [`RiseDecayPeak`](@ref).
"""
function Base.show(io::IO, 𝑝::RiseDecayPeak)
    println("peak type : risedecay")
    @printf("  center : %16.12f\n", 𝑝.c)
    @printf("  expon. : %16.12f\n", 𝑝.γ)
    @printf("  height : %16.12f  ", 𝑝.h)
end
