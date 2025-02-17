#
# Project : Lily
# Source  : types.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2024/09/19
#

#=
### *Customized Types*
=#

"Customized types. It is used to define the following dicts."
const DType = Any

"Customized types. It is used to define the following dicts."
const ADT = Array{DType,1}

#=
### *Customized Dictionaries*
=#

#=
*Remarks* :

The values in the following dictionaries are actually arrays, which
usually contain four elements:
* Element[1] -> Actually value.
* Element[2] -> If it is 1, this key-value pair is mandatory.
                If it is 0, this key-value pair is optional.
* Element[3] -> Numerical type (A julia Symbol).
* Element[4] -> Brief explanations.

The following dictionaries are used as global variables.
=#

"""
    PTEST

Dictionary for configuration parameters: general setup.
"""
const PTEST    = Dict{String,ADT}(
    "solver"  => [missing, 1, :String, "Analytic continuation solver"],
    "ptype"   => [missing, 1, :String, "Type of peaks in the spectrum"],
    "ktype"   => [missing, 1, :String, "Type of kernel function"],
    "grid"    => [missing, 1, :String, "Grid for correlation function"],
    "mesh"    => [missing, 1, :String, "Mesh for spectral function"],
    "ngrid"   => [missing, 1, :I64   , "Number of grid points"],
    "nmesh"   => [missing, 1, :I64   , "Number of mesh points"],
    "ntest"   => [missing, 1, :I64   , "Number of tests"],
    "wmax"    => [missing, 1, :F64   , "Right boundary (maximum value) of real mesh"],
    "wmin"    => [missing, 1, :F64   , "Left boundary (minimum value) of real mesh"],
    "pmax"    => [missing, 1, :F64   , "Right boundary (maximum value) for possible peaks"],
    "pmin"    => [missing, 1, :F64   , "Left boundary (minimum value) for possible peaks"],
    "beta"    => [missing, 1, :F64   , "Inverse temperature"],
    "noise"   => [missing, 1, :F64   , "Noise level"],
    "offdiag" => [missing, 1, :Bool  , "Is it the offdiagonal part in matrix-valued function"],
    "lpeak"   => [missing, 1, :Array , "Number of peaks in the spectrum"],
    "pmesh"   => [missing, 0, :Array , "Additional parameters for customizing the mesh"],
)

# Default parameters for PTEST
const _PTEST   = Dict{String,Any}(
    "solver"  => "MaxEnt",
    "ptype"   => "gauss",
    "ktype"   => "fermi",
    "grid"    => "ffreq",
    "mesh"    => "linear",
    "ngrid"   => 10,
    "nmesh"   => 501,
    "ntest"   => 100,
    "wmax"    => 5.0,
    "wmin"    => -5.0,
    "pmax"    => 4.0,
    "pmin"    => -4.0,
    "beta"    => 10.0,
    "noise"   => 1.0e-6,
    "offdiag" => false,
    "lpeak"   => [1,2,3],
)

#=
### *Customized Structs* : *Green's Function Grid*
=#

"""
    AbstractGrid

An abstract type representing the imaginary axis. It is used to build
the internal type system.
"""
abstract type AbstractGrid end

"""
    FermionicImaginaryTimeGrid

Mutable struct. It represents the fermionic imaginary time grid.

### Members
* ntime -> Number of time slices.
* β     -> Inverse temperature.
* τ     -> Vector of grid points， τᵢ.

See also: [`FermionicMatsubaraGrid`](@ref).
"""
mutable struct FermionicImaginaryTimeGrid <: AbstractGrid
    ntime :: I64
    β :: F64
    τ :: Vector{F64}
end

"""
    FermionicMatsubaraGrid

Mutable struct. It represents the fermionic Matsubara frequency grid.

### Members
* nfreq -> Number of Matsubara frequency points.
* β     -> Inverse temperature.
* ω     -> Vector of grid points, iωₙ.

See also: [`FermionicImaginaryTimeGrid`](@ref).
"""
mutable struct FermionicMatsubaraGrid <: AbstractGrid
    nfreq :: I64
    β :: F64
    ω :: Vector{F64}
end

"""
    BosonicImaginaryTimeGrid

Mutable struct. It represents the bosonic imaginary time grid.

### Members
* ntime -> Number of time slices.
* β     -> Inverse temperature.
* τ     -> Vector of grid points, τᵢ.

See also: [`BosonicMatsubaraGrid`](@ref).
"""
mutable struct BosonicImaginaryTimeGrid <: AbstractGrid
    ntime :: I64
    β :: F64
    τ :: Vector{F64}
end

"""
    BosonicMatsubaraGrid

Mutable struct. It represents the bosonic Matsubara frequency grid.

### Members
* nfreq -> Number of Matsubara frequency points.
* β     -> Inverse temperature.
* ω     -> Vector of grid points, iωₙ.

See also: [`BosonicImaginaryTimeGrid`](@ref).
"""
mutable struct BosonicMatsubaraGrid <: AbstractGrid
    nfreq :: I64
    β :: F64
    ω :: Vector{F64}
end

#=
### *Customized Structs* : *Spectral Function Mesh*
=#

"""
    AbstractMesh

An abstract type representing the real axis. It is used to build the
internal type system.
"""
abstract type AbstractMesh end

"""
    LinearMesh

Mutable struct. A linear and uniform mesh.

### Members
* nmesh  -> Number of mesh points.
* wmax   -> Right boundary (maximum value).
* wmin   -> Left boundary (minimum value).
* mesh   -> Mesh itself.
* weight -> Precomputed integration weights (composite trapezoidal rule).

See also: [`TangentMesh`](@ref).
"""
mutable struct LinearMesh <: AbstractMesh
    nmesh :: I64
    wmax :: F64
    wmin :: F64
    mesh :: Vector{F64}
    weight :: Vector{F64}
end

"""
    TangentMesh

Mutable struct. A non-linear and non-uniform mesh. Note that it should
be defined on both negative and positive half-axis.

### Members
* nmesh  -> Number of mesh points.
* wmax   -> Right boundary (maximum value).
* wmin   -> Left boundary (minimum value).
* mesh   -> Mesh itself.
* weight -> Precomputed integration weights (composite trapezoidal rule).

See also: [`LinearMesh`](@ref).
"""
mutable struct TangentMesh <: AbstractMesh
    nmesh :: I64
    wmax :: F64
    wmin :: F64
    mesh :: Vector{F64}
    weight :: Vector{F64}
end

"""
    LorentzMesh

Mutable struct. A non-linear and non-uniform mesh. Note that it should
be defined on both negative and positive half-axis.

### Members
* nmesh  -> Number of mesh points.
* wmax   -> Right boundary (maximum value).
* wmin   -> Left boundary (minimum value).
* mesh   -> Mesh itself.
* weight -> Precomputed integration weights (composite trapezoidal rule).

See also: [`HalfLorentzMesh`](@ref).
"""
mutable struct LorentzMesh <: AbstractMesh
    nmesh :: I64
    wmax :: F64
    wmin :: F64
    mesh :: Vector{F64}
    weight :: Vector{F64}
end

"""
    HalfLorentzMesh

Mutable struct. A non-linear and non-uniform mesh. Note that it should
be defined on positive half-axis only.

### Members
* nmesh  -> Number of mesh points.
* wmax   -> Right boundary (maximum value).
* wmin   -> Left boundary (minimum value). It must be 0.0.
* mesh   -> Mesh itself.
* weight -> Precomputed integration weights (composite trapezoidal rule).

See also: [`LorentzMesh`](@ref).
"""
mutable struct HalfLorentzMesh <: AbstractMesh
    nmesh :: I64
    wmax :: F64
    wmin :: F64
    mesh :: Vector{F64}
    weight :: Vector{F64}
end

"""
    DynamicMesh

Mutable struct. A mesh used internally in the `util/acplot.jl` script. It
supports both uniform and non-uniform meshes.

### Members
* nmesh  -> Number of mesh points.
* wmax   -> Right boundary (maximum value).
* wmin   -> Left boundary (minimum value).
* mesh   -> Mesh itself.
* weight -> Precomputed integration weights (composite trapezoidal rule).

See also: [`LinearMesh`](@ref).
"""
mutable struct DynamicMesh <: AbstractMesh
    nmesh :: I64
    wmax :: F64
    wmin :: F64
    mesh :: Vector{F64}
    weight :: Vector{F64}
end

#=
### *Customized Structs* : *Peak*
=#

"""
    AbstractPeak

An abstract type representing the peaks in the spectral functions. It is
used to build the internal type system.
"""
abstract type AbstractPeak end

"""
    GaussianPeak

Mutable struct. It represents the gaussian peak.

### Members
* A -> Amplitude of the peak.
* Γ -> Broadening of the peak.
* ϵ -> Shift of the peak.
"""
mutable struct GaussianPeak <: AbstractPeak
    A :: F64
    Γ :: F64
    ϵ :: F64
end

"""
    LorentzianPeak

Mutable struct. It represents the lorentzian peak.

### Members
* A -> Amplitude of the peak.
* Γ -> Broadening of the peak.
* ϵ -> Shift of the peak.
"""
mutable struct LorentzianPeak <: AbstractPeak
    A :: F64
    Γ :: F64
    ϵ :: F64
end

"""
    DeltaPeak

Mutable struct. It represents the δ-like peak. Actually, it is a variation
of the GaussianPeak with Γ = 0.02.

### Members
* A -> Amplitude of the peak.
* Γ -> Broadening of the peak.
* ϵ -> Shift of the peak.
"""
mutable struct DeltaPeak <: AbstractPeak
    A :: F64
    Γ :: F64
    ϵ :: F64
end

"""
    RectanglePeak

Mutable struct. It represents the rectangle peak.

### Members
* c -> Center of the peak.
* w -> Width of the peak.
* h -> Height of the peak.
"""
mutable struct RectanglePeak <: AbstractPeak
    c :: F64
    w :: F64
    h :: F64
end

"""
    RiseDecayPeak

Mutable struct. It represents the rise-and-decay peak.

### Members
* c -> Center of the peak.
* γ -> Exponential parameter. It is used to control the width of the peak.
* h -> Height of the peak.
"""
mutable struct RiseDecayPeak <: AbstractPeak
    c :: F64
    γ :: F64
    h :: F64
end

#=
### *Customized Structs* : *Function*
=#

"""
    AbstractFunction

An abstract type representing the general functions. It is used to build
the internal type system.
"""
abstract type AbstractFunction end

"""
    SpectralFunction

Mutable struct. It represents the spectral function A(ω).

### Members
* mesh -> Real frequency mesh, ω.
* image -> Spectral function, A(ω).
"""
mutable struct SpectralFunction <: AbstractFunction
    mesh  :: AbstractMesh
    image :: Vector{F64}
end

"""
    GreenFunction

Mutable struct. It represents the Green's function G(iωₙ) or G(τ).

### Members
* grid -> Imaginary axis grid, ωₙ or τ.
* green -> Green's function, G(iωₙ) or G(τ).
* error -> Standard deviation of the Green's function, σ.
"""
mutable struct GreenFunction <: AbstractFunction
    grid :: AbstractGrid
    green :: Vector{F64}
    error :: Vector{F64}
end
