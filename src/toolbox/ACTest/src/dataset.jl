#
# Project : Lily
# Source  : dataset.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2025/07/12
#

#=
### *Fermionic Systems: Continuum Spectra*
=#

"""
    STD_FG

Dictionary for standard spectral functions: fermionic + gaussian peaks.
"""
const STD_FG = Dict{String,Any}[
    # Test: 001 / Fermionic + Gaussian Peaks
    # single peak, central
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.8,0.0)
        ],
        "signs" => [1.0]
    ),
    #
    # Test: 002 / Fermionic + Gaussian Peaks
    # single peak, off-centered (left)
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.8,-2.0)
        ],
        "signs" => [1.0]
    ),
    #
    # Test: 003 / Fermionic + Gaussian Peaks
    # single peak, off-centered (right)
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.8,2.0)
        ],
        "signs" => [1.0]
    ),
    #
    # Test: 004 / Fermionic + Gaussian Peaks
    # two peaks, gapless
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.5,-0.8),
            GaussianPeak(1.0,0.5, 0.8)
        ],
        "signs" => [1.0,1.0]
    ),
    #
    # Test: 005 / Fermionic + Gaussian Peaks
    # two peaks, small gap
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.5,-1.8),
            GaussianPeak(1.0,0.5, 1.8)
        ],
        "signs" => [1.0,1.0]
    ),
    #
    # Test: 006 / Fermionic + Gaussian Peaks
    # two peaks, large gap
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.5,-3.0),
            GaussianPeak(1.0,0.5, 3.0)
        ],
        "signs" => [1.0,1.0]
    ),
    #
    # Test: 007 / Fermionic + Gaussian Peaks
    # two peaks, central + off-centered (left)
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.4, 0.0),
            GaussianPeak(1.0,0.4,-3.0)
        ],
        "signs" => [1.0,1.0]
    ),
    #
    # Test: 008 / Fermionic + Gaussian Peaks
    # two peaks, central + off-centered (right)
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.4,0.0),
            GaussianPeak(1.0,0.4,3.0)
        ],
        "signs" => [1.0,1.0]
    ),
    #
    # Test: 009 / Fermionic + Gaussian Peaks
    # two peaks, off-centered (left small + right large)
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(0.5,0.5,-3.0),
            GaussianPeak(1.5,0.5, 3.0)
        ],
        "signs" => [1.0,1.0]
    ),
    #
    # Test: 010 / Fermionic + Gaussian Peaks
    # two peaks, off-centered (left large + right small)
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.5,0.5,-3.0),
            GaussianPeak(0.5,0.5, 3.0)
        ],
        "signs" => [1.0,1.0]
    ),
    #
    # Test: 011 / Fermionic + Gaussian Peaks
    # three peaks, sharp quasiparticle peak + two shoulder peaks
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.1, 0.0),
            GaussianPeak(0.5,1.0,-1.0),
            GaussianPeak(0.5,1.0, 1.0)
        ],
        "signs" => [1.0,1.0,1.0]
    ),
    #
    # Test: 012 / Fermionic + Gaussian Peaks
    # three peaks, sharp quasiparticle peak + two shoulder peaks
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.1, 0.0),
            GaussianPeak(0.5,1.0,-2.0),
            GaussianPeak(0.5,1.0, 2.0)
        ],
        "signs" => [1.0,1.0,1.0]
    ),
    #
    # Test: 013 / Fermionic + Gaussian Peaks
    # three peaks, sharp quasiparticle peak + two Hubbard bands
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.05, 0.0),
            GaussianPeak(0.5,1.00,-3.0),
            GaussianPeak(0.5,1.00, 3.0)
        ],
        "signs" => [1.0,1.0,1.0]
    ),
    #
    # Test: 014 / Fermionic + Gaussian Peaks
    # three peaks, sharp quasiparticle peak + two Hubbard bands
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.05, 0.0),
            GaussianPeak(0.4,0.50,-3.5),
            GaussianPeak(0.4,0.50, 3.5)
        ],
        "signs" => [1.0,1.0,1.0]
    ),
    #
    # Test: 015 / Fermionic + Gaussian Peaks
    # three peaks, large gap
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.05,-1.0),
            GaussianPeak(0.1,0.80,-3.0),
            GaussianPeak(0.4,0.80, 3.0)
        ],
        "signs" => [1.0,1.0,1.0]
    ),
    #
    # Test: 016 / Fermionic + Gaussian Peaks
    # three peaks, large gap
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.05, 1.0),
            GaussianPeak(0.4,0.80,-3.0),
            GaussianPeak(0.1,0.80, 3.0)
        ],
        "signs" => [1.0,1.0,1.0]
    ),
    #
    # Test: 017 / Fermionic + Gaussian Peaks
    # three peaks, gapless, left half-axis
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.2, 0.0),
            GaussianPeak(0.6,0.6,-1.0),
            GaussianPeak(0.6,0.6,-3.0)
        ],
        "signs" => [1.0,1.0,1.0]
    ),
    #
    # Test: 018 / Fermionic + Gaussian Peaks
    # three peaks, gapless, right half-axis
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.2,0.0),
            GaussianPeak(0.6,0.6,1.0),
            GaussianPeak(0.6,0.6,3.0)
        ],
        "signs" => [1.0,1.0,1.0]
    ),
    #
    # Test: 019 / Fermionic + Gaussian Peaks
    # three peaks, gapless, left small + right large
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.2, 0.0),
            GaussianPeak(0.1,0.4,-2.0),
            GaussianPeak(0.6,1.0, 2.0)
        ],
        "signs" => [1.0,1.0,1.0]
    ),
    #
    # Test: 020 / Fermionic + Gaussian Peaks
    # three peaks, gapless, left large + right small
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.2, 0.0),
            GaussianPeak(0.1,0.4, 2.0),
            GaussianPeak(0.6,1.0,-2.0)
        ],
        "signs" => [1.0,1.0,1.0]
    ),
]

#=
### *Fermionic Systems: Discrete Spectra*
=#

"""
    STD_FD

Dictionary for standard spectral functions: fermionic + delta-like peaks.
"""
const STD_FD = Dict{String,Any}[
    # Test: 001 / Fermionic + Delta Peaks
    # single peak, central
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02,0.0)
        ],
        "signs" => [1.0]
    ),
    #
    # Test: 002 / Fermionic + Delta Peaks
    # single peak, off-centered (left)
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02,-2.0)
        ],
        "signs" => [1.0]
    ),
    #
    # Test: 003 / Fermionic + Delta Peaks
    # single peak, off-centered (right)
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02,2.0)
        ],
        "signs" => [1.0]
    ),
    #
    # Test: 004 / Fermionic + Delta Peaks
    # single peak, off-centered (left)
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02,-4.0)
        ],
        "signs" => [1.0]
    ),
    #
    # Test: 005 / Fermionic + Delta Peaks
    # single peak, off-centered (right)
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02,4.0)
        ],
        "signs" => [1.0]
    ),
    #
    # Test: 006 / Fermionic + Delta Peaks
    # two peaks, off-centered (left + right)
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02,-2.0),
            DeltaPeak(1.0,0.02, 2.0)
        ],
        "signs" => [1.0,1.0]
    ),
    #
    # Test: 007 / Fermionic + Delta Peaks
    # two peaks, off-centered (left + left)
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02,-1.0),
            DeltaPeak(1.0,0.02,-3.0)
        ],
        "signs" => [1.0,1.0]
    ),
    #
    # Test: 008 / Fermionic + Delta Peaks
    # two peaks, off-centered (right + right)
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02,1.0),
            DeltaPeak(1.0,0.02,3.0)
        ],
        "signs" => [1.0,1.0]
    ),
    #
    # Test: 009 / Fermionic + Delta Peaks
    # two peaks, off-centered (left near + right far)
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02,-0.5),
            DeltaPeak(1.0,0.02, 3.0)
        ],
        "signs" => [1.0,1.0]
    ),
    #
    # Test: 010 / Fermionic + Delta Peaks
    # two peaks, off-centered (left far + right near)
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02,-3.0),
            DeltaPeak(1.0,0.02, 0.5)
        ],
        "signs" => [1.0,1.0]
    ),
    #
    # Test: 011 / Fermionic + Delta Peaks
    # three peaks, small distance
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02, 0.0),
            DeltaPeak(1.0,0.02,-0.5),
            DeltaPeak(1.0,0.02, 0.5)
        ],
        "signs" => [1.0,1.0,1.0]
    ),
    #
    # Test: 012 / Fermionic + Delta Peaks
    # three peaks, large distance
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02, 0.0),
            DeltaPeak(1.0,0.02,-3.0),
            DeltaPeak(1.0,0.02, 3.0)
        ],
        "signs" => [1.0,1.0,1.0]
    ),
    #
    # Test: 013 / Fermionic + Delta Peaks
    # four peaks, small distance
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02,-1.5),
            DeltaPeak(1.0,0.02,-0.5),
            DeltaPeak(1.0,0.02, 0.5),
            DeltaPeak(1.0,0.02, 1.5)
        ],
        "signs" => [1.0,1.0,1.0,1.0]
    ),
    #
    # Test: 014 / Fermionic + Delta Peaks
    # four peaks, large distance
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02,-3.0),
            DeltaPeak(1.0,0.02,-0.5),
            DeltaPeak(1.0,0.02, 0.5),
            DeltaPeak(1.0,0.02, 3.0)
        ],
        "signs" => [1.0,1.0,1.0,1.0]
    ),
    #
    # Test: 015 / Fermionic + Delta Peaks
    # four peaks, large distance
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02,-3.0),
            DeltaPeak(1.0,0.02,-1.0),
            DeltaPeak(1.0,0.02, 1.0),
            DeltaPeak(1.0,0.02, 3.0)
        ],
        "signs" => [1.0,1.0,1.0,1.0]
    ),
    #
    # Test: 016 / Fermionic + Delta Peaks
    # five peaks, small distance
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02, 0.0),
            DeltaPeak(1.0,0.02,-3.0),
            DeltaPeak(1.0,0.02,-2.5),
            DeltaPeak(1.0,0.02, 2.5),
            DeltaPeak(1.0,0.02, 3.0)
        ],
        "signs" => [1.0,1.0,1.0,1.0,1.0]
    ),
    #
    # Test: 017 / Fermionic + Delta Peaks
    # five peaks, small distance
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02, 0.0),
            DeltaPeak(1.0,0.02,-3.0),
            DeltaPeak(1.0,0.02,-0.5),
            DeltaPeak(1.0,0.02, 0.5),
            DeltaPeak(1.0,0.02, 3.0)
        ],
        "signs" => [1.0,1.0,1.0,1.0,1.0]
    ),
    #
    # Test: 018 / Fermionic + Delta Peaks
    # five peaks, large distance
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02, 0.0),
            DeltaPeak(1.0,0.02,-4.0),
            DeltaPeak(1.0,0.02,-2.0),
            DeltaPeak(1.0,0.02, 2.0),
            DeltaPeak(1.0,0.02, 4.0)
        ],
        "signs" => [1.0,1.0,1.0,1.0,1.0]
    ),
    #
    # Test: 019 / Fermionic + Delta Peaks
    # six peaks, small distance
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02,-4.0),
            DeltaPeak(1.0,0.02,-3.5),
            DeltaPeak(1.0,0.02,-3.0),
            DeltaPeak(1.0,0.02, 3.0),
            DeltaPeak(1.0,0.02, 3.5),
            DeltaPeak(1.0,0.02, 4.0)
        ],
        "signs" => [1.0,1.0,1.0,1.0,1.0,1.0]
    ),
    #
    # Test: 020 / Fermionic + Delta Peaks
    # three peaks, small distance
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02,-4.0),
            DeltaPeak(1.0,0.02,-1.0),
            DeltaPeak(1.0,0.02,-0.5),
            DeltaPeak(1.0,0.02, 0.5),
            DeltaPeak(1.0,0.02, 1.0),
            DeltaPeak(1.0,0.02, 4.0)
        ],
        "signs" => [1.0,1.0,1.0,1.0,1.0,1.0]
    ),
]

#=
### *Fermionic Systems: Non-positive Definite Spectra*
=#

"""
    STD_FRD

Dictionary for standard spectral functions: fermionic + rise-and-decay peaks.
"""
const STD_FRD = Dict{String,Any}[
    # Test: 001 / Fermionic + Rise-And-Decay Peaks
    # single peak, central, down
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => true,
        "peaks" => [
            RiseDecayPeak(0.0,1.5,1.0)
        ],
        "signs" => [-1.0]
    ),
    #
    # Test: 002 / Fermionic + Rise-And-Decay Peaks
    # two peaks, left down + right up
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => true,
        "peaks" => [
            RiseDecayPeak(-1.0,1.5,0.5),
            RiseDecayPeak( 1.0,1.5,0.5)
        ],
        "signs" => [-1.0,1.0]
    ),
    #
    # Test: 003 / Fermionic + Rise-And-Decay Peaks
    # two peaks, left up + right down
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => true,
        "peaks" => [
            RiseDecayPeak(-1.0,1.5,0.5),
            RiseDecayPeak( 1.0,1.5,0.5)
        ],
        "signs" => [1.0,-1.0]
    ),
    #
    # Test: 004 / Fermionic + Rise-And-Decay Peaks
    # two peaks, left down small + right up large
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => true,
        "peaks" => [
            RiseDecayPeak(-3.0,2.0,0.2),
            RiseDecayPeak( 3.0,2.0,0.5)
        ],
        "signs" => [-1.0,1.0]
    ),
    #
    # Test: 005 / Fermionic + Rise-And-Decay Peaks
    # two peaks, left up large + right down small
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => true,
        "peaks" => [
            RiseDecayPeak(-3.0,2.0,0.5),
            RiseDecayPeak( 3.0,2.0,0.2)
        ],
        "signs" => [1.0,-1.0]
    ),
    #
    # Test: 006 / Fermionic + Rise-And-Decay Peaks
    # two peaks, left down + right down
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => true,
        "peaks" => [
            RiseDecayPeak(-2.0,1.5,0.5),
            RiseDecayPeak( 2.0,1.5,0.5)
        ],
        "signs" => [-1.0,-1.0]
    ),
    #
    # Test: 007 / Fermionic + Rise-And-Decay Peaks
    # three peaks, left down + central up + right down
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => true,
        "peaks" => [
            RiseDecayPeak(-2.5,1.5,0.2),
            RiseDecayPeak( 0.0,1.5,0.5),
            RiseDecayPeak( 2.5,1.5,0.2)
        ],
        "signs" => [-1.0,1.0,-1.0]
    ),
    #
    # Test: 008 / Fermionic + Rise-And-Decay Peaks
    # three peaks, left up + central down + right up
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => true,
        "peaks" => [
            RiseDecayPeak(-2.5,1.5,0.5),
            RiseDecayPeak( 0.0,1.5,0.2),
            RiseDecayPeak( 2.5,1.5,0.5)
        ],
        "signs" => [1.0,-1.0,1.0]
    ),
    #
    # Test: 009 / Fermionic + Rise-And-Decay Peaks
    # four peaks, left down + left up + right up + right down
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => true,
        "peaks" => [
            RiseDecayPeak(-3.0,0.5,0.20),
            RiseDecayPeak(-1.0,0.5,0.25),
            RiseDecayPeak( 1.0,0.5,0.25),
            RiseDecayPeak( 3.0,0.5,0.20)
        ],
        "signs" => [-1.0,1.0,1.0,-1.0]
    ),
    #
    # Test: 010 / Fermionic + Rise-And-Decay Peaks
    # five peaks, left up + left down + central up + right down + right up
    Dict(
        "ktype" => "fermi",
        "grid"  => "ffreq",
        "mesh"  => "linear",
        "fnpd"  => true,
        "peaks" => [
            RiseDecayPeak(-3.0,0.5,0.2),
            RiseDecayPeak(-1.0,0.5,0.5),
            RiseDecayPeak( 0.0,0.5,0.5),
            RiseDecayPeak( 1.0,0.5,0.5),
            RiseDecayPeak( 3.0,0.5,0.2)
        ],
        "signs" => [1.0,-1.0,1.0,-1.0,1.0]
    ),
]

#=
### *Bosonic Systems: Continuum Spectra*
=#

"""
    STD_BG

Dictionary for standard spectral functions: bosonic + gaussian peaks.
"""
const STD_BG = Dict{String,Any}[
    # Test: 001 / Bosonic + Gaussian Peaks
    # single peak, central
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.8,0.0)
        ],
        "signs" => [1.0]
    ),
    #
    # Test: 002 / Bosonic + Gaussian Peaks
    # single peak, off-centered (left)
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.8,-2.0)
        ],
        "signs" => [1.0]
    ),
    #
    # Test: 003 / Bosonic + Gaussian Peaks
    # single peak, off-centered (right)
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.8,2.0)
        ],
        "signs" => [1.0]
    ),
    #
    # Test: 004 / Bosonic + Gaussian Peaks
    # two peaks, gapless
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.5,-0.8),
            GaussianPeak(1.0,0.5, 0.8)
        ],
        "signs" => [1.0,1.0]
    ),
    #
    # Test: 005 / Bosonic + Gaussian Peaks
    # two peaks, small gap
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.5,-1.8),
            GaussianPeak(1.0,0.5, 1.8)
        ],
        "signs" => [1.0,1.0]
    ),
    #
    # Test: 006 / Bosonic + Gaussian Peaks
    # two peaks, large gap
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.5,-3.0),
            GaussianPeak(1.0,0.5, 3.0)
        ],
        "signs" => [1.0,1.0]
    ),
    #
    # Test: 007 / Bosonic + Gaussian Peaks
    # two peaks, central + off-centered (left)
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.4, 0.0),
            GaussianPeak(1.0,0.4,-3.0)
        ],
        "signs" => [1.0,1.0]
    ),
    #
    # Test: 008 / Bosonic + Gaussian Peaks
    # two peaks, central + off-centered (right)
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.4,0.0),
            GaussianPeak(1.0,0.4,3.0)
        ],
        "signs" => [1.0,1.0]
    ),
    #
    # Test: 009 / Bosonic + Gaussian Peaks
    # two peaks, off-centered (left small + right large)
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(0.5,0.5,-3.0),
            GaussianPeak(1.5,0.5, 3.0)
        ],
        "signs" => [1.0,1.0]
    ),
    #
    # Test: 010 / Bosonic + Gaussian Peaks
    # two peaks, off-centered (left large + right small)
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.5,0.5,-3.0),
            GaussianPeak(0.5,0.5, 3.0)
        ],
        "signs" => [1.0,1.0]
    ),
    #
    # Test: 011 / Bosonic + Gaussian Peaks
    # three peaks, sharp quasiparticle peak + two shoulder peaks
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.1, 0.0),
            GaussianPeak(0.5,1.0,-1.0),
            GaussianPeak(0.5,1.0, 1.0)
        ],
        "signs" => [1.0,1.0,1.0]
    ),
    #
    # Test: 012 / Bosonic + Gaussian Peaks
    # three peaks, sharp quasiparticle peak + two shoulder peaks
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.1, 0.0),
            GaussianPeak(0.5,1.0,-2.0),
            GaussianPeak(0.5,1.0, 2.0)
        ],
        "signs" => [1.0,1.0,1.0]
    ),
    #
    # Test: 013 / Bosonic + Gaussian Peaks
    # three peaks, sharp quasiparticle peak + two Hubbard bands
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.05, 0.0),
            GaussianPeak(0.5,1.00,-3.0),
            GaussianPeak(0.5,1.00, 3.0)
        ],
        "signs" => [1.0,1.0,1.0]
    ),
    #
    # Test: 014 / Bosonic + Gaussian Peaks
    # three peaks, sharp quasiparticle peak + two Hubbard bands
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.05, 0.0),
            GaussianPeak(0.4,0.50,-3.5),
            GaussianPeak(0.4,0.50, 3.5)
        ],
        "signs" => [1.0,1.0,1.0]
    ),
    #
    # Test: 015 / Bosonic + Gaussian Peaks
    # three peaks, large gap
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.05,-1.0),
            GaussianPeak(0.1,0.80,-3.0),
            GaussianPeak(0.4,0.80, 3.0)
        ],
        "signs" => [1.0,1.0,1.0]
    ),
    #
    # Test: 016 / Bosonic + Gaussian Peaks
    # three peaks, large gap
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.05, 1.0),
            GaussianPeak(0.4,0.80,-3.0),
            GaussianPeak(0.1,0.80, 3.0)
        ],
        "signs" => [1.0,1.0,1.0]
    ),
    #
    # Test: 017 / Bosonic + Gaussian Peaks
    # three peaks, gapless, left half-axis
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.2, 0.0),
            GaussianPeak(0.6,0.6,-1.0),
            GaussianPeak(0.6,0.6,-3.0)
        ],
        "signs" => [1.0,1.0,1.0]
    ),
    #
    # Test: 018 / Bosonic + Gaussian Peaks
    # three peaks, gapless, right half-axis
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.2,0.0),
            GaussianPeak(0.6,0.6,1.0),
            GaussianPeak(0.6,0.6,3.0)
        ],
        "signs" => [1.0,1.0,1.0]
    ),
    #
    # Test: 019 / Bosonic + Gaussian Peaks
    # three peaks, gapless, left small + right large
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.2, 0.0),
            GaussianPeak(0.1,0.4,-2.0),
            GaussianPeak(0.6,1.0, 2.0)
        ],
        "signs" => [1.0,1.0,1.0]
    ),
    #
    # Test: 020 / Bosonic + Gaussian Peaks
    # three peaks, gapless, left large + right small
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            GaussianPeak(1.0,0.2, 0.0),
            GaussianPeak(0.1,0.4, 2.0),
            GaussianPeak(0.6,1.0,-2.0)
        ],
        "signs" => [1.0,1.0,1.0]
    ),
]

#=
### *Bosonic Systems: Discrete Spectra*
=#

"""
    STD_BD

Dictionary for standard spectral functions: bosonic + delta-like peaks.
"""
const STD_BD = Dict{String,Any}[
    # Test: 001 / Bosonic + Delta Peaks
    # single peak, off-centered (left)
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02,-0.2)
        ],
        "signs" => [1.0]
    ),
    #
    # Test: 002 / Bosonic + Delta Peaks
    # single peak, off-centered (left)
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02,-2.0)
        ],
        "signs" => [1.0]
    ),
    #
    # Test: 003 / Bosonic + Delta Peaks
    # single peak, off-centered (right)
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02,2.0)
        ],
        "signs" => [1.0]
    ),
    #
    # Test: 004 / Bosonic + Delta Peaks
    # single peak, off-centered (left)
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02,-4.0)
        ],
        "signs" => [1.0]
    ),
    #
    # Test: 005 / Bosonic + Delta Peaks
    # single peak, off-centered (right)
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02,4.0)
        ],
        "signs" => [1.0]
    ),
    #
    # Test: 006 / Bosonic + Delta Peaks
    # two peaks, off-centered (left + right)
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02,-2.0),
            DeltaPeak(1.0,0.02, 2.0)
        ],
        "signs" => [1.0,1.0]
    ),
    #
    # Test: 007 / Bosonic + Delta Peaks
    # two peaks, off-centered (left + left)
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02,-1.0),
            DeltaPeak(1.0,0.02,-3.0)
        ],
        "signs" => [1.0,1.0]
    ),
    #
    # Test: 008 / Bosonic + Delta Peaks
    # two peaks, off-centered (right + right)
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02,1.0),
            DeltaPeak(1.0,0.02,3.0)
        ],
        "signs" => [1.0,1.0]
    ),
    #
    # Test: 009 / Bosonic + Delta Peaks
    # two peaks, off-centered (left near + right far)
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02,-0.5),
            DeltaPeak(1.0,0.02, 3.0)
        ],
        "signs" => [1.0,1.0]
    ),
    #
    # Test: 010 / Bosonic + Delta Peaks
    # two peaks, off-centered (left far + right near)
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02,-3.0),
            DeltaPeak(1.0,0.02, 0.5)
        ],
        "signs" => [1.0,1.0]
    ),
    #
    # Test: 011 / Bosonic + Delta Peaks
    # three peaks, small distance
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02, 0.0),
            DeltaPeak(1.0,0.02,-0.5),
            DeltaPeak(1.0,0.02, 0.5)
        ],
        "signs" => [1.0,1.0,1.0]
    ),
    #
    # Test: 012 / Bosonic + Delta Peaks
    # three peaks, large distance
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02, 0.0),
            DeltaPeak(1.0,0.02,-3.0),
            DeltaPeak(1.0,0.02, 3.0)
        ],
        "signs" => [1.0,1.0,1.0]
    ),
    #
    # Test: 013 / Bosonic + Delta Peaks
    # four peaks, small distance
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02,-1.5),
            DeltaPeak(1.0,0.02,-0.5),
            DeltaPeak(1.0,0.02, 0.5),
            DeltaPeak(1.0,0.02, 1.5)
        ],
        "signs" => [1.0,1.0,1.0,1.0]
    ),
    #
    # Test: 014 / Bosonic + Delta Peaks
    # four peaks, large distance
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02,-3.0),
            DeltaPeak(1.0,0.02,-0.5),
            DeltaPeak(1.0,0.02, 0.5),
            DeltaPeak(1.0,0.02, 3.0)
        ],
        "signs" => [1.0,1.0,1.0,1.0]
    ),
    #
    # Test: 015 / Bosonic + Delta Peaks
    # four peaks, large distance
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02,-3.0),
            DeltaPeak(1.0,0.02,-1.0),
            DeltaPeak(1.0,0.02, 1.0),
            DeltaPeak(1.0,0.02, 3.0)
        ],
        "signs" => [1.0,1.0,1.0,1.0]
    ),
    #
    # Test: 016 / Bosonic + Delta Peaks
    # five peaks, small distance
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02, 0.0),
            DeltaPeak(1.0,0.02,-3.0),
            DeltaPeak(1.0,0.02,-2.5),
            DeltaPeak(1.0,0.02, 2.5),
            DeltaPeak(1.0,0.02, 3.0)
        ],
        "signs" => [1.0,1.0,1.0,1.0,1.0]
    ),
    #
    # Test: 017 / Bosonic + Delta Peaks
    # five peaks, small distance
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02, 0.0),
            DeltaPeak(1.0,0.02,-3.0),
            DeltaPeak(1.0,0.02,-0.5),
            DeltaPeak(1.0,0.02, 0.5),
            DeltaPeak(1.0,0.02, 3.0)
        ],
        "signs" => [1.0,1.0,1.0,1.0,1.0]
    ),
    #
    # Test: 018 / Bosonic + Delta Peaks
    # five peaks, large distance
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02, 0.0),
            DeltaPeak(1.0,0.02,-4.0),
            DeltaPeak(1.0,0.02,-2.0),
            DeltaPeak(1.0,0.02, 2.0),
            DeltaPeak(1.0,0.02, 4.0)
        ],
        "signs" => [1.0,1.0,1.0,1.0,1.0]
    ),
    #
    # Test: 019 / Bosonic + Delta Peaks
    # six peaks, small distance
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02,-4.0),
            DeltaPeak(1.0,0.02,-3.5),
            DeltaPeak(1.0,0.02,-3.0),
            DeltaPeak(1.0,0.02, 3.0),
            DeltaPeak(1.0,0.02, 3.5),
            DeltaPeak(1.0,0.02, 4.0)
        ],
        "signs" => [1.0,1.0,1.0,1.0,1.0,1.0]
    ),
    #
    # Test: 020 / Bosonic + Delta Peaks
    # three peaks, small distance
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => false,
        "peaks" => [
            DeltaPeak(1.0,0.02,-4.0),
            DeltaPeak(1.0,0.02,-1.0),
            DeltaPeak(1.0,0.02,-0.5),
            DeltaPeak(1.0,0.02, 0.5),
            DeltaPeak(1.0,0.02, 1.0),
            DeltaPeak(1.0,0.02, 4.0)
        ],
        "signs" => [1.0,1.0,1.0,1.0,1.0,1.0]
    ),
]

#=
### *Bosonic Systems: Non-positive Definite Spectra*
=#

"""
    STD_BRD

Dictionary for standard spectral functions: bosonic + rise-and-decay peaks.
"""
const STD_BRD = Dict{String,Any}[
    # Test: 001 / Bosonic + Rise-And-Decay Peaks
    # single peak, central
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => true,
        "peaks" => [
            RiseDecayPeak(0.0,1.5,1.0)
        ],
        "signs" => [-1.0]
    ),
    #
    # Test: 002 / Bosonic + Rise-And-Decay Peaks
    # two peaks, left down + right up
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => true,
        "peaks" => [
            RiseDecayPeak(-0.5,1.5,0.5),
            RiseDecayPeak( 1.0,1.5,0.5)
        ],
        "signs" => [-1.0,1.0]
    ),
    #
    # Test: 003 / Bosonic + Rise-And-Decay Peaks
    # two peaks, left up + right down
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => true,
        "peaks" => [
            RiseDecayPeak(-1.0,1.5,0.5),
            RiseDecayPeak( 1.0,1.5,0.5)
        ],
        "signs" => [1.0,-1.0]
    ),
    #
    # Test: 004 / Bosonic + Rise-And-Decay Peaks
    # two peaks, left up small + right down large
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => true,
        "peaks" => [
            RiseDecayPeak(-3.0,2.0,0.2),
            RiseDecayPeak( 3.0,2.0,0.5)
        ],
        "signs" => [1.0,-1.0]
    ),
    #
    # Test: 005 / Bosonic + Rise-And-Decay Peaks
    # two peaks, left up large + right down small
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => true,
        "peaks" => [
            RiseDecayPeak(-3.0,2.0,0.5),
            RiseDecayPeak( 3.0,2.0,0.2)
        ],
        "signs" => [1.0,-1.0]
    ),
    #
    # Test: 006 / Bosonic + Rise-And-Decay Peaks
    # two peaks, left down + right down
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => true,
        "peaks" => [
            RiseDecayPeak(-2.0,1.5,0.5),
            RiseDecayPeak( 2.0,1.5,0.5)
        ],
        "signs" => [-1.0,-1.0]
    ),
    #
    # Test: 007 / Bosonic + Rise-And-Decay Peaks
    # three peaks, left down + central up + right down
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => true,
        "peaks" => [
            RiseDecayPeak(-2.5,1.5,0.2),
            RiseDecayPeak( 0.0,1.5,0.5),
            RiseDecayPeak( 2.5,1.5,0.2)
        ],
        "signs" => [-1.0,1.0,-1.0]
    ),
    #
    # Test: 008 / Bosonic + Rise-And-Decay Peaks
    # three peaks, left up + central down + right up
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => true,
        "peaks" => [
            RiseDecayPeak(-2.5,1.5,0.5),
            RiseDecayPeak( 0.0,1.5,0.2),
            RiseDecayPeak( 2.5,1.5,0.5)
        ],
        "signs" => [1.0,-1.0,1.0]
    ),
    #
    # Test: 009 / Bosonic + Rise-And-Decay Peaks
    # four peaks, left down + left up + right up + right down
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => true,
        "peaks" => [
            RiseDecayPeak(-3.0,0.5,0.20),
            RiseDecayPeak(-1.0,0.5,0.25),
            RiseDecayPeak( 1.0,0.5,0.25),
            RiseDecayPeak( 3.0,0.5,0.20)
        ],
        "signs" => [-1.0,1.0,1.0,-1.0]
    ),
    #
    # Test: 010 / Bosonic + Rise-And-Decay Peaks
    # five peaks, left up + left down + central up + right down + right up
    Dict(
        "ktype" => "boson",
        "grid"  => "bfreq",
        "mesh"  => "linear",
        "fnpd"  => true,
        "peaks" => [
            RiseDecayPeak(-3.0,0.5,0.2),
            RiseDecayPeak(-1.0,0.5,0.5),
            RiseDecayPeak( 0.0,0.5,0.5),
            RiseDecayPeak( 1.0,0.5,0.5),
            RiseDecayPeak( 3.0,0.5,0.2)
        ],
        "signs" => [1.0,-1.0,1.0,-1.0,1.0]
    ),
]
