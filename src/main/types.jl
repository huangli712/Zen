"""
    PCASE

Dictionary for parameters: case summary
"""
PCASE = Dict{String,Any}
        (
            "case"     => [nothing, "system's name"]
        )

"""
    PDFT

Dictionary for parameters: density functional theory calculations
"""
PDFT  = Dict{String,Any}
        (
            "engine"   => [nothing, "engine for density functional theory calculations"],
            "smear"    => [nothing, "scheme for smearing"],
            "kmesh"    => [nothing, "density of k-mesh sampling in the brillouin zone"],
            "magmom"   => [nothing, "initial magnetic moment"],
            "lsymm"    => [nothing, "whether the symmetry is considered"],
            "lspins"   => [nothing, "whether the spin polarization is considered"],
            "lspinorb" => [nothing, "whether the spin-orbit coupling is considered"],
            "window"   => [nothing, "energy window for generating optimal projectors"],
            "loptim"   => [nothing, "try to optimize the generated projectors"],
            "lproj"    => [nothing, "try to generate projectors"],
            "nproj"    => [nothing, "number of types of projectors"],
            "sproj"    => [nothing, "scheme for generating projectors"]
        )

"""
    PDMFT

Dictionary for parameters: dynamical mean-field theory calculations
"""
PDMFT = Dict{String,Any}
        (
            "mode"     => [],
            "axis"     => [],
            "beta"     => [],
            "niter"    => [],
            "mixer"    => [],
            "dcount"   => [],
            "nominal"  => [],
            "cc"       => [],
            "ec"       => [],
            "fc"       => [],
            "lforce"   => [],
            "lcharge"  => [],
            "lenergy"  => []
        )

"""
    PIMP

Dictionary for parameters: quantum impurity problems
"""
PIMP  = Dict{String,Any}
        (
            "nsite"    => [],
            "atoms"    => [],
            "equiv"    => [],
            "shell"    => [],
            "ising"    => [],
            "upara"    => [],
            "jpara"    => [],
            "lpara"    => []
        )

"""
    PSOLVER

Dictionary for parameters: quantum impurity solvers
"""
PSOLVER = Dict{String,Any}
        (
            "engine"   => [],
            "params"   => []
        )


mutable struct IterInfo
    total_iter
    dmft1_iter
    dmft2_iter
    dft_dmft_iter
end

function IterInfo(iter::Int64 = 0)
    IterInfo(iter, iter, iter, iter)
end

mutable struct Lattice
end

mutable struct DFTData
end
