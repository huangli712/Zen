"""
    PCASE

Dictionary for parameters: case summary
"""
PCASE = Dict{String,Any}
        (
            "case"     => [missing, "system's name"]
        )

"""
    PDFT

Dictionary for parameters: density functional theory calculations
"""
PDFT  = Dict{String,Any}
        (
            "engine"   => [missing, "engine for density functional theory calculations"],
            "smear"    => [missing, "scheme for smearing"],
            "kmesh"    => [missing, "density of k-mesh sampling in the brillouin zone"],
            "magmom"   => [missing, "initial magnetic moment"],
            "lsymm"    => [missing, "whether the symmetry is considered"],
            "lspins"   => [missing, "whether the spin polarization is considered"],
            "lspinorb" => [missing, "whether the spin-orbit coupling is considered"],
            "window"   => [missing, "energy window for generating optimal projectors"],
            "loptim"   => [missing, "try to optimize the generated projectors"],
            "lproj"    => [missing, "try to generate projectors"],
            "nproj"    => [missing, "number of types of projectors"],
            "sproj"    => [missing, "scheme for generating projectors"]
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
