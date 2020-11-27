"""
    PCASE

Dictionary for parameters: case
"""
PCASE = Dict{String,Any}
        (
            "case"     => ["unknown", "system's name"]
        )

"""
    PDFT

Dictionary for parameters: density functional theory calculations
"""
PDFT  = Dict{String,Any}
        (
            "engine"   => ["unknown", "engine for density functional theory calculations"],
            "smear"    => ["unknown", "scheme for smearing"],
            "kmesh"    => ["unknown", "density of k-mesh sampling in the brillouin zone"],
            "magmom"   => ["unknown", "initial magnetic moment"],
            "lsymm"    => ["unknown", "whether the symmetry is considered"],
            "lspins"   => ["unknown", "whether the spin polarization is considered"],
            "lspinorb" => ["unknown", "whether the spin-orbit coupling is considered"],
            "window"   => ["unknown", "energy window for generating optimal projectors"],
            "loptim"   => ["unknown", "try to optimize the generated projectors"],
            "lproj"    => ["unknown", "try to generate projectors"],
            "nproj"    => ["unknown", "number of types of projectors"],
            "sproj"    => ["unknown", "scheme for generating projectors"]
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
        ()

"""
    PSOLVER

Dictionary for parameters: quantum impurity solvers
"""
PSOLVER = Dict{String,Any}
        ()


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
