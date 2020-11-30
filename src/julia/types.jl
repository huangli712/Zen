"""
    PCASE

Dictionary for parameters: case summary    
"""
PCASE = Dict{String,Any}
        (
            "case"     => [missing, String, "system's name"]
        )

"""
    PDFT

Dictionary for parameters: density functional theory calculations
"""
PDFT  = Dict{String,Any}
        (
            "engine"   => [missing, String, "engine for density functional theory calculations"],
            "smear"    => [missing, String, "scheme for smearing"],
            "kmesh"    => [missing, String, "density of k-mesh sampling in the brillouin zone"],
            "magmom"   => [missing, String, "initial magnetic moment"],
            "lsymm"    => [missing, Bool  , "whether the symmetry is considered"],
            "lspins"   => [missing, Bool  , "whether the spin orientations are polarized"],
            "lspinorb" => [missing, Bool  , "whether the spin-orbit coupling is considered"],
            "window"   => [missing, Array , "energy window for generating optimal projectors"],
            "loptim"   => [missing, Bool  , "try to optimize the generated projectors"],
            "lproj"    => [missing, Bool  , "try to generate projectors"],
            "nproj"    => [missing, UInt  , "number of types of projectors"],
            "sproj"    => [missing, Array , "scheme for generating projectors"]
        )

"""
    PDMFT

Dictionary for parameters: dynamical mean-field theory calculations
"""
PDMFT = Dict{String,Any}
        (
            "mode"     => [missing, UInt  , ""],
            "axis"     => [missing, UInt  , ""],
            "beta"     => [missing, Real  , ""],
            "niter"    => [missing, UInt  , ""],
            "mixer"    => [missing, Real  , ""],
            "dcount"   => [missing, String, ""],
            "nominal"  => [missing, Real  , ""],
            "cc"       => [missing, Real  , ""],
            "ec"       => [missing, Real  , ""],
            "fc"       => [missing, Real  , ""],
            "lforce"   => [missing, Bool  , ""],
            "lcharge"  => [missing, Bool  , ""],
            "lenergy"  => [missing, Bool  , ""]
        )

"""
    PIMP

Dictionary for parameters: quantum impurity problems
"""
PIMP  = Dict{String,Any}
        (
            "nsite"    => [missing],
            "atoms"    => [missing],
            "equiv"    => [missing],
            "shell"    => [missing],
            "ising"    => [missing],
            "upara"    => [missing],
            "jpara"    => [missing],
            "lpara"    => [missing]
        )

"""
    PSOLVER

Dictionary for parameters: quantum impurity solvers
"""
PSOLVER = Dict{String,Any}
        (
            "engine"   => [missing, "name of quantum impurity solver"],
            "params"   => [missing, "parameters set of quantum impurity solver"]
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
