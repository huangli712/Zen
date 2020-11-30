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
            "mode"     => [missing, UInt  , "scheme of dynamical mean-field theory calculations"],
            "axis"     => [missing, UInt  , "imaginary-time axis or real-frequency axis"],
            "beta"     => [missing, Real  , "inverse system temperature"],
            "niter"    => [missing, UInt  , "number of iterations"],
            "mixer"    => [missing, Real  , "mixing factor"],
            "dcount"   => [missing, String, "scheme of double counting term"],
            "nominal"  => [missing, Array , "nominal impurity occupancy"],
            "cc"       => [missing, Real  , "convergence criterion of charge"],
            "ec"       => [missing, Real  , "convergence criterion of total energy"],
            "fc"       => [missing, Real  , "convergence criterion of force"],
            "lcharge"  => [missing, Bool  , "examine whether charge is converged"],
            "lenergy"  => [missing, Bool  , "examine whether total energy is converged"],
            "lforce"   => [missing, Bool  , "examine whether force is converged"]
        )

"""
    PIMP

Dictionary for parameters: quantum impurity problems
"""
PIMP  = Dict{String,Any}
        (
            "nsite"    => [missing, UInt  , "number of impurity sites"],
            "atoms"    => [missing, Array , "chemical symbols of impurity atoms"],
            "equiv"    => [missing, Array , "equivalency of quantum impurity atoms"],
            "shell"    => [missing, Array , "angular momentum of correlated orbitals"],
            "ising"    => [missing, Array , "interaction types of correlated orbitals"],
            "upara"    => [missing, Array , "Coulomb interaction parameter"],
            "jpara"    => [missing, Array , "Hund's coupling parameter"],
            "lpara"    => [missing, Array , "spin-orbit coupling parameter"]
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
