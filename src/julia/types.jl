"""
    PCASE

Dictionary for parameters: case summary    
"""
PCASE = Dict{String,Any}(
            "case"     => [missing, 0, String, "system's name"]
        )

"""
    PDFT

Dictionary for parameters: density functional theory calculations
"""
PDFT  = Dict{String,Any}(
            "engine"   => [missing, 0, String, "engine for density functional theory calculations"],
            "smear"    => [missing, 0, String, "scheme for smearing"],
            "kmesh"    => [missing, 0, String, "density of k-mesh sampling in the brillouin zone"],
            "magmom"   => [missing, 0, String, "initial magnetic moment"],
            "lsymm"    => [missing, 0, Bool  , "whether the symmetry is considered"],
            "lspins"   => [missing, 0, Bool  , "whether the spin orientations are polarized"],
            "lspinorb" => [missing, 0, Bool  , "whether the spin-orbit coupling is considered"],
            "window"   => [missing, 0, Array , "energy window for generating optimal projectors"],
            "loptim"   => [missing, 0, Bool  , "try to optimize the generated projectors"],
            "lproj"    => [missing, 0, Bool  , "try to generate projectors"],
            "nproj"    => [missing, 0, UInt  , "number of types of projectors"],
            "sproj"    => [missing, 0, Array , "scheme for generating projectors"]
        )

"""
    PDMFT

Dictionary for parameters: dynamical mean-field theory calculations
"""
PDMFT = Dict{String,Any}(
            "mode"     => [missing, 0, UInt  , "scheme of dynamical mean-field theory calculations"],
            "axis"     => [missing, 0, UInt  , "imaginary-time axis or real-frequency axis"],
            "beta"     => [missing, 0, Real  , "inverse system temperature"],
            "niter"    => [missing, 0, UInt  , "number of iterations"],
            "mixer"    => [missing, 0, Real  , "mixing factor"],
            "dcount"   => [missing, 0, String, "scheme of double counting term"],
            "cc"       => [missing, 0, Real  , "convergence criterion of charge"],
            "ec"       => [missing, 0, Real  , "convergence criterion of total energy"],
            "fc"       => [missing, 0, Real  , "convergence criterion of force"],
            "lcharge"  => [missing, 0, Bool  , "examine whether charge is converged"],
            "lenergy"  => [missing, 0, Bool  , "examine whether total energy is converged"],
            "lforce"   => [missing, 0, Bool  , "examine whether force is converged"]
        )

"""
    PIMP

Dictionary for parameters: quantum impurity problems
"""
PIMP  = Dict{String,Any}(
            "nsite"    => [missing, 0, UInt  , "number of impurity sites"],
            "atoms"    => [missing, 0, Array , "chemical symbols of impurity atoms"],
            "equiv"    => [missing, 0, Array , "equivalency of quantum impurity atoms"],
            "shell"    => [missing, 0, Array , "angular momentum of correlated orbitals"],
            "ising"    => [missing, 0, Array , "interaction types of correlated orbitals"],
            "occup"    => [missing, 0, Array , "nominal impurity occupancy"],
            "upara"    => [missing, 0, Array , "Coulomb interaction parameter"],
            "jpara"    => [missing, 0, Array , "Hund's coupling parameter"],
            "lpara"    => [missing, 0, Array , "spin-orbit coupling parameter"]
        )

"""
    PSOLVER

Dictionary for parameters: quantum impurity solvers
"""
PSOLVER = Dict{String,Any}(
            "engine"   => [missing, 0, String, "name of quantum impurity solver"],
            "params"   => [missing, 0, Array , "parameters set of quantum impurity solver"]
        )

mutable struct IterInfo
    total_iter :: Integer
    dmft1_iter :: Integer
    dmft2_iter :: Integer
    dft_dmft_iter :: Integer
end

function IterInfo(iter::Integer = 0)
    IterInfo(iter, iter, iter, iter)
end

mutable struct Lattice
end

mutable struct DFTData
end
