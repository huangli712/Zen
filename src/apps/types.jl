"""
    PCASE

Dictionary for configuration parameters: case summary    
"""
PCASE = Dict{String,Any}(
            "case"     => [missing, 1, String, "system's name"]
        )

"""
    PDFT

Dictionary for configuration parameters: density functional theory calculations
"""
PDFT  = Dict{String,Any}(
            "engine"   => [missing, 1, String, "engine for density functional theory calculations"],
            "smear"    => [missing, 1, String, "scheme for smearing"],
            "kmesh"    => [missing, 1, String, "density of k-mesh sampling in the brillouin zone"],
            "magmom"   => [missing, 0, String, "initial magnetic moment"],
            "lsymm"    => [missing, 1, Bool  , "whether the symmetry is considered"],
            "lspins"   => [missing, 1, Bool  , "whether the spin orientations are polarized"],
            "lspinorb" => [missing, 1, Bool  , "whether the spin-orbit coupling is considered"],
            "window"   => [missing, 0, Array , "energy window for generating optimal projectors"],
            "loptim"   => [missing, 0, Bool  , "try to optimize the generated projectors"],
            "lproj"    => [missing, 1, Bool  , "try to generate projectors"],
            "nproj"    => [missing, 1, I64   , "number of types of projectors"],
            "sproj"    => [missing, 1, Array , "scheme for generating projectors"]
        )

"""
    PDMFT

Dictionary for configuration parameters: dynamical mean-field theory calculations
"""
PDMFT = Dict{String,Any}(
            "mode"     => [missing, 1, I64   , "scheme of dynamical mean-field theory calculations"],
            "axis"     => [missing, 1, I64   , "imaginary-time axis or real-frequency axis"],
            "beta"     => [missing, 1, Real  , "inverse system temperature"],
            "niter"    => [missing, 1, I64   , "number of iterations"],
            "mixer"    => [missing, 1, Real  , "mixing factor"],
            "dcount"   => [missing, 1, String, "scheme of double counting term"],
            "cc"       => [missing, 1, Real  , "convergence criterion of charge"],
            "ec"       => [missing, 1, Real  , "convergence criterion of total energy"],
            "fc"       => [missing, 0, Real  , "convergence criterion of force"],
            "lcharge"  => [missing, 1, Bool  , "examine whether charge is converged"],
            "lenergy"  => [missing, 1, Bool  , "examine whether total energy is converged"],
            "lforce"   => [missing, 0, Bool  , "examine whether force is converged"]
        )

"""
    PIMP

Dictionary for configuration parameters: quantum impurity problems
"""
PIMP  = Dict{String,Any}(
            "nsite"    => [missing, 1, I64   , "number of impurity sites"],
            "atoms"    => [missing, 1, Array , "chemical symbols of impurity atoms"],
            "equiv"    => [missing, 1, Array , "equivalency of quantum impurity atoms"],
            "shell"    => [missing, 1, Array , "angular momentum of correlated orbitals"],
            "ising"    => [missing, 1, Array , "interaction types of correlated orbitals"],
            "occup"    => [missing, 1, Array , "nominal impurity occupancy"],
            "upara"    => [missing, 1, Array , "Coulomb interaction parameter"],
            "jpara"    => [missing, 1, Array , "Hund's coupling parameter"],
            "lpara"    => [missing, 1, Array , "spin-orbit coupling parameter"]
        )

"""
    PSOLVER

Dictionary for configuration parameters: quantum impurity solvers
"""
PSOLVER = Dict{String,Any}(
            "engine"   => [missing, 1, String, "name of quantum impurity solver"],
            "params"   => [missing, 1, Array , "parameters set of quantum impurity solver"]
        )

"""
    IterInfo

Record the runtime information
"""
mutable struct IterInfo
    total_iter :: I64
    dmft1_iter :: I64
    dmft2_iter :: I64
    dft_dmft_iter :: I64
end

"""
    IterInfo(iter::I64 = 0)

Outer constructor for IterInfo struct
"""
function IterInfo(iter::I64 = 0)
    IterInfo(iter, iter, iter, iter)
end
