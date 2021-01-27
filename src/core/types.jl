#
# Project : pansy
# Source  : types.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : unstable
# Comment :
#
# Last modified: 2021/01/28
#

#
# Customized Dictionaries
#

#
# Remarks:
#
# The values in the following dictionaries are actually arrays, which
# usually contain four elements:
#     Element[1] -> Actually value.
#     Element[2] -> If it is 1, this key-value pair is mandatory.
#                   If it is 0, this key-value pair is optional.
#     Element[3] -> Numerical type (A julia Symbol).
#     Element[4] -> Brief explanations.
#

"""
    PCASE

Dictionary for configuration parameters: case summary.
"""
const PCASE = Dict{String,Array{Any,1}}(
          "case"     => [missing, 1, :String, "System's name"]
      )

"""
    PDFT

Dictionary for configuration parameters: density functional theory calculations.
"""
const PDFT  = Dict{String,Array{Any,1}}(
          "engine"   => [missing, 1, :String, "Engine for density functional theory calculations"],
          "smear"    => [missing, 0, :String, "Scheme for smearing"],
          "kmesh"    => [missing, 0, :String, "K-mesh for brillouin zone sampling / integration"],
          "magmom"   => [missing, 0, :String, "Initial magnetic moments"],
          "lsymm"    => [missing, 0, :Bool  , "The symmetry is turned on or off"],
          "lspins"   => [missing, 0, :Bool  , "The spin orientations are polarized or not"],
          "lspinorb" => [missing, 0, :Bool  , "The spin-orbit coupling is considered or not"],
          "loptim"   => [missing, 0, :Bool  , "The generated projectors are optimized or not"],
          "lproj"    => [missing, 1, :Bool  , "The projectors are generated or not"],
          "sproj"    => [missing, 1, :Array , "Specifications for generating projectors"],
          "window"   => [missing, 0, :Array , "Band / energy window for orthogonalizing projectors"],
      )

"""
    PDMFT

Dictionary for configuration parameters: dynamical mean-field theory calculations.
"""
const PDMFT = Dict{String,Array{Any,1}}(
          "mode"     => [missing, 1, :I64   , "Scheme of dynamical mean-field theory calculations"],
          "axis"     => [missing, 1, :I64   , "Imaginary-time axis or real-frequency axis"],
          "niter"    => [missing, 1, :I64   , "Maximum number of iterations"],
          "dcount"   => [missing, 1, :String, "Scheme of double counting term"],
          "beta"     => [missing, 1, :F64   , "Inverse system temperature"],
          "mixer"    => [missing, 0, :F64   , "Mixing factor"],
          "cc"       => [missing, 0, :F64   , "Convergence criterion of charge"],
          "ec"       => [missing, 0, :F64   , "Convergence criterion of total energy"],
          "fc"       => [missing, 0, :F64   , "Convergence criterion of force"],
          "lcharge"  => [missing, 0, :Bool  , "Test whether charge is converged"],
          "lenergy"  => [missing, 0, :Bool  , "Test whether total energy is converged"],
          "lforce"   => [missing, 0, :Bool  , "Test whether force is converged"],
      )

"""
    PIMP

Dictionary for configuration parameters: quantum impurity problems.
"""
const PIMP  = Dict{String,Array{Any,1}}(
          "nsite"    => [missing, 1, :I64   , "Number of (correlated) impurity sites"],
          "atoms"    => [missing, 1, :Array , "Chemical symbols of impurity atoms"],
          "equiv"    => [missing, 1, :Array , "Equivalency of quantum impurity atoms"],
          "shell"    => [missing, 1, :Array , "Angular momenta of correlated orbitals"],
          "ising"    => [missing, 1, :Array , "Interaction types of correlated orbitals"],
          "occup"    => [missing, 1, :Array , "Nominal impurity occupancy"],
          "upara"    => [missing, 1, :Array , "Coulomb interaction parameter"],
          "jpara"    => [missing, 1, :Array , "Hund's coupling parameter"],
          "lpara"    => [missing, 1, :Array , "Spin-orbit coupling parameter"],
      )

"""
    PSOLVER

Dictionary for configuration parameters: quantum impurity solvers.
"""
const PSOLVER= Dict{String,Array{Any,1}}(
          "engine"   => [missing, 1, :String, "Name of quantum impurity solver"],
          "params"   => [missing, 1, :Array , "Parameter sets of quantum impurity solver"],
      )

#
# Customized Structs
#

"""
    Logger

Store the IOStreams for case.log and case.cycle files.

.log   -> IOStream for case.log file
.cycle -> IOStream for case.cycle file
"""
mutable struct Logger
    # Members
    log   :: IOStream
    cycle :: IOStream

    # Default constructor
    Logger() = new()
end

"""
    IterInfo

Record the DFT + DMFT iteration information.

.dmft1_iter -> Number of iterations between dmft1 and quantum impurity solver
.dmft2_iter -> Number of iterations between dmft2 and DFT engine
.dmft_cycle -> Number of DFT + DMFT iterations
.full_cycle -> Counter for each iteration
._dft_fermi -> Fermi level obtained by DFT engine
.dmft_fermi -> Fermi level obtained by DMFT engine (dmft1)
"""
mutable struct IterInfo
    dmft1_iter :: I64
    dmft2_iter :: I64
    dmft_cycle :: I64
    full_cycle :: I64
    _dft_fermi :: F64
    dmft_fermi :: F64
end

"""
    Lattice

Contain the crystallography information.

._case -> The name of system
.scale -> Universal scaling factor (lattice constant), which is used to
          scale all lattice vectors and all atomic coordinates.
.lvect -> Three lattice vectors defining the unit cell of the system. Its
          shape must be (3, 3).
.nsort -> Number of sorts of atoms
.natom -> Number of atoms
.sorts -> Sorts of atoms. Its shape must be (nsort, 2).
.atoms -> Lists of atoms. Its shape must be (natom).
.coord -> Atomic positions are provided in cartesian coordinates or in
          direct coordinates (respectively fractional coordinates). Its
          shape must be (natom, 3).
"""
mutable struct Lattice
    _case :: String
    scale :: F64
    lvect :: Array{F64,2}
    nsort :: I64
    natom :: I64
    sorts :: Array{Union{String,I64},2}
    atoms :: Array{String,1}
    coord :: Array{F64,2}
end

"""
    PrTrait

Essential information of a given projector.

.site -> Site in which the projector is defined
.l    -> Quantum number l
.m    -> Quantum number m
.desc -> Projector's specification
"""
mutable struct PrTrait
    site  :: I64
    l     :: I64
    m     :: I64
    desc  :: String
end

"""
    PrGroup

Essential information of group of projectors.

.site   -> Site in which the projectors are defined. In principle, the
           projectors included in the same group should be defined at
           the same site (or equivalently atom).
.l      -> Quantum number l. In principle, the projectors included in
           the same group should have the same quantum number l (but
           with different m).
.corr   -> Test if the projectors in this group are correlated
.shell  -> Type of correlated orbitals. It is infered from quantum number l.
.Pr     -> Array. It contains the indices of projectors.
.Tr     -> Array. It contains the transformation matrix. This parameter
           could be useful to select certain subset of orbitals or perform
           a simple global rotation.
.window -> Tuple. It is the band window or energy window, which is used
           to filter the Kohn-Sham band structure.
"""
mutable struct PrGroup
    site  :: I64
    l     :: I64
    corr  :: Bool
    shell :: String
    Pr    :: Array{I64,1}
    Tr    :: Array{C64,2}
    window:: Tuple{R64,R64}
end

"""
    PrUnion

Essential information of group of projectors (be transformed or rotated).

Note that the `site`, `l`, `corr`, and `shell` arguments are borrowed
from the PrGroup struct. In other words, they are just copies of those
in the PrGroup struct. On the other hand, the `ndim` is new, and the
`Pr` is different from the one in the PrGroup struct.

.site  -> Site in which the projectors are defined. In principle, the
          projectors included in the same group should be defined at
          the same site (or equivalently atom).
.l     -> Quantum number l. In principle, the projectors included in
          the same group should have the same quantum number l (but
          with different m).
.ndim  -> How many projectors are actually included in this group, which
          should be equal to the length of vector Pr.
.corr  -> Test if the projectors in this group are correlated
.shell -> Type of correlated orbitals. It is infered from the PIMP dict.
.Pr    -> Array. It contains the indices of projectors.
"""
mutable struct PrUnion
    site  :: I64
    l     :: I64
    ndim  :: I64
    corr  :: Bool
    shell :: String
    Pr    :: Array{I64,1}
end

#
# Customized Constructors
#

"""
    IterInfo(iter::I64 = 0, fermi::F64 = 0.0)

Outer constructor for IterInfo struct.
"""
function IterInfo(iter::I64 = 0, fermi::F64 = 0.0)
    IterInfo(iter, iter, iter, iter, fermi, fermi)
end

"""
    Lattice(_case::String, scale::F64, nsort::I64, natom::I64)

Outer constructor for Lattice struct.
"""
function Lattice(_case::String, scale::F64, nsort::I64, natom::I64)
    # Initialize the arrays
    lvect = zeros(F64, 3, 3)
    sorts = Array{Union{String,I64}}(undef, nsort, 2)
    atoms = fill("", natom)
    coord = zeros(F64, natom, 3)

    # Call the default constructor
    Lattice(_case, scale, lvect, nsort, natom, sorts, atoms, coord)
end

"""
    PrTrait(site::I64, sort::String, desc::String)

Outer constructor for PrTrait struct.
"""
function PrTrait(site::I64, desc::String)

#
# Remarks:
#
# Please go to the following webpage for more details about the original
# specifications of projectors in the vasp code:
#     https://www.vasp.at/wiki/index.php/LOCPROJ
#

    # Angular character of the local functions on the specified sites
    orb_labels = ("s",
                  "py", "pz", "px",
                  "dxy", "dyz", "dz2", "dxz", "dx2-y2",
                  "fz3", "fxz2", "fyz2", "fz(x2-y2)", "fxyz", "fx(x2-3y2)", "fy(3x2-y2)")

    # To make sure the specified desc is valid
    @assert desc in orb_labels

    # Determine quantum numbers l and m according to desc
    lm = findfirst(x -> x === desc, orb_labels) - 1
    l = convert(I64, floor(sqrt(lm)))
    m = lm - l * l

    # Call the default constructor
    PrTrait(site, l, m, desc)
end

"""
    PrGroup(site::I64, l::I64)

Outer constructor for PrGroup struct.
"""
function PrGroup(site::I64, l::I64)
    # The lshell defines a mapping from l (integer) to shell (string)
    lshell = Dict{I64,String}(
                 0 => "s",
                 1 => "p",
                 2 => "d",
                 3 => "f",
             )

    # Setup initial parameters
    # They will be further initialized in vaspio_projs() and plo_group()
    corr  = false
    shell = lshell[l]

    # Allocate memory for Pr and Tr
    # They will be further initialized in vaspio_projs() and plo_group()
    max_dim = 2 * l + 1
    Pr = zeros(I64, max_dim)
    Tr = zeros(C64, max_dim, max_dim)

    # Call the default constructor
    PrGroup(site, l, corr, shell, Pr, Tr)
end

"""
    PrUnion(site::I64, l::I64, ndim::I64, corr::Bool, shell::String)

Outer constructor for PrUnion struct.
"""
function PrUnion(site::I64, l::I64, ndim::I64, corr::Bool, shell::String)
    # Allocate memory for Pr
    Pr = zeros(I64, ndim)

    # Call the default constructor
    PrUnion(site, l, ndim, corr, shell, Pr)
end

"""
    PrUnion(PG::PrGroup)

Outer constructor for PrUnion struct.
"""
function PrUnion(PG::PrGroup)
    # Determine ndim
    ndim = size(PG.Tr)[1]

    # Allocate memory for Pr
    Pr = zeros(I64, ndim)

    # Call the default constructor
    PrUnion(PG.site, PG.l, ndim, PG.corr, PG.shell, Pr)
end

#
# Customized Base.show() Functions
#

"""
    Base.show(io::IO, it::IterInfo)

Base.show() function for Logger struct.
"""
function Base.show(io::IO, logger::Logger)
    println(io, "Logger struct")
    println(io, ".log   : ", logger.log  )
    println(io, ".cycle : ", logger.cycle)
end

"""
    Base.show(io::IO, it::IterInfo)

Base.show() function for IterInfo struct.
"""
function Base.show(io::IO, it::IterInfo)
    println(io, "IterInfo struct")
    println(io, ".dmft1_iter : ", it.dmft1_iter)
    println(io, ".dmft2_iter : ", it.dmft2_iter)
    println(io, ".dmft_cycle : ", it.dmft_cycle)
    println(io, ".full_cycle : ", it.full_cycle)
    println(io, "._dft_fermi : ", it._dft_fermi)
    println(io, ".dmft_fermi : ", it.dmft_fermi)
end

"""
    Base.show(io::IO, latt::Lattice)

Base.show() function for Lattice struct.
"""
function Base.show(io::IO, latt::Lattice)
    println(io, "Lattice struct")
    println(io, "._case : ", latt._case)
    println(io, ".scale : ", latt.scale)
    println(io, ".lvect : ", latt.lvect)
    println(io, ".nsort : ", latt.nsort)
    println(io, ".natom : ", latt.natom)
    println(io, ".sorts : ", latt.sorts)
    println(io, ".atoms : ", latt.atoms)
    println(io, ".coord : ", latt.coord)
end

"""
    Base.show(io::IO, PT::PrTrait)

Base.show() function for PrTrait struct.
"""
function Base.show(io::IO, PT::PrTrait)
    println(io, "PrTrait struct")
    println(io, ".site : ", PT.site)
    println(io, ".l    : ", PT.l)
    println(io, ".m    : ", PT.m)
    println(io, ".desc : ", PT.desc)
end

"""
    Base.show(io::IO, PG::PrGroup)

Base.show() function for PrGroup struct.
"""
function Base.show(io::IO, PG::PrGroup)
    println(io, "PrGroup struct")
    println(io, ".site  : ", PG.site)
    println(io, ".l     : ", PG.l)
    println(io, ".corr  : ", PG.corr)
    println(io, ".shell : ", PG.shell)
    println(io, ".Pr    : ", PG.Pr)
    println(io, ".Tr    : ", PG.Tr)
end

"""
    Base.show(io::IO, PU::PrUnion)

Base.show() function for PrUnion struct.
"""
function Base.show(io::IO, PU::PrUnion)
    println(io, "PrUnion struct")
    println(io, ".site  : ", PU.site)
    println(io, ".l     : ", PU.l)
    println(io, ".ndim  : ", PU.ndim)
    println(io, ".corr  : ", PU.corr)
    println(io, ".shell : ", PU.shell)
    println(io, ".Pr    : ", PU.Pr)
end
