#
# project : pansy
# source  : types.jl
# author  : Li Huang (lihuang.dmft@gmail.com)
# status  : unstable
# comment :
#
# last modified: 2021/01/12
#

#
# Customized Dictionaries
#

#
# remarks:
#
# the values in the following dictionaries are actually arrays, which
# contain four elements
#     value[1] -> actually value
#     value[2] -> if it is 1, this key-value pair is mandatory
#                 if it is 0, this key-value pair is optional
#     value[3] -> numerical type (a julia Symbol)
#     value[4] -> brief explanations
#

"""
    PCASE

Dictionary for configuration parameters: case summary
"""
const PCASE = Dict{String,Array{Any,1}}(
          "case"     => [missing, 1, :String, "system's name"]
      )

"""
    PDFT

Dictionary for configuration parameters: density functional theory calculations
"""
const PDFT  = Dict{String,Array{Any,1}}(
          "engine"   => [missing, 1, :String, "engine for density functional theory calculations"],
          "smear"    => [missing, 0, :String, "scheme for smearing"],
          "kmesh"    => [missing, 0, :String, "kmesh for brillouin zone sampling / integration"],
          "magmom"   => [missing, 0, :String, "initial magnetic moments"],
          "lsymm"    => [missing, 0, :Bool  , "the symmetry is turned on or off"],
          "lspins"   => [missing, 0, :Bool  , "the spin orientations are polarized or not"],
          "lspinorb" => [missing, 0, :Bool  , "the spin-orbit coupling is considered or not"],
          "loptim"   => [missing, 0, :Bool  , "the generated projectors are optimized or not"],
          "lproj"    => [missing, 1, :Bool  , "the projectors are generated or not"],
          "ewidth"   => [missing, 0, :F64   , "half-width of energy window for generating optimal projectors"],
          "sproj"    => [missing, 1, :Array , "strings / descriptions for generating projectors"],
      )

"""
    PDMFT

Dictionary for configuration parameters: dynamical mean-field theory calculations
"""
const PDMFT = Dict{String,Array{Any,1}}(
          "mode"     => [missing, 1, :I64   , "scheme of dynamical mean-field theory calculations"],
          "axis"     => [missing, 1, :I64   , "imaginary-time axis or real-frequency axis"],
          "niter"    => [missing, 0, :I64   , "maximum number of iterations"],
          "dcount"   => [missing, 0, :String, "scheme of double counting term"],
          "beta"     => [missing, 1, :F64   , "inverse system temperature"],
          "mixer"    => [missing, 0, :F64   , "mixing factor"],
          "cc"       => [missing, 0, :F64   , "convergence criterion of charge"],
          "ec"       => [missing, 0, :F64   , "convergence criterion of total energy"],
          "fc"       => [missing, 0, :F64   , "convergence criterion of force"],
          "lcharge"  => [missing, 0, :Bool  , "test whether charge is converged"],
          "lenergy"  => [missing, 0, :Bool  , "test whether total energy is converged"],
          "lforce"   => [missing, 0, :Bool  , "test whether force is converged"],
      )

"""
    PIMP

Dictionary for configuration parameters: quantum impurity problems
"""
const PIMP  = Dict{String,Array{Any,1}}(
          "nsite"    => [missing, 1, :I64   , "number of (correlated) impurity sites"],
          "atoms"    => [missing, 1, :Array , "chemical symbols of impurity atoms"],
          "equiv"    => [missing, 1, :Array , "equivalency of quantum impurity atoms"],
          "shell"    => [missing, 1, :Array , "angular momenta of correlated orbitals"],
          "ising"    => [missing, 1, :Array , "interaction types of correlated orbitals"],
          "occup"    => [missing, 1, :Array , "nominal impurity occupancy"],
          "upara"    => [missing, 1, :Array , "Coulomb interaction parameter"],
          "jpara"    => [missing, 1, :Array , "Hund's coupling parameter"],
          "lpara"    => [missing, 1, :Array , "spin-orbit coupling parameter"],
      )

"""
    PSOLVER

Dictionary for configuration parameters: quantum impurity solvers
"""
const PSOLVER= Dict{String,Array{Any,1}}(
          "engine"   => [missing, 1, :String, "name of quantum impurity solver"],
          "params"   => [missing, 1, :Array , "parameter sets of quantum impurity solver"],
      )

#
# Customized Structs
#

"""
    IterInfo

Record the runtime information

.dmft1_iter -> number of iterations between dmft1 and quantum impurity solver
.dmft2_iter -> number of iterations between dmft2 and dft engine
.dmft_cycle -> number of dft + dmft iterations
.full_cycle -> counter for each iteration
._dft_fermi -> fermi level obtained by dft engine
.dmft_fermi -> fermi level obtained by dmft engine (dmft1)
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

Contain the crystallography information

._case -> the name of system
.scale -> universal scaling factor (lattice constant), which is used to
          scale all lattice vectors and all atomic coordinates
.lvect -> three lattice vectors defining the unit cell of the system. its
          shape must be (3, 3)
.nsort -> number of sorts of atoms
.natom -> number of atoms
.sorts -> sorts of atoms. its shape must be (nsort, 2)
.atoms -> lists of atoms. its shape must be (natom)
.coord -> atomic positions are provided in cartesian coordinates or in
          direct coordinates (respectively fractional coordinates). its
          shape must be (natom, 3)
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

Essential information of a given projector

.site -> site in which the projector is defined
.l    -> quantum number l
.m    -> quantum number m
.desc -> projector's specification
"""
mutable struct PrTrait
    site  :: I64
    l     :: I64
    m     :: I64
    desc  :: String
end

"""
    PrGroup

Essential information of group of projectors

.site  -> site in which the projectors are defined. in principle, the
          projectors included in the same group should be defined at
          the same site (or equivalently atom)
.l     -> quantum number l. in principle, the projectors included in
          the same group should have the same quantum number l (but
          with different m)
.corr  -> if the projectors in this group are correlated
.shell -> type of correlated orbitals. it is infered from quantum number l
.Pr    -> array. it contains the indices of projectors
.Tr    -> array. it contains the transformation matrix. this parameter
          can be useful to select certain subset of orbitals or perform
          a simple global rotation
"""
mutable struct PrGroup
    site  :: I64
    l     :: I64
    corr  :: Bool
    shell :: String
    Pr    :: Array{I64,1}
    Tr    :: Array{C64,2}
end

"""
    PrGroupT

Essential information of group of projectors (be transformed or rotated)

.site  -> site in which the projectors are defined. in principle, the
          projectors included in the same group should be defined at
          the same site (or equivalently atom)
.l     -> quantum number l. in principle, the projectors included in
          the same group should have the same quantum number l (but
          with different m)
.ndim  -> how many projectors are actually included in this group, which
          should be equal to the length of vector Pr
.corr  -> if the projectors in this group are correlated
.shell -> type of correlated orbitals. it is infered from the PIMP dict
.Pr    -> array. it contains the indices of projectors
"""
mutable struct PrGroupT
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

Outer constructor for IterInfo struct
"""
function IterInfo(iter::I64 = 0, fermi::F64 = 0.0)
    IterInfo(iter, iter, iter, iter, fermi, fermi)
end

"""
    Lattice(_case::String, scale::F64, nsort::I64, natom::I64)

Outer constructor for Lattice struct
"""
function Lattice(_case::String, scale::F64, nsort::I64, natom::I64)
    # initialize the arrays
    lvect = zeros(F64, 3, 3)
    sorts = Array{Union{String,I64}}(undef, nsort, 2)
    atoms = fill("", natom)
    coord = zeros(F64, natom, 3)

    # call the default constructor
    Lattice(_case, scale, lvect, nsort, natom, sorts, atoms, coord)
end

"""
    PrTrait(site::I64, sort::String, desc::String)

Outer constructor for PrTrait struct
"""
function PrTrait(site::I64, desc::String)
    # angular character of the local functions on the specified sites
    # see the following webpage for more details
    #     https://www.vasp.at/wiki/index.php/LOCPROJ
    orb_labels = ("s",
                  "py", "pz", "px",
                  "dxy", "dyz", "dz2", "dxz", "dx2-y2",
                  "fz3", "fxz2", "fyz2", "fz(x2-y2)", "fxyz", "fx(x2-3y2)", "fy(3x2-y2)")

    # to make sure the specified desc is valid
    @assert desc in orb_labels

    # determine quantum numbers l and m according to desc
    lm = findfirst(x -> x === desc, orb_labels) - 1
    l = convert(I64, floor(sqrt(lm)))
    m = lm - l * l

    # call the default constructor
    PrTrait(site, l, m, desc)
end

"""
    PrGroup(site::I64, l::I64)

Outer constructor for PrGroup struct
"""
function PrGroup(site::I64, l::I64)
    # lshell defines a mapping from l (integer) to shell (string)
    lshell = Dict{I64,String}(
                 0 => "s",
                 1 => "p",
                 2 => "d",
                 3 => "f",
             )

    # setup initial parameters
    # they will be further initialized in vaspio_projs() and plo_group()
    corr  = false
    shell = lshell[l]

    # allocate memory for Pr and Tr
    # they will be further initialized in vaspio_projs() and plo_group()
    max_dim = 2 * l + 1
    Pr = zeros(I64, max_dim)
    Tr = zeros(C64, max_dim, max_dim)

    # call the default constructor
    PrGroup(site, l, corr, shell, Pr, Tr)
end

"""
    PrGroupT(site::I64, l::I64, ndim::I64, corr::Bool, shell::String)

Outer constructor for PrGroupT struct
"""
function PrGroupT(site::I64, l::I64, ndim::I64, corr::Bool, shell::String)
    # allocate memory for Pr
    Pr = zeros(I64, ndim)

    # call the default constructor
    PrGroupT(site, l, ndim, corr, shell, Pr)
end

"""
    PrGroupT(PG::PrGroup)

Outer constructor for PrGroupT struct
"""
function PrGroupT(PG::PrGroup)
    # determine ndim
    ndim = size(PG.Tr)[1]

    # allocate memory for Pr
    Pr = zeros(I64, ndim)

    # call the default constructor
    PrGroupT(PG.site, PG.l, ndim, PG.corr, PG.shell, Pr)
end

#
# Customized Base.show() functions
#

"""
    Base.show(io::IO, it::IterInfo)

Base.show() function for IterInfo struct
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

Base.show() function for Lattice struct
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

Base.show() function for PrTrait struct
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

Base.show() function for PrGroup struct
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
    Base.show(io::IO, PGT::PrGroupT)

Base.show() function for PrGroupT struct
"""
function Base.show(io::IO, PGT::PrGroupT)
    println(io, "PrGroupT struct")
    println(io, ".site  : ", PGT.site)
    println(io, ".l     : ", PGT.l)
    println(io, ".ndim  : ", PGT.ndim)
    println(io, ".corr  : ", PGT.corr)
    println(io, ".shell : ", PGT.shell)
    println(io, ".Pr    : ", PGT.Pr)
end
