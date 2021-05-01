#
# Project : Pansy
# Source  : types.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/05/01
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

See also: [`PDFT`](@ref), [`PDMFT`](@ref), [`PIMP`](@ref), [`PSOLVER`](@ref).
"""
const PCASE = Dict{String,Array{Any,1}}(
          "case"     => [missing, 1, :String, "System's name"]
      )

"""
    PDFT

Dictionary for configuration parameters: density functional theory calculations.

See also: [`PCASE`](@ref), [`PDMFT`](@ref), [`PIMP`](@ref), [`PSOLVER`](@ref).
"""
const PDFT  = Dict{String,Array{Any,1}}(
          "engine"   => [missing, 1, :String, "Engine for density functional theory calculations"],
          "projtype" => [missing, 1, :String, "Types of projectors"],
          "smear"    => [missing, 1, :String, "Scheme for smearing"],
          "kmesh"    => [missing, 1, :String, "K-mesh for brillouin zone sampling / integration"],
          "magmom"   => [missing, 0, :String, "Initial magnetic moments"],
          "lsymm"    => [missing, 1, :Bool  , "The symmetry is turned on or off"],
          "lspins"   => [missing, 1, :Bool  , "The spin orientations are polarized or not"],
          "lspinorb" => [missing, 1, :Bool  , "The spin-orbit coupling is considered or not"],
          "loptim"   => [missing, 1, :Bool  , "The generated projectors are optimized or not"],
          "lproj"    => [missing, 1, :Bool  , "The projectors are generated or not"],
          "sproj"    => [missing, 1, :Array , "Specifications for generating projectors"],
          "window"   => [missing, 1, :Array , "Band / energy window for orthogonalizing projectors"],
      )

"""
    PDMFT

Dictionary for configuration parameters: dynamical mean-field theory calculations.

See also: [`PCASE`](@ref), [`PDFT`](@ref), [`PIMP`](@ref), [`PSOLVER`](@ref).
"""
const PDMFT = Dict{String,Array{Any,1}}(
          "mode"     => [missing, 1, :I64   , "Scheme of dynamical mean-field theory calculations"],
          "axis"     => [missing, 1, :I64   , "Imaginary-time axis or real-frequency axis"],
          "niter"    => [missing, 1, :I64   , "Maximum number of all iterations"],
          "nmesh"    => [missing, 1, :I64   , "Number of frequency points"],
          "dcount"   => [missing, 1, :String, "Scheme of double counting term"],
          "beta"     => [missing, 1, :F64   , "Inverse system temperature"],
          "mixer"    => [missing, 1, :F64   , "Mixing factor"],
          "mc"       => [missing, 0, :F64   , "Convergence criterion of chemical potential"],
          "cc"       => [missing, 0, :F64   , "Convergence criterion of charge"],
          "ec"       => [missing, 0, :F64   , "Convergence criterion of total energy"],
          "fc"       => [missing, 0, :F64   , "Convergence criterion of force"],
          "lfermi"   => [missing, 0, :Bool  , "Test whether chemical potential is updated"],
          "lcharge"  => [missing, 0, :Bool  , "Test whether charge is converged"],
          "lenergy"  => [missing, 0, :Bool  , "Test whether total energy is converged"],
          "lforce"   => [missing, 0, :Bool  , "Test whether force is converged"],
      )

"""
    PIMP

Dictionary for configuration parameters: quantum impurity problems.

See also: [`PCASE`](@ref), [`PDFT`](@ref), [`PDMFT`](@ref), [`PSOLVER`](@ref).
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

See also: [`PCASE`](@ref), [`PDFT`](@ref), [`PDMFT`](@ref), [`PIMP`](@ref).
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

Mutable struct. Store the IOStreams for case.log and case.cycle files.

## Members

log   -> IOStream for case.log file.\n
cycle -> IOStream for case.cycle file.

See also: [`IterInfo`](@ref).
"""
mutable struct Logger
    log   :: IOStream
    cycle :: IOStream
end

"""
    IterInfo

Mutable struct. Record the DFT + DMFT iteration information.

## Members

dmft1_iter -> Number of iterations between dmft1 and quantum impurity solver.\n
dmft2_iter -> Number of iterations between dmft2 and DFT engine.\n
dmft_cycle -> Number of DFT + DMFT iterations.\n
full_cycle -> Counter for each iteration.\n
dft_fermi  -> Fermi level obtained by DFT engine.\n
dmft_fermi -> Fermi level obtained by DMFT engine (dmft1).

See also: [`Logger`](@ref).
"""
mutable struct IterInfo
    dmft1_iter :: I64
    dmft2_iter :: I64
    dmft_cycle :: I64
    full_cycle :: I64
     dft_fermi :: F64
    dmft_fermi :: F64
end

"""
    Lattice

Mutable struct. Contain the crystallography information. This struct is
designed for the `POSCAR` file used by the vasp code.

## Members

_case -> The name of system.\n
scale -> Universal scaling factor (lattice constant), which is used to
         scale all lattice vectors and all atomic coordinates.\n
lvect -> Three lattice vectors defining the unit cell of the system. Its
         shape must be (3, 3).\n
nsort -> Number of sorts of atoms.\n
natom -> Number of atoms.\n
sorts -> Sorts of atoms. Its shape must be (nsort, 2).\n
atoms -> Lists of atoms. Its shape must be (natom).\n
coord -> Atomic positions are provided in cartesian coordinates or in
         direct coordinates (respectively fractional coordinates). Its
         shape must be (natom, 3).

See also: [`vaspio_lattice`](@ref).
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
    Mapping

Mutable struct. Mapping between quantum impurity problems and groups
of projectors.

## Members

i_grp -> Mapping from impurity problems to groups of projectors.\n
i_wnd -> Mapping from impurity problems to windows of dft bands.\n
g_imp -> Mapping from groups of projectors to impurity problems.\n
w_imp -> Mapping from windows of dft bands to impurity problems.

See also: [`PrGroup`](@ref).
"""
mutable struct Mapping
    i_grp :: Array{I64,1}
    i_wnd :: Array{I64,1}
    g_imp :: Array{I64,1}
    w_imp :: Array{I64,1}
end

"""
    PrTrait

Mutable struct. Essential information of a given projector.

## Members

site -> Site in which the projector is defined.\n
l    -> Quantum number l.\n
m    -> Quantum number m.\n
desc -> Projector's specification.

See also: [`PrGroup`](@ref), [`PrWindow`](@ref).
"""
mutable struct PrTrait
    site  :: I64
    l     :: I64
    m     :: I64
    desc  :: String
end

"""
    PrGroup

Mutable struct. Essential information of group of projectors.

## Members

site   -> Site in which the projectors are defined. In principle, the
          projectors included in the same group should be defined at
          the same site (or equivalently atom).\n
l      -> Quantum number l. In principle, the projectors included in
          the same group should have the same quantum number l (but
          with different m).\n
corr   -> Test if the projectors in this group are correlated.\n
shell  -> Type of correlated orbitals. It is infered from quantum number l.\n
Pr     -> Array. It contains the indices of projectors.\n
Tr     -> Array. It contains the transformation matrix. This parameter
          could be useful to select certain subset of orbitals or perform
          a simple global rotation.

See also: [`PrTrait`](@ref), [`PrWindow`](@ref), [`Mapping`](@ref).
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
    PrWindow

Mutable struct. Define the band window for group of projectors.

## Members

bmin -> Minimum band index.\n
bmax -> Maximum band index.\n
nbnd -> Maximum number of bands in the current window (= bmax-bmin+1).\n
kwin -> Momentum-dependent and spin-dependent band window.\n
bwin -> Tuple. It is the band window or energy window, which is used
        to filter the Kohn-Sham band structure. The mesh for calculating
        density of states is also deduced from `bwin`.

See also: [`PrTrait`](@ref), [`PrGroup`](@ref).
"""
mutable struct PrWindow
    bmin  :: I64
    bmax  :: I64
    nbnd  :: I64
    kwin  :: Array{I64,3}
    bwin  :: Tuple{R64,R64}
end

#
# Customized Constructors
#

"""
    Logger(case::String = "case")

Outer constructor for Logger struct.
"""
function Logger(case::String = "case")
    # Open the streams for IO
    log = open("$case.log", "a")
    cycle = open("$case.cycle", "a")

    # Call the default constructor
    Logger(log, cycle)
end

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
    Mapping(nsite::I64, ngrp::I64)

Outer constructor for Mapping struct.
"""
function Mapping(nsite::I64, ngrp::I64)
    # Sanity check
    @assert ngrp >= nsite

    # Initialize the arrays
    i_grp = zeros(I64, nsite)
    g_imp = zeros(I64, ngrp)

    # Call the default constructor
    Mapping(i_grp, g_imp)
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
    # The dict lshell defines a mapping from l (integer) to shell (string)
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
    PrWindow(kwin::Array{I64,3}, bwin::Tuple{R64,R64})

Outer constructor for PrWindow struct.
"""
function PrWindow(kwin::Array{I64,3}, bwin::Tuple{R64,R64})
    # Determine the boundaries of the window
    bmin = minimum(kwin[:, :, 1])
    bmax = maximum(kwin[:, :, 2])

    # Determine the size of the window
    nbnd = bmax - bmin + 1

    # Call the default constructor
    PrWindow(bmin, bmax, nbnd, kwin, bwin)
end

#
# Customized Base.show() Functions
#

"""
    Base.show(io::IO, it::IterInfo)

Base.show() function for Logger struct.

See also: [`Logger`](@ref).
"""
function Base.show(io::IO, logger::Logger)
    # To make sure these IOStreams are alive
    @assert isopen(logger.log) && isopen(logger.cycle)

    println(io, "Logger struct")
    println(io, "log   : ", logger.log  )
    println(io, "cycle : ", logger.cycle)
end

"""
    Base.show(io::IO, it::IterInfo)

Base.show() function for IterInfo struct.

See also: [`IterInfo`](@ref).
"""
function Base.show(io::IO, it::IterInfo)
    println(io, "IterInfo struct")
    println(io, "dmft1_iter : ", it.dmft1_iter)
    println(io, "dmft2_iter : ", it.dmft2_iter)
    println(io, "dmft_cycle : ", it.dmft_cycle)
    println(io, "full_cycle : ", it.full_cycle)
    println(io, "dft_fermi  : ", it.dft_fermi)
    println(io, "dmft_fermi : ", it.dmft_fermi)
end

"""
    Base.show(io::IO, latt::Lattice)

Base.show() function for Lattice struct.

See also: [`Lattice`](@ref).
"""
function Base.show(io::IO, latt::Lattice)
    println(io, "Lattice struct")
    println(io, "_case : ", latt._case)
    println(io, "scale : ", latt.scale)
    println(io, "lvect : ", latt.lvect)
    println(io, "nsort : ", latt.nsort)
    println(io, "natom : ", latt.natom)
    println(io, "sorts : ", latt.sorts)
    println(io, "atoms : ", latt.atoms)
    println(io, "coord : ", latt.coord)
end

"""
    Base.show(io::IO, map::Mapping)

Base.show() function for Mapping struct.

See also: [`Mapping`](@ref).
"""
function Base.show(io::IO, map::Mapping)
    println(io, "Mapping struct")
    println(io, "i_grp : ", map.i_grp)
    println(io, "g_imp : ", map.g_imp)
end

"""
    Base.show(io::IO, PT::PrTrait)

Base.show() function for PrTrait struct.

See also: [`PrTrait`](@ref).
"""
function Base.show(io::IO, PT::PrTrait)
    println(io, "PrTrait struct")
    println(io, "site : ", PT.site)
    println(io, "l    : ", PT.l)
    println(io, "m    : ", PT.m)
    println(io, "desc : ", PT.desc)
end

"""
    Base.show(io::IO, PG::PrGroup)

Base.show() function for PrGroup struct.

See also: [`PrGroup`](@ref).
"""
function Base.show(io::IO, PG::PrGroup)
    println(io, "PrGroup struct")
    println(io, "site  : ", PG.site)
    println(io, "l     : ", PG.l)
    println(io, "corr  : ", PG.corr)
    println(io, "shell : ", PG.shell)
    println(io, "Pr    : ", PG.Pr)
    println(io, "Tr    : ", PG.Tr)
end

"""
    Base.show(io::IO, PW::PrWindow)

Base.show() function for PrWindow struct.

See also: [`PrWindow`](@ref).
"""
function Base.show(io::IO, PW::PrWindow)
    println(io, "PrWindow struct")
    println(io, "bmin : ", PW.bmin)
    println(io, "bmax : ", PW.bmax)
    println(io, "nbnd : ", PW.nbnd)
    println(io, "kwin : ", PW.kwin)
    println(io, "bwin : ", PW.bwin)
end
