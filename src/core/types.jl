#
# Project : Pansy
# Source  : types.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/06/24
#

#=
### *Customized Dictionaries*
=#

#=
*Remarks*:

The values in the following dictionaries are actually arrays, which
usually contain four elements:
* Element[1] -> Actually value.
* Element[2] -> If it is 1, this key-value pair is mandatory.
                If it is 0, this key-value pair is optional.
* Element[3] -> Numerical type (A julia Symbol).
* Element[4] -> Brief explanations.

The following dictionaries are used as global variables.
=#

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
          "window"   => [missing, 1, :Array , "Band / energy window for normalizing projectors"],
      )

"""
    PDMFT

Dictionary for configuration parameters: dynamical mean-field theory calculations.

See also: [`PCASE`](@ref), [`PDFT`](@ref), [`PIMP`](@ref), [`PSOLVER`](@ref).
"""
const PDMFT = Dict{String,Array{Any,1}}(
          "mode"     => [missing, 1, :I64   , "Scheme of dynamical mean-field theory calculations"],
          "axis"     => [missing, 1, :I64   , "Imaginary-time axis or real-frequency axis"],
          "niter"    => [missing, 1, :Array , "Maximum number of all iterations"],
          "nmesh"    => [missing, 1, :I64   , "Number of frequency points"],
          "dcount"   => [missing, 1, :String, "Scheme of double counting term"],
          "beta"     => [missing, 1, :F64   , "Inverse system temperature"],
          "mixer"    => [missing, 1, :F64   , "Mixing factor"],
          "mc"       => [missing, 0, :F64   , "Convergence criterion of chemical potential"],
          "cc"       => [missing, 0, :F64   , "Convergence criterion of charge"],
          "ec"       => [missing, 0, :F64   , "Convergence criterion of total energy"],
          "sc"       => [missing, 0, :F64   , "Convergence criterion of self-energy function"],
          "lfermi"   => [missing, 0, :Bool  , "Test whether chemical potential is updated"],
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
          "params"   => [missing, 1, :Array , "Extra parameter sets of quantum impurity solver"],
      )

#=
### *Customized Structs*
=#

"""
    Logger

Mutable struct. Store the IOStreams for case.log and case.cycle files.

### Members

* log   -> IOStream for case.log file.
* cycle -> IOStream for case.cycle file.

See also: [`IterInfo`](@ref).
"""
mutable struct Logger
    log   :: IOStream
    cycle :: IOStream
end

"""
    IterInfo

Mutable struct. Record the DFT + DMFT iteration information.

### Members

* I₁ -> Number of iterations between `dmft1` and quantum impurity solver.
* I₂ -> Number of iterations between `dmft2` and DFT engine.
* I₃ -> Number of DFT + DMFT iterations.
* I₄ -> Counter for each iteration.
* M₁ -> Maximum allowed number of iterations (between `dmft1` and solver).
* M₂ -> Maximum allowed number of iterations (between `dmft2` and DFT).
* M₃ -> Maximum allowed number of DFT + DMFT iterations.
* sc -> Self-consistent mode.
* μ₀ -> Fermi level obtained by DFT engine.
* μ₁ -> Fermi level obtained by DMFT engine (`dmft1`).
* μ₂ -> Fermi level obtained by DMFT engine (`dmft2`).
* dc -> Double counting terms.
* n₁ -> Number of lattice occupancy obtained by DMFT engine (`dmft1`).
* n₂ -> Number of lattice occupancy obtained by DMFT engine (`dmft2`).
* nf -> Number of impurity occupancy obtained by impurity solver.
* et -> Total DFT + DMFT energy.
* cc -> Convergence flag for charge density.
* ce -> Convergence flag for total energy.
* cs -> Convergence flag for self-energy functions.

See also: [`Logger`](@ref).
"""
mutable struct IterInfo
    I₁ :: I64
    I₂ :: I64
    I₃ :: I64
    I₄ :: I64
    M₁ :: I64
    M₂ :: I64
    M₃ :: I64
    sc :: I64
    μ₀ :: F64
    μ₁ :: F64
    μ₂ :: F64
    dc :: Vector{F64}
    n₁ :: F64
    n₂ :: F64
    nf :: Vector{F64}
    et :: F64
    cc :: Bool
    ce :: Bool
    cs :: Bool
end

"""
    Lattice

Mutable struct. Contain the crystallography information. This struct is
designed for the `POSCAR` file used by the `vasp` code.

### Members

* _case -> The name of system.
* scale -> Universal scaling factor (lattice constant), which is used to
           scale all lattice vectors and all atomic coordinates.
* lvect -> Three lattice vectors defining the unit cell of the system. Its
           shape must be (3, 3).
* nsort -> Number of sorts of atoms.
* natom -> Number of atoms.
* sorts -> Sorts of atoms. Its shape must be (nsort, 2).
* atoms -> Lists of atoms. Its shape must be (natom).
* coord -> Atomic positions are provided in cartesian coordinates or in
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
of projectors (or band windows).

### Members

* i_grp -> Mapping from quntum impurity problems to groups of projectors.
* i_wnd -> Mapping from quantum impurity problems to windows of dft bands.
* g_imp -> Mapping from groups of projectors to quantum impurity problems.
* w_imp -> Mapping from windows of dft bands to quantum impurity problems.

See also: [`Impurity`](@ref), [`PrGroup`](@ref), [`PrWindow`](@ref).
"""
mutable struct Mapping
    i_grp :: Array{I64,1}
    i_wnd :: Array{I64,1}
    g_imp :: Array{I64,1}
    w_imp :: Array{I64,1}
end

"""
    Impurity

Mutable struct. Essential information of quantum impurity problem.

### Members

* index -> Index of the quantum impurity problem.
* atoms -> Chemical symbol of impurity atom.
* sites -> Index of impurity atom.
* shell -> Angular momentum of correlated orbitals.
* ising -> Interaction type of correlated orbitals.
* occup -> Impurity occupancy.
* upara -> Coulomb interaction parameter.
* jpara -> Hund's coupling parameter.
* lpara -> Spin-orbit coupling parameter.
* beta  -> Inverse temperature.
* nband -> Number of correlated orbitals (spin is not included).

See also: [`Mapping`](@ref), [`PrGroup`](@ref), [`PrWindow`](@ref).
"""
mutable struct Impurity
    index :: I64
    atoms :: String
    sites :: I64
    shell :: String
    ising :: String
    occup :: F64
    upara :: F64
    jpara :: F64
    lpara :: F64
    beta  :: F64
    nband :: I64
end

"""
    PrTrait

Mutable struct. Essential information of a given projector.

### Members

* site -> Site in which the projector is defined.
* l    -> Quantum number 𝑙.
* m    -> Quantum number 𝑚.
* desc -> Projector's specification.

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

### Members

* site   -> Site in which the projectors are defined. In principle, the
            projectors included in the same group should be defined at
            the same site (or equivalently atom).
* l      -> Quantum number 𝑙. In principle, the projectors included in
            the same group should have the same quantum number 𝑙 (but
            with different 𝑚).
* corr   -> Test if the projectors in this group are correlated.
* shell  -> Type of correlated orbitals. It is infered from quantum number 𝑙.
* Pr     -> Array. It contains the indices of projectors.
* Tr     -> Array. It contains the transformation matrix. This parameter
            could be useful to select certain subset of orbitals or perform
            a simple global rotation.

See also: [`PrTrait`](@ref), [`PrWindow`](@ref), [`Mapping`](@ref), [`Impurity`](@ref).
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

Mutable struct. Define the band window for a given group of projectors.

### Members

* bmin -> Minimum band index.
* bmax -> Maximum band index.
* nbnd -> Maximum number of bands in the current window (≡ `bmax-bmin+1`).
* kwin -> Momentum-dependent and spin-dependent band window.
* bwin -> Tuple. It is the band window or energy window, which is used
          to filter the Kohn-Sham band structure. The mesh for calculating
          density of states is also deduced from `bwin`.

See also: [`PrTrait`](@ref), [`PrGroup`](@ref), [`Mapping`](@ref), [`Impurity`](@ref).
"""
mutable struct PrWindow
    bmin  :: I64
    bmax  :: I64
    nbnd  :: I64
    kwin  :: Array{I64,3}
    bwin  :: Tuple{R64,R64}
end

#=
### *Customized Constructors*
=#

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
    IterInfo()

Outer constructor for IterInfo struct.
"""
function IterInfo()
    # Extract the parameter `nsite` and `niter`
    nsite = get_i("nsite")
    _M₃, _M₁, _M₂ = get_m("niter")
    @assert nsite ≥ 1
    @assert _M₃ ≥ 1 && _M₁ ≥ 1 && _M₂ ≥ 1

    # Initialize key fields
    #
    # Note that sc = 1 means one-shot DFT + DMFT calculations,
    # while sc = 2 means fully self-consistent DFT + DMFT calculations.
    I  = 0
    M₁ = _M₁
    M₂ = _M₂
    M₃ = _M₃
    sc = 1
    μ  = 0.0
    dc = fill(0.0, nsite)
    n₁ = 0.0
    n₂ = 0.0
    nf = fill(0.0, nsite)
    et = 0.0
    cc = false
    ce = false
    cs = false

    # Call the default constructor
    IterInfo(I, I, I, I,
             M₁, M₂, M₃,
             sc,
             μ, μ, μ,
             dc,
             n₁, n₂, nf,
             et,
             cc, ce, cs)
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
    Mapping(nsite::I64, ngrp::I64, nwnd::I64)

Outer constructor for Mapping struct.
"""
function Mapping(nsite::I64, ngrp::I64, nwnd::I64)
    # Sanity check
    @assert ngrp ≥ nsite
    @assert nwnd ≤ ngrp

    # Initialize the arrays
    i_grp = zeros(I64, nsite)
    i_wnd = zeros(I64, nsite)
    g_imp = zeros(I64, ngrp)
    w_imp = zeros(I64, nwnd)

    # Call the default constructor
    Mapping(i_grp, i_wnd, g_imp, w_imp)
end

"""
    Impurity(index::I64, ...)

Outer constructor for Impurity struct.
"""
function Impurity(index::I64,
                  atoms::String, sites::I64, shell::String, ising::String,
                  occup::F64, upara::F64, jpara::F64, lpara::F64, beta::F64)
    # Define the mapping between `shell` and number of orbitals
    shell_to_dim = Dict{String,I64}(
                 "s"     => 1,
                 "p"     => 3,
                 "d"     => 5,
                 "f"     => 7,
                 "d_t2g" => 3, # Only a subset of d orbitals
                 "d_eg"  => 2, # Only a subset of d orbitals
             )

    # Determine number of orbitals of the quantum impurity problem
    nband = shell_to_dim[shell]

    # Call the default constructor
    Impurity(index, atoms, sites, shell, ising, occup, upara, jpara, lpara, beta, nband)
end

#=
*Remarks*:

Please go to the following webpage for more details about the original
specifications of projectors in the `vasp` code:
* https://www.vasp.at/wiki/index.php/LOCPROJ

May need to be fixed for the other DFT codes.
=#

"""
    PrTrait(site::I64, sort::String, desc::String)

Outer constructor for PrTrait struct.
"""
function PrTrait(site::I64, desc::String)
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

#=
### *Customized Base.show() Functions*
=#

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
    println(io, "I₁ : ", it.I₁)
    println(io, "I₂ : ", it.I₂)
    println(io, "I₃ : ", it.I₃)
    println(io, "I₄ : ", it.I₄)
    println(io, "M₁ : ", it.M₁)
    println(io, "M₂ : ", it.M₂)
    println(io, "M₃ : ", it.M₃)
    println(io, "sc : ", it.sc)
    println(io, "μ₀ : ", it.μ₀)
    println(io, "μ₁ : ", it.μ₁)
    println(io, "μ₂ : ", it.μ₂)
    println(io, "dc : ", it.dc)
    println(io, "n₁ : ", it.n₁)
    println(io, "n₂ : ", it.n₂)
    println(io, "nf : ", it.nf)
    println(io, "et : ", it.et)
    println(io, "cc : ", it.cc)
    println(io, "ce : ", it.ce)
    println(io, "cs : ", it.cs)
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
    println(io, "i_wnd : ", map.i_wnd)
    println(io, "g_imp : ", map.g_imp)
    println(io, "w_imp : ", map.w_imp)
end

"""
    Base.show(io::IO, imp::Impurity)

Base.show() function for Impurity struct.

See also: [`Impurity`](@ref).
"""
function Base.show(io::IO, imp::Impurity)
    println(io, "Impurity struct")
    println(io, "index : ", imp.index)
    println(io, "atoms : ", imp.atoms)
    println(io, "sites : ", imp.sites)
    println(io, "shell : ", imp.shell)
    println(io, "ising : ", imp.ising)
    println(io, "occup : ", imp.occup)
    println(io, "upara : ", imp.upara)
    println(io, "jpara : ", imp.jpara)
    println(io, "lpara : ", imp.lpara)
    println(io, "beta  : ", imp.beta)
    println(io, "nband : ", imp.nband)
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
