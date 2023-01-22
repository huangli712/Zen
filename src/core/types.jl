#
# Project : Pansy
# Source  : types.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/11/20
#

#=
### *Customized Types*
=#

"Customized types. It is used to define the following dicts."
const DType = Any

"Customized types. It is used to define the following dicts."
const ADT = Array{DType,1}

#=
### *Customized Dictionaries*
=#

#=
*Remarks* :

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
const PCASE  = Dict{String,ADT}(
    "case"     => [missing, 1, :String, "System's name or seedname"]
)

"""
    PDFT

Dictionary for configuration parameters: density functional theory
calculations. Its entries can not be updated once they are initialized.

See also: [`PCASE`](@ref), [`PDMFT`](@ref), [`PIMP`](@ref), [`PSOLVER`](@ref).
"""
const PDFT   = Dict{String,ADT}(
    "engine"   => [missing, 1, :String, "Engine for density functional theory calculations"],
    "projtype" => [missing, 1, :String, "Types of projectors"],
    "smear"    => [missing, 1, :String, "Scheme for smearing"],
    "kmesh"    => [missing, 1, :String, "K-mesh for brillouin zone sampling / integration"],
    "magmom"   => [missing, 0, :String, "Initial magnetic moments"],
    "ncycle"   => [missing, 1, :I64   , "Number of DFT iterations per DFT + DMFT cycle"],
    "lsymm"    => [missing, 1, :Bool  , "The symmetry is turned on or off"],
    "lspins"   => [missing, 1, :Bool  , "The spin orientations are polarized or not"],
    "lspinorb" => [missing, 1, :Bool  , "The spin-orbit coupling is considered or not"],
    "lproj"    => [missing, 1, :Bool  , "The projectors are generated or not"],
    "sproj"    => [missing, 1, :Array , "Specifications for generating projectors"],
    "window"   => [missing, 1, :Array , "Band / energy window for normalizing projectors"],
)

"""
    PDMFT

Dictionary for configuration parameters: dynamical mean-field theory
calculations. Be careful, its entries (such as `dcount`, `mixer`, `mc`,
`cc`, `ec`, `sc`, and `lfermi`) can be updated during the iteration.

See also: [`PCASE`](@ref), [`PDFT`](@ref), [`PIMP`](@ref), [`PSOLVER`](@ref).
"""
const PDMFT  = Dict{String,ADT}(
    "mode"     => [missing, 1, :I64   , "Scheme of dynamical mean-field theory calculations"],
    "axis"     => [missing, 1, :I64   , "Imaginary-time axis or real-frequency axis"],
    "niter"    => [missing, 1, :I64   , "Maximum allowed number of DFT + DMFT iterations"],
    "nmesh"    => [missing, 1, :I64   , "Number of frequency points"],
    "dcount"   => [missing, 1, :String, "Scheme of double counting term"],
    "beta"     => [missing, 1, :F64   , "Inverse system temperature"],
    "mixer"    => [missing, 1, :F64   , "Mixing factor"],
    "mc"       => [missing, 0, :F64   , "Convergence criterion of chemical potential"],
    "cc"       => [missing, 0, :F64   , "Convergence criterion of charge"],
    "ec"       => [missing, 0, :F64   , "Convergence criterion of total energy"],
    "sc"       => [missing, 0, :F64   , "Convergence criterion of self-energy function"],
    "lfermi"   => [missing, 0, :Bool  , "Whether chemical potential should be updated"],
)

"""
    PIMP

Dictionary for configuration parameters: quantum impurity problems. Be
careful, its entries (such as `ising`, `occup`, `upara`, `jpara`, and
`lpara`) can be updated during the iteration. 

See also: [`PCASE`](@ref), [`PDFT`](@ref), [`PDMFT`](@ref), [`PSOLVER`](@ref).
"""
const PIMP   = Dict{String,ADT}(
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

Dictionary for configuration parameters: quantum impurity solvers. Be
careful, its entry (only `params`) can be updated during the iteration.

See also: [`PCASE`](@ref), [`PDFT`](@ref), [`PDMFT`](@ref), [`PIMP`](@ref).
"""
const PSOLVER= Dict{String,ADT}(
    "engine"   => [missing, 1, :String, "Name of quantum impurity solver"],
    "ncycle"   => [missing, 1, :I64   , "Number of solver iterations per DFT + DMFT cycle"],
    "params"   => [missing, 1, :Array , "Extra parameter sets of quantum impurity solver"],
)

#=
### *Customized Structs* : *DFT Engine*
=#

"""
    AbstractEngine

An abstract type representing the DFT engine. It is used to build the
internal type system.
"""
abstract type AbstractEngine end

"""
    NULLEngine

It represents a null DFT engine. It is the default engine and will be
replaced by the realistic engine.

See also: [`_engine_`](@ref).
"""
struct NULLEngine    <: AbstractEngine end

"""
    VASPEngine

It represents a vasp engine, which is used to perform DFT calculations.

See also: [`QEEngine`](@ref).
"""
struct VASPEngine    <: AbstractEngine end

"""
    QEEngine

It represents a quantum espresso (actually the PWSCF code) engine, which
is used to perform DFT calculations.

See also: [`VASPEngine`](@ref).
"""
struct QEEngine      <: AbstractEngine end

"""
    WANNIEREngine

It represents a wannier90 engine, which is used to generate the maximally
localized wannier functions.
"""
struct WANNIEREngine <: AbstractEngine end

"Set up the default density functional theory calculation engine."
_engine_ = NULLEngine()

"Get name of density functional theory calculation engine."
Base.nameof(::NULLEngine) = "null"
Base.nameof(::VASPEngine) = "vasp"
Base.nameof(::QEEngine) = "qe"
Base.nameof(::WANNIEREngine) = "wan90"

#=
### *Customized Structs* : *Quantum Impurity Solver*
=#

"""
    AbstractSolver

An abstract type representing the quantum impurity solver. It is used to
build the internal type system.
"""
abstract type AbstractSolver end

"""
    NULLSolver

It represents a null quantum impurity solver. It is the default solver
and will be replaced by the realistic solver.

See also: [`_solver_`](@ref).
"""
struct NULLSolver   <: AbstractSolver end

"""
    CTHYB₁Solver

It represents a quantum impurity solver based on the CTHYB₁ algorithm.

See also: [`CTHYB₂Solver`](@ref).
"""
struct CTHYB₁Solver <: AbstractSolver end

"""
    CTHYB₂Solver

It represents a quantum impurity solver based on the CTHYB₂ algorithm.

See also: [`CTHYB₁Solver`](@ref).
"""
struct CTHYB₂Solver <: AbstractSolver end

"""
    HIASolver

It represents a quantum impurity solver based on the HIA algorithm.
"""
struct HIASolver    <: AbstractSolver end

"""
    NORGSolver

It represents a quantum impurity solver based on the NORG algorithm.
"""
struct NORGSolver   <: AbstractSolver end

"""
    ATOMSolver

It represents a solver for atomic eigenvalue problems.
"""
struct ATOMSolver   <: AbstractSolver end

"Set up the default quantum impurity solver."
_solver_ = NULLSolver()

"Get name of quantum impurity solver."
Base.nameof(::NULLSolver) = "null"
Base.nameof(::CTHYB₁Solver) = "ct_hyb1"
Base.nameof(::CTHYB₂Solver) = "ct_hyb2"
Base.nameof(::HIASolver) = "hia"
Base.nameof(::NORGSolver) = "norg"
Base.nameof(::ATOMSolver) = "atomic"

#=
### *Customized Structs* : *Adaptor*
=#

"""
    AbstractAdaptor

An abstract type representing the DFT-DMFT adaptor. It is used to build
the internal type system.
"""
abstract type AbstractAdaptor end

"""
    NULLAdaptor

It represents a null DFT-DMFT adaptor. It is the default adaptor and will
be replaced by the realistic adaptor.

See also: [`_adaptor_`](@ref).
"""
struct NULLAdaptor    <: AbstractAdaptor end

"""
    PLOAdaptor

It represents a DFT-DMFT adaptor based on the projected local orbitals.

See also: [`WANNIERAdaptor`](@ref).
"""
struct PLOAdaptor     <: AbstractAdaptor end

"""
    PLOAdaptor

It represents a DFT-DMFT adaptor based on the maximally localized
wannier functions.

See also: [`PLOAdaptor`](@ref).
"""
struct WANNIERAdaptor <: AbstractAdaptor end

"Set up the default DFT-DMFT adaptor."
_adaptor_ = NULLAdaptor()

"Get name of DFT-DMFT adaptor."
Base.nameof(::NULLAdaptor) = "null"
Base.nameof(::PLOAdaptor) = "plo"
Base.nameof(::WANNIERAdaptor) = "wannier"

#=
### *Customized Structs* : *Sigma Engine*
=#

"""
    AbstractMode

An abstract type representing various operations on the self-energy
functions. It is used to build the internal type system.
"""
abstract type AbstractMode end

"""
    NULLMode

It represents a null operation. It is the default operation and will
be replaced by the realistic operation.

See also: [`_mode_`](@ref).
"""
struct NULLMode   <: AbstractMode end

"""
    RESETMode

It represents a reset operation for resetting the self-energy functions.
"""
struct RESETMode  <: AbstractMode end

"""
    DCOUNTMode

It represents a dcount operation for creating the double counting terms.
"""
struct DCOUNTMode <: AbstractMode end

"""
    SPLITMode

It represents a split operation for splitting the hybridization functions.
"""
struct SPLITMode  <: AbstractMode end

"""
    GATHERMode

It represents a gather operation for gathering the self-energy functions.
"""
struct GATHERMode <: AbstractMode end

"Set up the default operation."
_mode_ = NULLMode()

"Get name of operation."
Base.nameof(::NULLMode) = "null"
Base.nameof(::RESETMode) = "reset"
Base.nameof(::DCOUNTMode) = "dcount"
Base.nameof(::SPLITMode) = "split"
Base.nameof(::GATHERMode) = "gather"

#=
### *Customized Structs* : *Mixer Engine*
=#

"""
    AbstractMixer

An abstract type representing various mixers for self-energy functions Σ,
hybridization function Δ, effective impurity level ϵ, and density matrix
Γ. It is used to build the internal type system.
"""
abstract type AbstractMixer end

"""
    NULLMixer

It represents a null mixer. It is the default mixer and will be replaced
by the realistic mixer.

See also: [`_mixer_`](@ref).
"""
struct NULLMixer <: AbstractMixer end

"""
    ΣMixer

It represents a mixer for self-energy function Σ.
"""
struct ΣMixer    <: AbstractMixer end

"""
    ΔMixer

It represents a mixer for hybridization function Δ.
"""
struct ΔMixer    <: AbstractMixer end

"""
    EMixer

It represents a mixer for effective impurity level ϵ.
"""
struct EMixer    <: AbstractMixer end

"""
    ΓMixer

It represents a mixer for density matrix Γ.
"""
struct ΓMixer    <: AbstractMixer end

"Set up the default mixer."
_mixer_ = NULLMixer()

"Get name of mixer."
Base.nameof(::NULLMixer) = "null"
Base.nameof(::ΣMixer) = "Σ"
Base.nameof(::ΔMixer) = "Δ"
Base.nameof(::EMixer) = "E"
Base.nameof(::ΓMixer) = "Γ"

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
    Energy

Mutable struct. Store decomposition of the total DFT + DMFT energy.

### Members

* dft  -> DFT band energy.
* dmft -> DMFT interaction energy Tr(ΣG).
* corr -> DMFT Correction to the DFT band energy.
* dc   -> Energy contributed by the double counting term.

See also: [`IterInfo`](@ref).
"""
mutable struct Energy
    dft  :: F64
    dmft :: F64
    corr :: F64
    dc   :: F64
end

"""
    IterInfo

Mutable struct. Record the DFT + DMFT iteration information.

### Members

* I₁ -> Number of cycles between `dmft1` and quantum impurity solver.
* I₂ -> Number of cycles between `dmft2` and DFT engine.
* I₃ -> Number of DFT + DMFT iterations.
* I₄ -> Counter for all the internal cycles.
* M₁ -> Maximum allowed number of cycles (between `dmft1` and solver).
* M₂ -> Maximum allowed number of cycles (between `dmft2` and DFT).
* M₃ -> Maximum allowed number of DFT + DMFT iterations.
* sc -> Self-consistent mode.
* μ₀ -> Fermi level obtained by DFT engine.
* μ₁ -> Fermi level obtained by DMFT engine (`dmft1`).
* μ₂ -> Fermi level obtained by DMFT engine (`dmft2`).
* dc -> Double counting terms.
* n₁ -> Number of lattice occupancy obtained by DMFT engine (`dmft1`).
* n₂ -> Number of lattice occupancy obtained by DMFT engine (`dmft2`).
* nf -> Number of impurity occupancy obtained by quantum impurity solver.
* et -> Total DFT + DMFT energy (for current iteration).
* ep -> Total DFT + DMFT energy (for previous iteration).
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
    dc :: Vector{F64} # dc is site-dependent
    n₁ :: F64
    n₂ :: F64
    nf :: Vector{F64} # nf is site-dependent
    et :: Energy
    ep :: Energy
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

Mutable struct. It denotes a mapping between quantum impurity problems
and groups of projectors (or windows of Kohn-Sham states).

### Members

* i_grp -> Mapping from quntum impurity problems to groups of projectors.
* i_wnd -> Mapping from quantum impurity problems to windows of Kohn-Sham states.
* g_imp -> Mapping from groups of projectors to quantum impurity problems.
* w_imp -> Mapping from windows of Kohn-Sham states to quantum impurity problems.

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

Mutable struct. It includes key information of quantum impurity problems.

### Members

* index -> Index of the quantum impurity problem.
* atoms -> Chemical symbol of impurity atom.
* sites -> Index of impurity atom.
* equiv -> Equivalence of quantum impurity problem.
* shell -> Angular momentum of correlated orbitals.
* ising -> Interaction type of correlated orbitals.
* occup -> Impurity occupancy 𝑛.
* nup   -> Impurity occupancy 𝑛↑ (spin up).
* ndown -> Impurity occupancy 𝑛↓ (spin down).
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
    equiv :: I64
    shell :: String
    ising :: String
    occup :: F64
    nup   :: F64
    ndown :: F64
    upara :: F64
    jpara :: F64
    lpara :: F64
    beta  :: F64
    nband :: I64
end

"""
    PrTrait

Mutable struct. It defines a local orbital projection (a projector).

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

Mutable struct. It defines a group of projectors. There are quite a lot
of local projectors. We always gather, reorganize, and split them into
various groups according to their labels (such as site and 𝑙). Each group
is connected with a quantum impurity problem.

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

Mutable struct. It defines a window for the Kohn-Sham states (DFT bands).
Each window is connected with a quantum impurity problem.

### Members

* bmin -> Minimum band index.
* bmax -> Maximum band index.
* nbnd -> Maximum number of bands in the current window (≡ `bmax - bmin + 1`).
* kwin -> Momentum-dependent and spin-dependent band window.
* bwin -> Tuple. It is the band window or energy window, which is used
          to filter the Kohn-Sham states (i.e DFT bands). The mesh for
          calculating density of states is also deduced from `bwin`.

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

    # Sanity check
    @assert isopen(log) && isopen(cycle)

    # Call the default constructor
    Logger(log, cycle)
end

"""
    Energy()

Outer constructor for Energy struct.
"""
function Energy()
    # Set the decomposition of energy to zero
    dft  = 0.0
    dmft = 0.0
    corr = 0.0
    dc   = 0.0

    # Call the default constructor
    Energy(dft, dmft, corr, dc)
end

"""
    IterInfo()

Outer constructor for IterInfo struct.
"""
function IterInfo()
    # Extract the parameter `nsite`
    nsite = get_i("nsite")
    @assert nsite ≥ 1

    # Extract the parameters `niter` and `ncycle`
    _M₁ = get_s("ncycle")
    _M₂ = get_d("ncycle")
    _M₃ = get_m("niter")
    @assert _M₁ ≥ 1 && _M₂ ≥ 1 && _M₃ ≥ 1

    # Initialize key fields
    #
    # sc = 0 means in preparation stage;
    # sc = 1 means one-shot DFT + DMFT calculations;
    # sc = 2 means fully self-consistent DFT + DMFT calculations.
    I  = 0
    M₁ = _M₁
    M₂ = _M₂
    M₃ = _M₃
    sc = 0
    μ  = 0.0
    dc = fill(0.0, nsite)
    n₁ = 0.0
    n₂ = 0.0
    nf = fill(0.0, nsite)
    et = Energy()
    ep = Energy()
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
             et, ep,
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
    @assert nwnd == ngrp ≥ nsite

    # Initialize the arrays
    i_grp = zeros(I64, nsite)
    i_wnd = zeros(I64, nsite)
    g_imp = zeros(I64, ngrp)
    w_imp = zeros(I64, nwnd)

    # Call the default constructor
    Mapping(i_grp, i_wnd, g_imp, w_imp)
end

"""
    Impurity(index::I64, atoms::String, sites::I64,
             equiv::I64, shell::String, ising::String,
             occup::F64, upara::F64, jpara::F64, lpara::F64,
             beta::F64)

Outer constructor for Impurity struct.
"""
function Impurity(index::I64, atoms::String, sites::I64,
                  equiv::I64, shell::String, ising::String,
                  occup::F64, upara::F64, jpara::F64, lpara::F64,
                  beta::F64)
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

    # Determine impurity occupancy
    nup = occup / 2.0
    ndown = occup / 2.0

    # Call the default constructor
    Impurity(index, atoms, sites,
             equiv, shell, ising,
             occup, nup, ndown,
             upara, jpara, lpara, beta, nband)
end

#=
*Remarks* :

Please go to the following webpage for more details about the original
specifications of projectors in the `vasp` code:

* https://www.vasp.at/wiki/index.php/LOCPROJ

May need to be fixed for the other DFT codes.
=#

"""
    PrTrait(site::I64, desc::String)

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
    lm = findfirst(x -> x == desc, orb_labels) - 1
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
### *Customized Binary Operations*
=#

import Base.==

"""
    ==(PT₁::PrTrait, PT₂::PrTrait)

Compare two PrTrait objects.

See also: [`PrTrait`](@ref).
"""
function ==(PT₁::PrTrait, PT₂::PrTrait)
    PT₁.site == PT₂.site &&
    PT₁.l    == PT₂.l    &&
    PT₁.m    == PT₂.m    &&
    PT₁.desc == PT₂.desc
end

"""
    ==(PG₁::PrGroup, PG₂::PrGroup)

Compare two PrGroup objects.

See also: [`PrGroup`](@ref).
"""
function ==(PG₁::PrGroup, PG₂::PrGroup)
    PG₁.site  == PG₂.site  &&
    PG₁.l     == PG₂.l     &&
    PG₁.corr  == PG₂.corr  &&
    PG₁.shell == PG₂.shell &&
    PG₁.Pr    == PG₂.Pr    &&
    PG₁.Tr    == PG₂.Tr
end

"""
    ==(PW₁::PrWindow, PW₂::PrWindow)

Compare two PrWindow objects.

See also: [`PrWindow`](@ref).
"""
function ==(PW₁::PrWindow, PW₂::PrWindow)
    PW₁.bmin == PW₂.bmin &&
    PW₁.bmax == PW₂.bmax &&
    PW₁.nbnd == PW₂.nbnd &&
    PW₁.kwin == PW₂.kwin &&
    PW₁.bwin == PW₂.bwin
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
    println(io, repeat("=", 20))
    #
    println(io, "log   : ", logger.log)
    println(io, "cycle : ", logger.cycle)
    #
    println(io, repeat("=", 20))
end

"""
    Base.show(io::IO, ene::Energy)

Base.show() function for Energy struct. Note that `total` is not a real
field of the `Energy` struct.

See also: [`Energy`](@ref).
"""
function Base.show(io::IO, ene::Energy)
    println(io, "Energy struct")
    println(io, repeat("=", 20))
    #
    println(io, "dft   : ", ene.dft)
    println(io, "dmft  : ", ene.dmft)
    println(io, "corr  : ", ene.corr)
    println(io, "dc    : ", ene.dc)
    println(io, "total : ", ene.total)
    #
    println(io, repeat("=", 20))
end

"""
    Base.show(io::IO, it::IterInfo)

Base.show() function for IterInfo struct.

See also: [`IterInfo`](@ref).
"""
function Base.show(io::IO, it::IterInfo)
    println(io, "IterInfo struct")
    println(io, repeat("=", 20))
    #
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
    println(io, "ep : ", it.ep)
    println(io, "cc : ", it.cc)
    println(io, "ce : ", it.ce)
    println(io, "cs : ", it.cs)
    #
    println(io, repeat("=", 20))
end

"""
    Base.show(io::IO, latt::Lattice)

Base.show() function for Lattice struct.

See also: [`Lattice`](@ref).
"""
function Base.show(io::IO, latt::Lattice)
    println(io, "Lattice struct")
    println(io, repeat("=", 20))
    #
    println(io, "_case : ", latt._case)
    println(io, "scale : ", latt.scale)
    println(io, "lvect : ", latt.lvect)
    println(io, "nsort : ", latt.nsort)
    println(io, "natom : ", latt.natom)
    println(io, "sorts : ", latt.sorts)
    println(io, "atoms : ", latt.atoms)
    println(io, "coord : ", latt.coord)
    #
    println(io, repeat("=", 20))
end

"""
    Base.show(io::IO, map::Mapping)

Base.show() function for Mapping struct.

See also: [`Mapping`](@ref).
"""
function Base.show(io::IO, map::Mapping)
    println(io, "Mapping struct")
    println(io, repeat("=", 20))
    #
    println(io, "i_grp : ", map.i_grp)
    println(io, "i_wnd : ", map.i_wnd)
    println(io, "g_imp : ", map.g_imp)
    println(io, "w_imp : ", map.w_imp)
    #
    println(io, repeat("=", 20))
end

"""
    Base.show(io::IO, imp::Impurity)

Base.show() function for Impurity struct.

See also: [`Impurity`](@ref).
"""
function Base.show(io::IO, imp::Impurity)
    println(io, "Impurity struct")
    println(io, repeat("=", 20))
    #
    println(io, "index : ", imp.index)
    println(io, "atoms : ", imp.atoms)
    println(io, "sites : ", imp.sites)
    println(io, "equiv : ", imp.equiv)
    println(io, "shell : ", imp.shell)
    println(io, "ising : ", imp.ising)
    println(io, "occup : ", imp.occup)
    println(io, "nup   : ", imp.nup)
    println(io, "ndown : ", imp.ndown)
    println(io, "upara : ", imp.upara)
    println(io, "jpara : ", imp.jpara)
    println(io, "lpara : ", imp.lpara)
    println(io, "beta  : ", imp.beta)
    println(io, "nband : ", imp.nband)
    #
    println(io, repeat("=", 20))
end

"""
    Base.show(io::IO, PT::PrTrait)

Base.show() function for PrTrait struct.

See also: [`PrTrait`](@ref).
"""
function Base.show(io::IO, PT::PrTrait)
    println(io, "PrTrait struct")
    println(io, repeat("=", 20))
    #
    println(io, "site : ", PT.site)
    println(io, "l    : ", PT.l)
    println(io, "m    : ", PT.m)
    println(io, "desc : ", PT.desc)
    #
    println(io, repeat("=", 20))
end

"""
    Base.show(io::IO, PG::PrGroup)

Base.show() function for PrGroup struct.

See also: [`PrGroup`](@ref).
"""
function Base.show(io::IO, PG::PrGroup)
    println(io, "PrGroup struct")
    println(io, repeat("=", 20))
    #
    println(io, "site  : ", PG.site)
    println(io, "l     : ", PG.l)
    println(io, "corr  : ", PG.corr)
    println(io, "shell : ", PG.shell)
    println(io, "Pr    : ", PG.Pr)
    println(io, "Tr    : ", PG.Tr)
    #
    println(io, repeat("=", 20))
end

"""
    Base.show(io::IO, PW::PrWindow)

Base.show() function for PrWindow struct.

See also: [`PrWindow`](@ref).
"""
function Base.show(io::IO, PW::PrWindow)
    println(io, "PrWindow struct")
    println(io, repeat("=", 20))
    #
    println(io, "bmin : ", PW.bmin)
    println(io, "bmax : ", PW.bmax)
    println(io, "nbnd : ", PW.nbnd)
    println(io, "kwin : ", PW.kwin)
    println(io, "bwin : ", PW.bwin)
    #
    println(io, repeat("=", 20))
end

#=
### *Customized Base.getproperty() Functions*
=#

"""
    Base.getproperty(et::Energy, sym::Symbol)

Implement the calculation of total DFT + DMFT energy. The `Energy` struct
does not really contains the `total` field. This function will implement
`et.total` by overriding the Base.getproperty() function.

See also: [`Energy`](@ref).
"""
function Base.getproperty(et::Energy, sym::Symbol)
    if sym == :total
        return et.dft + et.dmft + et.corr - et.dc
    else
        return getfield(et, sym)
    end
end
