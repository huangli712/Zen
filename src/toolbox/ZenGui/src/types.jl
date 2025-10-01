#
# Project : Camellia
# Source  : types.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2025/07/16
#

#=
### *Customized Colors*
=#

const COL_MAGENTA    = ImVec4(1.00, 0.00, 1.00, 1.00)
const COL_ORANGE     = ImVec4(1.00, 0.50, 0.00, 1.00)
const COL_LIGHTGREEN = ImVec4(0.50, 0.90, 0.50, 1.00)
const COL_PINK       = ImVec4(1.00, 0.75, 0.80, 1.00)
const COL_TEAL       = ImVec4(0.00, 0.50, 0.50, 1.00)
const COL_INDIGO     = ImVec4(0.29, 0.00, 0.51, 1.00)
const COL_PURPLE     = ImVec4(0.50, 0.00, 0.50, 1.00)
const COL_SALMON     = ImVec4(0.98, 0.50, 0.45, 1.00)
const COL_TAN        = ImVec4(0.82, 0.71, 0.55, 1.00)
const COL_TOMATO     = ImVec4(1.00, 0.39, 0.28, 1.00)
const COL_HOTPINK    = ImVec4(1.00, 0.41, 0.71, 1.00)
const COL_SKYBLUE    = ImVec4(0.50, 0.80, 1.00, 1.00)

#=
### *Customized Structs* : *Active Window*
=#

"""
    CURRENT_WINDOW

A struct used to keep the name of the current (activate) window. Note that
its global instance is `CWIN`

### Members
* name -> Name of the current window.
"""
mutable struct CURRENT_WINDOW
    name :: String
end

"""
    CWIN

An instance for the `CURRENT_WINDOW` struct.

See also: [`CURRENT_WINDOW`](@ref).
"""
const CWIN = CURRENT_WINDOW(
    "nothing"
)

#=
### *Customized Structs* : *Menu Flags*
=#

"""
    MenuFlags

A struct used to track the status of all the menu items. The renderloop
should respond the mouse events according to this struct. Note that its
global instance is `FMENU`

The meun items are created at `src/menu.jl`.

### Members
* F_SAVE     -> File | Save.
* F_EXIT     -> File | Exit.
* E_ZEN      -> Edit | Integrated Package | Zen.
* E_DYSON    -> Edit | Quantum Many-Body Theory Engines | Dyson.
* E_DFERMION -> Edit | Quantum Many-Body Theory Engines | DFermion.
* E_CTSEG    -> Edit | Quantum Impurity Solvers | iQIST | ctseg.
* E_CTHYB    -> Edit | Quantum Impurity Solvers | iQIST | cthyb.
* E_ATOMIC   -> Edit | Quantum Impurity Solvers | iQIST | atomic.
* E_ACFLOW   -> Edit | Analytic Continuation Tools | ACFlow.
* E_ACTEST   -> Edit | Analytic Continuation Tools | ACTest.
* S_BGIMAGE  -> Style | Change Background.
* S_CLASSIC  -> Style | Classic.
* S_DARK     -> Style | Dark.
* S_LIGHT    -> Style | Light.
* H_ZEN      -> Help | Documentation | Zen.
* H_DYSON    -> Help | Documentation | Dyson.
* H_DFERMION -> Help | Documentation | DFermion.
* H_IQIST    -> Help | Documentation | iQIST.
* H_ACFLOW   -> Help | Documentation | ACFlow.
* H_ACTEST   -> Help | Documentation | ACTest.
* H_ZENGUI   -> Help | User's Manual.
* H_ABOUT    -> Help | About ZenGui.
"""
mutable struct MenuFlags
    F_SAVE     :: Bool
    F_EXIT     :: Bool
    #
    E_ZEN      :: Bool
    E_DYSON    :: Bool
    E_DFERMION :: Bool
    E_CTSEG    :: Bool
    E_CTHYB    :: Bool
    E_ATOMIC   :: Bool
    E_ACFLOW   :: Bool
    E_ACTEST   :: Bool
    #
    S_BGIMAGE  :: Bool
    S_CLASSIC  :: Bool
    S_DARK     :: Bool
    S_LIGHT    :: Bool
    #
    H_ZEN      :: Bool
    H_DYSON    :: Bool
    H_DFERMION :: Bool
    H_IQIST    :: Bool
    H_ACFLOW   :: Bool
    H_ACTEST   :: Bool
    H_ZENGUI   :: Bool
    H_ABOUT    :: Bool
end

"""
    FMENU

An instance for the `MenuFlags` struct. Initially, all members are set to
be false.

See also: [`MenuFlags`](@ref).
"""
const FMENU = MenuFlags(
    false, # F_SAVE
    false, # F_EXIT
    false, # E_ZEN
    false, # E_DYSON
    false, # E_DFERMION
    false, # E_CTSEG
    false, # E_CTHYB
    false, # E_ATOMIC
    false, # E_ACFLOW
    false, # E_ACTEST
    false, # S_BGIMAGE
    false, # S_CLASSIC
    false, # S_DARK
    false, # S_LIGHT
    false, # H_ZEN
    false, # H_DYSON
    false, # H_DFERMION
    false, # H_IQIST
    false, # H_ACFLOW
    false, # H_ACTEST
    false, # H_ZENGUI
    false  # H_ABOUT
)

#=
### *Customized Structs* : *Zen Package*
=#

#=
*Remarks* :

The official configuration file for the Zen package is `case.toml`. It
includes `[case]`, `[dft]`, `[dmft]`, `[impurity]`, and `[solver]` blocks.
=#

"""
    ZEN_PCASE

This struct represents the `[case]` block in the `case.toml` file.
"""
mutable struct ZEN_PCASE
    case :: String
end

"""
    ZEN_PDFT

This struct represents the `[dft]` block in the `case.toml` file.
"""
mutable struct ZEN_PDFT
    engine   :: String
    projtype :: String
    smear    :: String
    kmesh    :: String
    magmom   :: String
    ncycle   :: I64
    lsymm    :: Bool
    lspins   :: Bool
    lspinorb :: Bool
    lproj    :: Bool
    sproj    :: Array
    window   :: Array
end

"""
    ZEN_PDMFT

This struct represents the `[dmft]` block in the `case.toml` file.
"""
mutable struct ZEN_PDMFT
    mode   :: I64
    axis   :: I64
    niter  :: I64
    nmesh  :: I64
    dcount :: String
    beta   :: F64
    mixer  :: F64
    mc     :: F64
    cc     :: F64
    ec     :: F64
    sc     :: F64
    lfermi :: Bool
end

"""
    ZEN_PIMPURITY

This struct represents the `[impurity]` block in the `case.toml` file.
"""
mutable struct ZEN_PIMPURITY
    nsite :: I64
    atoms :: Array
    equiv :: Array
    shell :: Array
    ising :: Array
    occup :: Array
    upara :: Array
    jpara :: Array
    lpara :: Array
end

"""
    ZEN_PSOLVER

This struct represents the `[solver]` block in the `case.toml` file.
"""
mutable struct ZEN_PSOLVER
    engine :: String
    ncycle :: I64
    params :: Array
end

"""
    PCASE

An instance for the `ZEN_PCASE` struct.

See also: [`ZEN_PCASE`](@ref).
"""
PCASE = ZEN_PCASE(
    "SrVO3" # case
)

"""
    PDFT

An instance for the `ZEN_PDFT` struct.

See also: [`ZEN_PDFT`](@ref).
"""
PDFT = ZEN_PDFT(
    "vasp",         # engine
    "plo",          # projtype
    "tetra",        # smear
    "medium",       # kmesh
    "1.0",          # magmom
    8,              # ncycle
    false,          # lsymm
    false,          # lspins
    false,          # lspinorb
    true,           # lproj
    ["2 : d : Pr"], # sproj
    [-1.4, 6.0]     # window
)

"""
    PDMFT

An instance for the `ZEN_PDMFT` struct.

See also: [`ZEN_PDMFT`](@ref).
"""
PDMFT = ZEN_PDMFT(
    1,      # mode
    1,      # axis
    60,     # niter
    8193,   # nmesh
    "fll2", # dcount
    40.0,   # beta
    0.1,    # mixer
    1.0e-4, # mc
    1.0e-6, # cc
    1.0e-4, # ec
    1.0e-4, # sc
    true    # lfermi
)

"""
    PIMPURITY

An instance for the `ZEN_PIMPURITY` struct.

See also: [`ZEN_PIMPURITY`](@ref).
"""
PIMPURITY = ZEN_PIMPURITY(
    1,         # nsite
    ["V : 2"], # atoms
    [1],       # equiv
    ["d"],     # shell
    ["ising"], # ising
    [1.0],     # occup
    [4.0],     # upara
    [0.7],     # jpara
    [0.0]      # lpara
)

"""
    PSOLVER

An instance for the `ZEN_PSOLVER` struct.

See also: [`ZEN_PSOLVER`](@ref).
"""
PSOLVER = ZEN_PSOLVER(
    "ctseg", # engine
    2,       # ncycle
    ["isbnd = 2", "isort = 2", "nsweep = 100000000"] # params
)

"""
    struct_to_dict(s::ZEN_PCASE)

Convert a struct to an ordered dictionary (for `ZEN_PCASE`).

See [`ZEN_PCASE`](@ref).
"""
function struct_to_dict(s::ZEN_PCASE)
    return OrderedDict{String,Any}(
        "case" => s.case,
    )
end

"""
    struct_to_dict(s::ZEN_PDFT)

Convert a struct to an ordered dictionary (for `ZEN_PDFT`).

See [`ZEN_PDFT`](@ref).
"""
function struct_to_dict(s::ZEN_PDFT)
    return OrderedDict{String,Any}(
        "engine"   => s.engine,
        "projtype" => s.projtype,
        "smear"    => s.smear,
        "kmesh"    => s.kmesh,
        "magmom"   => s.magmom,
        "ncycle"   => s.ncycle,
        "lsymm"    => s.lsymm,
        "lspins"   => s.lspins,
        "lspinorb" => s.lspinorb,
        "lproj"    => s.lproj,
        "sproj"    => s.sproj,
        "window"   => s.window,
    )
end

"""
    struct_to_dict(s::ZEN_PDMFT)

Convert a struct to an ordered dictionary (for `ZEN_PDMFT`).

See [`ZEN_PDMFT`](@ref).
"""
function struct_to_dict(s::ZEN_PDMFT)
    return OrderedDict{String,Any}(
        "mode"     => s.mode,
        "axis"     => s.axis,
        "niter"    => s.niter,
        "nmesh"    => s.nmesh,
        "dcount"   => s.dcount,
        "beta"     => s.beta,
        "mixer"    => s.mixer,
        "mc"       => s.mc,
        "cc"       => s.cc,
        "ec"       => s.ec,
        "sc"       => s.sc,
        "lfermi"   => s.lfermi,
    )
end

"""
    struct_to_dict(s::ZEN_PIMPURITY)

Convert a struct to an ordered dictionary (for `ZEN_PIMPURITY`).

See [`ZEN_PIMPURITY`](@ref).
"""
function struct_to_dict(s::ZEN_PIMPURITY)
    return OrderedDict{String,Any}(
        "nsite"    => s.nsite,
        "atoms"    => s.atoms,
        "equiv"    => s.equiv,
        "shell"    => s.shell,
        "ising"    => s.ising,
        "occup"    => s.occup,
        "upara"    => s.upara,
        "jpara"    => s.jpara,
        "lpara"    => s.lpara,
    )
end

"""
    struct_to_dict(s::ZEN_PSOLVER)

Convert a struct to an ordered dictionary (for `ZEN_PSOLVER`).

See [`ZEN_PSOLVER`](@ref).
"""
function struct_to_dict(s::ZEN_PSOLVER)
    return OrderedDict{String,Any}(
        "engine"   => s.engine,
        "ncycle"   => s.ncycle,
        "params"   => s.params,
    )
end

"""
    build_zen_dict()

Assemble the ordered dictionary, which is then converted into `case.toml`,
for the `Zen` package.
"""
function build_zen_dict()
    return OrderedDict{String,Any}(
        "case" => struct_to_dict(PCASE)["case"],
        "dft" => struct_to_dict(PDFT),
        "dmft" => struct_to_dict(PDMFT),
        "impurity" => struct_to_dict(PIMPURITY),
        "solver" => struct_to_dict(PSOLVER)
    )
end

#=
### *Customized Structs* : *Dyson Code*
=#

#=
*Remarks* :

The official configuration file for the Dyson code is `dmft.in`, which is
actually a ini-like file.
=#

"""
    DYSON_PDYSON

This struct encapsulates the parameters in the `dmft.in` file.

See also: [`PDYSON`](@ref).
"""
mutable struct DYSON_PDYSON
    task   :: I64
    axis   :: I64
    beta   :: F64
    mc     :: F64
    lfermi :: String
    ltetra :: String
end

"""
    _DYSON

This set records the names of modified parameters, which will be presented
in the `dmft.in` file. Note that not all the parameters in `DYSON_PDYSON`
should be presented in the `dmft.in` file.

See also: [`DYSON_PDYSON`](@ref).
"""
_DYSON = Set{String}()

"""
    PDYSON

An instance for the `DYSON_PDYSON` struct.

See also: [`DYSON_PDYSON`](@ref).
"""
PDYSON = DYSON_PDYSON(
    1,        # task
    1,        # axis
    8.0,      # beta
    0.0001,   # mc
    ".true.", # lfermi
    ".true."  # ltetra
)

"""
    struct_to_dict(s::DYSON_PDYSON)

Convert a struct to an ordered dictionary (for `DYSON_PDYSON`).

See also: [`DYSON_PDYSON`](@ref).
"""
function struct_to_dict(s::DYSON_PDYSON)
    OD = OrderedDict{String,Any}()
    #
    "task"   ∈ _DYSON && ( OD["task"]   = s.task   )
    "axis"   ∈ _DYSON && ( OD["axis"]   = s.axis   )
    "beta"   ∈ _DYSON && ( OD["beta"]   = s.beta   )
    "mc"     ∈ _DYSON && ( OD["mc"]     = s.mc     )
    "lfermi" ∈ _DYSON && ( OD["lfermi"] = s.lfermi )
    "ltetra" ∈ _DYSON && ( OD["ltetra"] = s.ltetra )
    #
    return OD
end

"""
    build_dyson_dict()

Assemble the ordered dictionary, which is then converted into `dmft.in`,
for the `Dyson` code.

See also: [`DYSON_PDYSON`](@ref).
"""
function build_dyson_dict()
    return struct_to_dict(PDYSON)
end

#=
### *Customized Structs* : *DFermion Code*
=#

#=
*Remarks* :

The official configuration file for the DFermion code is `dfa.in`, which
is actually a ini-like file.
=#

"""
    DFERMION_PDFERMION

This struct encapsulates the parameters in the `dfa.in` file.

See also: [`PDFERMION`](@ref).
"""
mutable struct DFERMION_PDFERMION
    isdia :: I64 # Cycle
    nband :: I64 # Model
    nspin :: I64 # Model
    norbs :: I64 # Model
    nffrq :: I64 # Dimension
    nbfrq :: I64 # Dimension
    nkpts :: I64 # K-mesh
    nkp_x :: I64 # K-mesh
    nkp_y :: I64 # K-mesh
    nkp_z :: I64 # K-mesh
    ndfit :: I64 # Cycle
    nbsit :: I64 # Cycle
    mune  :: F64 # Model
    beta  :: F64 # Model
    part  :: F64 # Model
    dfmix :: F64 # Cycle
    bsmix :: F64 # Cycle
end

"""
    _DFERMION

This set records the names of modified parameters, which will be presented
in the `dfa.in` file. Note that not all the parameters in `DFERMION_PDFERMION`
should be presented in the `dfa.in` file.

See also: [`DFERMION_PDFERMION`](@ref).
"""
_DFERMION = Set{String}()

"""
    PDFERMION

An instance for the `DFERMION_PDFERMION` struct.

See also: [`DFERMION_PDFERMION`](@ref).
"""
PDFERMION = DFERMION_PDFERMION(
    2,   # isdia
    1,   # nband
    2,   # nspin
    2,   # norbs
    16,  # nffrq
    7,   # nbfrq
    64,  # nkpts
    8,   # nkp_x
    8,   # nkp_y
    8,   # nkp_z
    10,  # ndfit
    10,  # nbsit
    0.0, # mune
    1.0, # beta
    1.0, # part
    1.0, # dfmix
    0.7  # bsmix
)

"""
    struct_to_dict(s::DFERMION_PDFERMION)

Convert a struct to an ordered dictionary (for `DFERMION_PDFERMION`).

See also: [`DFERMION_PDFERMION`](@ref).
"""
function struct_to_dict(s::DFERMION_PDFERMION)
    OD = OrderedDict{String,Any}()
    #
    "isdia" ∈ _DFERMION && ( OD["isdia"] = s.isdia )
    "nband" ∈ _DFERMION && ( OD["nband"] = s.nband )
    "nspin" ∈ _DFERMION && ( OD["nspin"] = s.nspin )
    "norbs" ∈ _DFERMION && ( OD["norbs"] = s.norbs )
    "nffrq" ∈ _DFERMION && ( OD["nffrq"] = s.nffrq )
    "nbfrq" ∈ _DFERMION && ( OD["nbfrq"] = s.nbfrq )
    "nkpts" ∈ _DFERMION && ( OD["nkpts"] = s.nkpts )
    "nkp_x" ∈ _DFERMION && ( OD["nkp_x"] = s.nkp_x )
    "nkp_y" ∈ _DFERMION && ( OD["nkp_y"] = s.nkp_y )
    "nkp_z" ∈ _DFERMION && ( OD["nkp_z"] = s.nkp_z )
    "ndfit" ∈ _DFERMION && ( OD["ndfit"] = s.ndfit )
    "nbsit" ∈ _DFERMION && ( OD["nbsit"] = s.nbsit )
    "mune"  ∈ _DFERMION && ( OD["mune"]  = s.mune  )
    "beta"  ∈ _DFERMION && ( OD["beta"]  = s.beta  )
    "part"  ∈ _DFERMION && ( OD["part"]  = s.part  )
    "dfmix" ∈ _DFERMION && ( OD["dfmix"] = s.dfmix )
    "bsmix" ∈ _DFERMION && ( OD["bsmix"] = s.bsmix )
    #
    return OD
end

"""
    build_dfermion_dict()

Assemble the ordered dictionary, which is then converted into `dfa.in`,
for the `DFermion` code.

See also: [`DFERMION_PDFERMION`](@ref).
"""
function build_dfermion_dict()
    return struct_to_dict(PDFERMION)
end

#=
### *Customized Structs* : *iQIST Package*
=#

#=
*Remarks* :

There are three codes in the iQIST package, namely `ctseg`, `cthyb`, and
`atomic`. The official configuration files for the `ctseg` and `cthyb`
codes are `solver.ctqmc.in`, while the one for the `atomic` code is just
`solver.atomic.in`. These files are actually ini-like files.
=#

"""
    IQIST_PCTSEG

This struct encapsulates the parameters in the `solver.ctqmc.in` file. It
is for the `iQIST/ctseg` code only.

See also: [`PCTSEG`](@ref).
"""
mutable struct IQIST_PCTSEG
    isscf  :: I64 # Cycle
    isscr  :: I64 # Model
    isbnd  :: I64 # Symmetry
    isspn  :: I64 # Symmetry
    iswor  :: I64 # Measurement
    isort  :: I64 # Representation
    isobs  :: I64 # Measurement
    issus  :: I64 # Measurement
    isvrt  :: I64 # Measurement
    nband  :: I64 # Model
    nspin  :: I64 # Model
    norbs  :: I64 # Model
    ncfgs  :: I64 # Model
    niter  :: I64 # Cycle
    lemax  :: I64 # Representation
    legrd  :: I64 # Representation
    svmax  :: I64 # Representation
    svgrd  :: I64 # Representation
    mkink  :: I64 # Monte Carlo
    mfreq  :: I64 # Dimension
    nffrq  :: I64 # Dimension
    nbfrq  :: I64 # Dimension
    nfreq  :: I64 # Dimension
    ntime  :: I64 # Dimension
    nflip  :: I64 # Monte Carlo
    ntherm :: I64 # Monte Carlo
    nsweep :: I64 # Monte Carlo
    nwrite :: I64 # Monte Carlo
    nclean :: I64 # Monte Carlo
    nmonte :: I64 # Monte Carlo
    ncarlo :: I64 # Monte Carlo
    Uc     :: F64 # Model
    Jz     :: F64 # Model
    lc     :: F64 # Model
    wc     :: F64 # Model
    mune   :: F64 # Model
    beta   :: F64 # Model
    part   :: F64 # Model
    alpha  :: F64 # Cycle
end

"""
    IQIST_PCTHYB

This struct encapsulates the parameters in the `solver.ctqmc.in` file. It
is for the `iQIST/cthyb` code only.

See also: [`PCTHYB`](@ref).
"""
mutable struct IQIST_PCTHYB
end

"""
    IQIST_PATOMIC

This struct encapsulates the parameters in the `solver.atomic.in` file. It
is for the `iQIST/atomic` code only.

See also: [`PATOMIC`](@ref).
"""
mutable struct IQIST_PATOMIC
    ibasis :: I64 # Natural eigenbasis
    ictqmc :: I64 # Algorithm
    icu    :: I64 # Interaction
    icf    :: I64 # Natural eigenbasis
    isoc   :: I64 # Natural eigenbasis
    nband  :: I64 # Model
    nspin  :: I64 # Model
    norbs  :: I64 # Model
    ncfgs  :: I64 # Model
    nmini  :: I64 # Algorithm
    nmaxi  :: I64 # Algorithm
    Uc     :: F64 # Interaction
    Uv     :: F64 # Interaction
    Jz     :: F64 # Interaction
    Js     :: F64 # Interaction
    Jp     :: F64 # Interaction
    Ud     :: F64 # Interaction
    Jh     :: F64 # Interaction
    mune   :: F64 # Natural eigenbasis
    lambda :: F64 # Natural eigenbasis
end

"""
    _CTSEG

This set records the names of modified parameters, which will be presented
in the `solver.ctqmc.in` file. It is for the `iQIST/ctseg` code only.

Note that not all the parameters in `IQIST_PCTSEG` will be presented in
the `solver.ctqmc.in` file.

See also: [`IQIST_PCTSEG`](@ref).
"""
_CTSEG = Set{String}()

"""
    _CTHYB

This set records the names of modified parameters, which will be presented
in the `solver.ctqmc.in` file. It is for the `iQIST/cthyb` code only.

Note that not all the parameters in `IQIST_PCTHYB` will be presented in
the `solver.ctqmc.in` file.

See also: [`IQIST_PCTHYB`](@ref).
"""
_CTHYB = Set{String}()

"""
    _ATOMIC

This set records the names of modified parameters, which will be presented
in the `solver.atomic.in` file. It is for the `iQIST/atomic` code only.

Note that not all the parameters in `IQIST_PATOMIC` will be presented in
the `solver.atomic.in` file.

See also: [`IQIST_PATOMIC`](@ref).
"""
_ATOMIC = Set{String}()

"""
    PCTSEG

An instance for the `IQIST_PCTSEG` struct.

See also: [`IQIST_PCTSEG`](@ref).
"""
PCTSEG = IQIST_PCTSEG(
    1,        # isscf
    1,        # isscr
    1,        # isbnd
    1,        # isspn
    1,        # iswor
    1,        # isort
    1,        # isobs
    1,        # issus
    1,        # isvrt
    1,        # nband
    2,        # nspin
    2,        # norbs
    4,        # ncfgs
    20,       # niter
    32,       # lemax
    20001,    # legrd
    32,       # svmax
    2001,     # svgrd
    1024,     # mkink
    8193,     # mfreq
    32,       # nffrq
    8,        # nbfrq
    128,      # nfreq
    1024,     # ntime
    20000,    # nflip
    200000,   # ntherm
    20000000, # nsweep
    2000000,  # nwrite
    100000,   # nclean
    10,       # nmonte
    10,       # ncarlo
    4.0,      # Uc
    0.0,      # Jz
    1.0,      # lc
    1.0,      # wc
    2.0,      # mune
    8.0,      # beta
    0.5,      # part
    0.7       # alpha
)

"""
    PCTHYB

An instance for the `IQIST_PCTHYB` struct.

See also: [`IQIST_PCTHYB`](@ref).
"""
PCTHYB = IQIST_PCTHYB()

"""
    PATOMIC

An instance for the `IQIST_PATOMIC` struct.

See also: [`IQIST_PATOMIC`](@ref).
"""
PATOMIC = IQIST_PATOMIC(
    1,   # ibasis
    1,   # ictqmc
    1,   # icu
    0,   # icf
    0,   # isoc
    1,   # nband
    2,   # nspin
    2,   # norbs
    4,   # ncfgs
    0,   # nmini
    2,   # nmaxi
    2.0, # Uc
    2.0, # Uv
    0.0, # Jz
    0.0, # Js
    0.0, # Jp
    2.0, # Ud
    0.0, # Jh
    0.0, # mune
    0.0  # lambda
)

"""
    struct_to_dict(s::IQIST_PCTSEG)

Convert a struct to an ordered dictionary (for `IQIST_PCTSEG`).

See [`IQIST_PCTSEG`](@ref).
"""
function struct_to_dict(s::IQIST_PCTSEG)
    OD = OrderedDict{String,Any}()
    #
    "isscf"  ∈ _CTSEG && ( OD["isscf"]  = s.isscf  )
    "isscr"  ∈ _CTSEG && ( OD["isscr"]  = s.isscr  )
    "isbnd"  ∈ _CTSEG && ( OD["isbnd"]  = s.isbnd  )
    "isspn"  ∈ _CTSEG && ( OD["isspn"]  = s.isspn  )
    "iswor"  ∈ _CTSEG && ( OD["iswor"]  = s.iswor  )
    "isort"  ∈ _CTSEG && ( OD["isort"]  = s.isort  )
    "isobs"  ∈ _CTSEG && ( OD["isobs"]  = s.isobs  )
    "issus"  ∈ _CTSEG && ( OD["issus"]  = s.issus  )
    "isvrt"  ∈ _CTSEG && ( OD["isvrt"]  = s.isvrt  )
    "nband"  ∈ _CTSEG && ( OD["nband"]  = s.nband  )
    "nspin"  ∈ _CTSEG && ( OD["nspin"]  = s.nspin  )
    "norbs"  ∈ _CTSEG && ( OD["norbs"]  = s.norbs  )
    "ncfgs"  ∈ _CTSEG && ( OD["ncfgs"]  = s.ncfgs  )
    "niter"  ∈ _CTSEG && ( OD["niter"]  = s.niter  )
    "lemax"  ∈ _CTSEG && ( OD["lemax"]  = s.lemax  )
    "legrd"  ∈ _CTSEG && ( OD["legrd"]  = s.legrd  )
    "svmax"  ∈ _CTSEG && ( OD["svmax"]  = s.svmax  )
    "svgrd"  ∈ _CTSEG && ( OD["svgrd"]  = s.svgrd  )
    "mkink"  ∈ _CTSEG && ( OD["mkink"]  = s.mkink  )
    "mfreq"  ∈ _CTSEG && ( OD["mfreq"]  = s.mfreq  )
    "nffrq"  ∈ _CTSEG && ( OD["nffrq"]  = s.nffrq  )
    "nbfrq"  ∈ _CTSEG && ( OD["nbfrq"]  = s.nbfrq  )
    "nfreq"  ∈ _CTSEG && ( OD["nfreq"]  = s.nfreq  )
    "ntime"  ∈ _CTSEG && ( OD["ntime"]  = s.ntime  )
    "nflip"  ∈ _CTSEG && ( OD["nflip"]  = s.nflip  )
    "ntherm" ∈ _CTSEG && ( OD["ntherm"] = s.ntherm )
    "nsweep" ∈ _CTSEG && ( OD["nsweep"] = s.nsweep )
    "nwrite" ∈ _CTSEG && ( OD["nwrite"] = s.nwrite )
    "nclean" ∈ _CTSEG && ( OD["nclean"] = s.nclean )
    "nmonte" ∈ _CTSEG && ( OD["nmonte"] = s.nmonte )
    "ncarlo" ∈ _CTSEG && ( OD["ncarlo"] = s.ncarlo )
    "Uc"     ∈ _CTSEG && ( OD["Uc"]     = s.Uc     )
    "Jz"     ∈ _CTSEG && ( OD["Jz"]     = s.Jz     )
    "lc"     ∈ _CTSEG && ( OD["lc"]     = s.lc     )
    "wc"     ∈ _CTSEG && ( OD["wc"]     = s.wc     )
    "mune"   ∈ _CTSEG && ( OD["mune"]   = s.mune   )
    "beta"   ∈ _CTSEG && ( OD["beta"]   = s.beta   )
    "part"   ∈ _CTSEG && ( OD["part"]   = s.part   )
    "alpha"  ∈ _CTSEG && ( OD["alpha"]  = s.alpha  )
    #
    return OD
end

"""
    struct_to_dict(s::IQIST_PCTHYB)

Convert a struct to an ordered dictionary (for `IQIST_PCTHYB`).

See [`IQIST_PCTHYB`](@ref).
"""
function struct_to_dict(s::IQIST_PCTHYB)
    return OrderedDict{String,Any}(
        "key" => "value",
    )
end

"""
    struct_to_dict(s::IQIST_PATOMIC)

Convert a struct to an ordered dictionary (for `IQIST_PATOMIC`).

See [`IQIST_PATOMIC`](@ref).
"""
function struct_to_dict(s::IQIST_PATOMIC)
    OD = OrderedDict{String,Any}()
    #
    "ibasis" ∈ _ATOMIC && ( OD["ibasis"] = s.ibasis )
    "ictqmc" ∈ _ATOMIC && ( OD["ictqmc"] = s.ictqmc )
    "icu"    ∈ _ATOMIC && ( OD["icu"]    = s.icu    )
    "icf"    ∈ _ATOMIC && ( OD["icf"]    = s.icf    )
    "isoc"   ∈ _ATOMIC && ( OD["isoc"]   = s.isoc   )
    "nband"  ∈ _ATOMIC && ( OD["nband"]  = s.nband  )
    "nspin"  ∈ _ATOMIC && ( OD["nspin"]  = s.nspin  )
    "norbs"  ∈ _ATOMIC && ( OD["norbs"]  = s.norbs  )
    "ncfgs"  ∈ _ATOMIC && ( OD["ncfgs"]  = s.ncfgs  )
    "nmini"  ∈ _ATOMIC && ( OD["nmini"]  = s.nmini  )
    "nmaxi"  ∈ _ATOMIC && ( OD["nmaxi"]  = s.nmaxi  )
    "Uc"     ∈ _ATOMIC && ( OD["Uc"]     = s.Uc     )
    "Uv"     ∈ _ATOMIC && ( OD["Uv"]     = s.Uv     )
    "Jz"     ∈ _ATOMIC && ( OD["Jz"]     = s.Jz     )
    "Js"     ∈ _ATOMIC && ( OD["Js"]     = s.Js     )
    "Jp"     ∈ _ATOMIC && ( OD["Jp"]     = s.Jp     )
    "Ud"     ∈ _ATOMIC && ( OD["Ud"]     = s.Ud     )
    "Jh"     ∈ _ATOMIC && ( OD["Jh"]     = s.Jh     )
    "mune"   ∈ _ATOMIC && ( OD["mune"]   = s.mune   )
    "lambda" ∈ _ATOMIC && ( OD["lambda"] = s.lambda )
    #
    return OD
end

"""
    build_iqist_dict()

Assemble the ordered dictionary, which is then converted into
`solver.ctqmc.in` or `solver.atomic.in` file for the `ctseg`, `cthyb`,
and `atomic` codes.
"""
function build_iqist_dict(solver::String)
    @cswitch solver begin

        @case "ctseg"
            return struct_to_dict(PCTSEG)
            break

        @case "cthyb"
            return struct_to_dict(PCTHYB)
            break

        @case "atomic"
            return struct_to_dict(PATOMIC)
            break

        @default
            sorry()
            break

    end
end

#=
### *Customized Structs* : *ACFlow Toolkit*
=#

#=
*Remarks* :

The official configuration file for the ACFlow toolkit is `ac.toml`. It
includes `[BASE]`, `[MaxEnt]`, `[BarRat]`, `[NevanAC]`, `[StochAC]`,
`[StochSK]`, `[StochOM]`, and `[StochPX]` blocks.
=#

"""
    ACFLOW_PBASE

This struct represents the `[BASE]` block in the `ac.toml` file.
"""
mutable struct ACFLOW_PBASE
    finput  :: String
    solver  :: String
    ktype   :: String
    mtype   :: String
    grid    :: String
    mesh    :: String
    ngrid   :: I64
    nmesh   :: I64
    wmax    :: F64
    wmin    :: F64
    beta    :: F64
    offdiag :: Bool
    fwrite  :: Bool
end

"""
    ACFLOW_PMaxEnt

This struct represents the `[MaxEnt]` block in the `ac.toml` file.
"""
mutable struct ACFLOW_PMaxEnt
    method :: String
    stype  :: String
    nalph  :: I64
    alpha  :: F64
    ratio  :: F64
    blur   :: F64
end

"""
    ACFLOW_PBarRat

This struct represents the `[BarRat]` block in the `ac.toml` file.
"""
mutable struct ACFLOW_PBarRat
    atype   :: String
    denoise :: String
    epsilon :: F64
    pcut    :: F64
    eta     :: F64
end

"""
    ACFLOW_PNevanAC

This struct represents the `[NevanAC]` block in the `ac.toml` file.
"""
mutable struct ACFLOW_PNevanAC
    pick  :: Bool
    hardy :: Bool
    hmax  :: I64
    alpha :: F64
    eta   :: F64
end

"""
    ACFLOW_PStochAC

This struct represents the `[StochAC]` block in the `ac.toml` file.
"""
mutable struct ACFLOW_PStochAC
    nfine :: I64
    ngamm :: I64
    nwarm :: I64
    nstep :: I64
    ndump :: I64
    nalph :: I64
    alpha :: F64
    ratio :: F64
end

"""
    ACFLOW_PStochSK

This struct represents the `[StochSK]` block in the `ac.toml` file.
"""
mutable struct ACFLOW_PStochSK
    method :: String
    nfine  :: I64
    ngamm  :: I64
    nwarm  :: I64
    nstep  :: I64
    ndump  :: I64
    retry  :: I64
    theta  :: F64
    ratio  :: F64
end

"""
    ACFLOW_PStochOM

This struct represents the `[StochOM]` block in the `ac.toml` file.
"""
mutable struct ACFLOW_PStochOM
    ntry  :: I64
    nstep :: I64
    nbox  :: I64
    sbox  :: F64
    wbox  :: F64
    norm  :: F64
end

"""
    ACFLOW_PStochPX

This struct represents the `[StochPX]` block in the `ac.toml` file.
"""
mutable struct ACFLOW_PStochPX
    method :: String
    nfine  :: I64
    npole  :: I64
    ntry   :: I64
    nstep  :: I64
    theta  :: F64
    eta    :: F64
end

"""
    PBASE

An instance for the `ACFLOW_PBASE` struct.

See also: [`ACFLOW_PBASE`](@ref).
"""
PBASE = ACFLOW_PBASE(
    "giw.data", # finput
    "MaxEnt",   # solver
    "fermi",    # ktype
    "flat",     # mtype
    "ffreq",    # grid
    "linear",   # mesh
    10,         # ngrid
    501,        # nmesh
    5.0,        # wmax
    -5.0,       # wmin
    10.0,       # beta
    false,      # offdiag
    true        # fwrite
)

"""
    PMaxEnt

An instance for the `ACFLOW_PMaxEnt` struct.

See also: [`ACFLOW_PMaxEnt`](@ref).
"""
PMaxEnt = ACFLOW_PMaxEnt(
    "chi2kink", # method
    "sj",       # stype
    12,         # nalph
    1e9,        # alpha
    10.0,       # ratio
    -1.0        # blur
)

"""
    PBarRat

An instance for the `ACFLOW_PBarRat` struct.

See also: [`ACFLOW_PBarRat`](@ref).
"""
PBarRat = ACFLOW_PBarRat(
    "cont",  # atype
    "prony", # denoise
    1e-10,   # epsilon
    1e-3,    # pcut
    1e-2     # eta
)

"""
    PNevanAC

An instance for the `ACFLOW_PNevanAC` struct.

See also: [`ACFLOW_PNevanAC`](@ref).
"""
PNevanAC = ACFLOW_PNevanAC(
    true, # pick
    true, # hardy
    50,   # hmax
    1e-4, # alpha
    1e-2  # eta
)

"""
    PStochAC

An instance for the `ACFLOW_PStochAC` struct.

See also: [`ACFLOW_PStochAC`](@ref).
"""
PStochAC = ACFLOW_PStochAC(
    10000,   # nfine
    512,     # ngamm
    4000,    # nwarm
    4000000, # nstep
    40000,   # ndump
    20,      # nalph
    1.00,    # alpha
    1.20     # ratio
)

"""
    PStochSK

An instance for the `ACFLOW_PStochSK` struct.

See also: [`ACFLOW_PStochSK`](@ref).
"""
PStochSK = ACFLOW_PStochSK(
    "chi2min", # method
    100000,    # nfine
    1000,      # ngamm
    1000,      # nwarm
    20000,     # nstep
    200,       # ndump
    10,        # retry
    1e+6,      # theta
    0.90       # ratio
)

"""
    PStochOM

An instance for the `ACFLOW_PStochOM` struct.

See also: [`ACFLOW_PStochOM`](@ref).
"""
PStochOM = ACFLOW_PStochOM(
    2000,  # ntry
    1000,  # nstep
    100,   # nbox
    0.005, # sbox
    0.02,  # wbox
    -1.0   # norm
)

"""
    PStochPX

An instance for the `ACFLOW_PStochPX` struct.

See also: [`ACFLOW_PStochPX`](@ref).
"""
PStochPX = ACFLOW_PStochPX(
    "mean",  # method
    100000,  # nfine
    200,     # npole
    1000,    # ntry
    1000000, # nstep
    1e+6,    # theta
    1e-4     # eta
)

"""
    struct_to_dict(s::ACFLOW_PBASE)

Convert a struct to an ordered dictionary (for `ACFLOW_PBASE`).

See [`ACFLOW_PBASE`](@ref).
"""
function struct_to_dict(s::ACFLOW_PBASE)
    return OrderedDict{String,Any}(
        "finput"  => s.finput,
        "solver"  => s.solver,
        "ktype"   => s.ktype,
        "mtype"   => s.mtype,
        "grid"    => s.grid,
        "mesh"    => s.mesh,
        "ngrid"   => s.ngrid,
        "nmesh"   => s.nmesh,
        "wmax"    => s.wmax,
        "wmin"    => s.wmin,
        "beta"    => s.beta,
        "offdiag" => s.offdiag,
        "fwrite"  => s.fwrite,
    )
end

"""
    struct_to_dict(s::ACFLOW_PMaxEnt)

Convert a struct to an ordered dictionary (for `ACFLOW_PMaxEnt`).

See [`ACFLOW_PMaxEnt`](@ref).
"""
function struct_to_dict(s::ACFLOW_PMaxEnt)
    return OrderedDict{String,Any}(
        "method" => s.method,
        "stype"  => s.stype,
        "nalph"  => s.nalph,
        "alpha"  => s.alpha,
        "ratio"  => s.ratio,
        "blur"   => s.blur,
    )
end

"""
    struct_to_dict(s::ACFLOW_PBarRat)

Convert a struct to an ordered dictionary (for `ACFLOW_PBarRat`).

See [`ACFLOW_PBarRat`](@ref).
"""
function struct_to_dict(s::ACFLOW_PBarRat)
    return OrderedDict{String,Any}(
        "atype"   => s.atype,
        "denoise" => s.denoise,
        "epsilon" => s.epsilon,
        "pcut"    => s.pcut,
        "eta"     => s.eta,
    )
end

"""
    struct_to_dict(s::ACFLOW_PNevanAC)

Convert a struct to an ordered dictionary (for `ACFLOW_PNevanAC`).

See [`ACFLOW_PNevanAC`](@ref).
"""
function struct_to_dict(s::ACFLOW_PNevanAC)
    return OrderedDict{String,Any}(
        "pick"    => s.pick,
        "hardy"   => s.hardy,
        "hmax"    => s.hmax,
        "alpha"   => s.alpha,
        "eta"     => s.eta,
    )
end

"""
    struct_to_dict(s::ACFLOW_PStochAC)

Convert a struct to an ordered dictionary (for `ACFLOW_PStochAC`).

See [`ACFLOW_PStochAC`](@ref).
"""
function struct_to_dict(s::ACFLOW_PStochAC)
    return OrderedDict{String,Any}(
        "nfine"   => s.nfine,
        "ngamm"   => s.ngamm,
        "nwarm"   => s.nwarm,
        "nstep"   => s.nstep,
        "ndump"   => s.ndump,
        "nalph"   => s.nalph,
        "alpha"   => s.alpha,
        "ratio"   => s.ratio,
    )
end

"""
    struct_to_dict(s::ACFLOW_PStochSK)

Convert a struct to an ordered dictionary (for `ACFLOW_PStochSK`).

See [`ACFLOW_PStochSK`](@ref).
"""
function struct_to_dict(s::ACFLOW_PStochSK)
    return OrderedDict{String,Any}(
        "method"  => s.method,
        "nfine"   => s.nfine,
        "ngamm"   => s.ngamm,
        "nwarm"   => s.nwarm,
        "nstep"   => s.nstep,
        "ndump"   => s.ndump,
        "retry"   => s.retry,
        "theta"   => s.theta,
        "ratio"   => s.ratio,
    )
end

"""
    struct_to_dict(s::ACFLOW_PStochOM)

Convert a struct to an ordered dictionary (for `ACFLOW_PStochOM`).

See [`ACFLOW_PStochOM`](@ref).
"""
function struct_to_dict(s::ACFLOW_PStochOM)
    return OrderedDict{String,Any}(
        "ntry"    => s.ntry,
        "nstep"   => s.nstep,
        "nbox"    => s.nbox,
        "sbox"    => s.sbox,
        "wbox"    => s.wbox,
        "norm"    => s.norm,
    )
end

"""
    struct_to_dict(s::ACFLOW_PStochPX)

Convert a struct to an ordered dictionary (for `ACFLOW_PStochPX`).

See [`ACFLOW_PStochPX`](@ref).
"""
function struct_to_dict(s::ACFLOW_PStochPX)
    return OrderedDict{String,Any}(
        "method"  => s.method,
        "nfine"   => s.nfine,
        "npole"   => s.npole,
        "ntry"    => s.ntry,
        "nstep"   => s.nstep,
        "theta"   => s.theta,
        "eta"     => s.eta,
    )
end

"""
    build_acflow_dict()

Assemble the ordered dictionary, which is then converted into `ac.toml`,
for the `ACFlow` toolkit.
"""
function build_acflow_dict()
    @cswitch PBASE.solver begin

        @case "MaxEnt"
            return OrderedDict{String,Any}(
                "BASE" => struct_to_dict(PBASE),
                "MaxEnt" => struct_to_dict(PMaxEnt)
            )
            break

        @case "BarRat"
            return OrderedDict{String,Any}(
                "BASE" => struct_to_dict(PBASE),
                "BarRat" => struct_to_dict(PBarRat)
            )
            break

        @case "NevanAC"
            return OrderedDict{String,Any}(
                "BASE" => struct_to_dict(PBASE),
                "NevanAC" => struct_to_dict(PNevanAC)
            )
            break

        @case "StochAC"
            return OrderedDict{String,Any}(
                "BASE" => struct_to_dict(PBASE),
                "StochAC" => struct_to_dict(PStochAC)
            )
            break

        @case "StochSK"
            return OrderedDict{String,Any}(
                "BASE" => struct_to_dict(PBASE),
                "StochSK" => struct_to_dict(PStochSK)
            )
            break

        @case "StochOM"
            return OrderedDict{String,Any}(
                "BASE" => struct_to_dict(PBASE),
                "StochOM" => struct_to_dict(PStochOM)
            )
            break

        @case "StochPX"
            return OrderedDict{String,Any}(
                "BASE" => struct_to_dict(PBASE),
                "StochPX" => struct_to_dict(PStochPX)
            )
            break

        @default
            sorry()
            break

    end
end

#=
### *Customized Structs* : *ACTest Toolkit*
=#

#=
*Remarks* :

The official configuration file for the ACTest toolkit is `act.toml`. It
includes two blocks: `[Test]` and `[Solver]`.
=#

"""
    ACTEST_PTEST

This struct represents the `[Test]` block in the `act.toml` file.
"""
mutable struct ACTEST_PTEST
    solver :: String
    ptype  :: String
    ktype  :: String
    grid   :: String
    mesh   :: String
    ngrid  :: I64
    nmesh  :: I64
    ntest  :: I64
    nbins  :: I64
    wmax   :: F64
    wmin   :: F64
    pmax   :: F64
    pmin   :: F64
    beta   :: F64
    noise  :: F64
    lcorr  :: F64
    tcorr  :: Bool
    fnpd   :: Bool
    fpbc   :: Bool
    lpeak  :: Array
end

"""
    PTEST

An instance for the `ACTEST_PTEST` struct.

See also: [`ACTEST_PTEST`](@ref).
"""
PTEST = ACTEST_PTEST(
    "MaxEnt", # solver
    "gauss",  # ptype
    "fermi",  # ktype
    "ffreq",  # grid
    "linear", # mesh
    10,       # ngrid
    501,      # nmesh
    100,      # ntest
    1,        # nbins
    5.0,      # wmax
    -5.0,     # wmin
    4.0,      # pmax
    -4.0,     # pmin
    10.0,     # beta
    1.0e-6,   # noise
    0.5,      # lcorr
    false,    # tcorr
    false,    # fnpd
    false,    # fpbc
    [1,2,3]   # lpeak
)

"""
    struct_to_dict(s::ACTEST_PTEST)

Convert a struct to an ordered dictionary (for `ACTEST_PTEST`).

See [`ACTEST_PTEST`](@ref).
"""
function struct_to_dict(s::ACTEST_PTEST)
    return OrderedDict{String,Any}(
        "solver" => s.solver,
        "ptype"  => s.ptype,
        "ktype"  => s.ktype,
        "grid"   => s.grid,
        "mesh"   => s.mesh,
        "ngrid"  => s.ngrid,
        "nmesh"  => s.nmesh,
        "ntest"  => s.ntest,
        "nbins"  => s.nbins,
        "wmax"   => s.wmax,
        "wmin"   => s.wmin,
        "pmax"   => s.pmax,
        "pmin"   => s.pmin,
        "beta"   => s.beta,
        "noise"  => s.noise,
        "lcorr"  => s.lcorr,
        "tcorr"  => s.tcorr,
        "fnpd"   => s.fnpd,
        "fpbc"   => s.fpbc,
        "lpeak"  => s.lpeak,
    )
end

"""
    build_actest_dict()

Assemble the ordered dictionary, which is then converted into `act.toml`,
for the `ACTest` toolkit.
"""
function build_actest_dict()
    @cswitch PTEST.solver begin

        @case "MaxEnt"
            return OrderedDict{String,Any}(
                "Test" => struct_to_dict(PTEST),
                "Solver" => struct_to_dict(PMaxEnt)
            )
            break

        @case "BarRat"
            return OrderedDict{String,Any}(
                "Test" => struct_to_dict(PTEST),
                "Solver" => struct_to_dict(PBarRat)
            )
            break

        @case "NevanAC"
            return OrderedDict{String,Any}(
                "Test" => struct_to_dict(PTEST),
                "Solver" => struct_to_dict(PNevanAC)
            )
            break

        @case "StochAC"
            return OrderedDict{String,Any}(
                "Test" => struct_to_dict(PTEST),
                "Solver" => struct_to_dict(PStochAC)
            )
            break

        @case "StochSK"
            return OrderedDict{String,Any}(
                "Test" => struct_to_dict(PTEST),
                "Solver" => struct_to_dict(PStochSK)
            )
            break

        @case "StochOM"
            return OrderedDict{String,Any}(
                "Test" => struct_to_dict(PTEST),
                "Solver" => struct_to_dict(PStochOM)
            )
            break

        @case "StochPX"
            return OrderedDict{String,Any}(
                "Test" => struct_to_dict(PTEST),
                "Solver" => struct_to_dict(PStochPX)
            )
            break

        @default
            sorry()
            break

    end
end
