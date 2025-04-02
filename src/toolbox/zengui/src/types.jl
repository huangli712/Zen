#
# Project : Camellia
# Source  : types.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2025/04/01
#

mutable struct MenuFlags
    _OPEN     :: Bool
    _SAVE     :: Bool
    _EXIT     :: Bool
    #
    ZEN       :: Bool
    DYSON     :: Bool
    DFERMION  :: Bool
    CTSEG     :: Bool
    CTHYB     :: Bool
    ATOMIC    :: Bool
    ACFLOW    :: Bool
    ACTEST    :: Bool
    #
    _CLASSIC  :: Bool
    _DARK     :: Bool
    _LIGHT    :: Bool
    #
    _ZEN      :: Bool
    _DYSON    :: Bool
    _DFERMION :: Bool
    _IQIST    :: Bool
    _ACFLOW   :: Bool
    _ACTEST   :: Bool
    _ZENGUI   :: Bool
    _ABOUT    :: Bool
end

FMENU = MenuFlags(
    false,
    false,
    false,
    false,
    false,
    false,
    false,
    false,
    false,
    false,
    false,
    false,
    false,
    false,
    false,
    false,
    false,
    false,
    false,
    false,
    false,
    false
)

mutable struct ACFLOW_PBASE
    finput :: String
    solver :: String
end

PBASE = ACFLOW_PBASE(
    "giw.data",
    "MaxEnt"
)
