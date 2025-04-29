#
# Project : Lily
# Source  : config.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2025/04/28
#

"""
    inp_toml(f::String, necessary::Bool)

Parse the configuration file (in toml format). It reads the whole file.

### Arguments
* f -> Filename of configuration.
* necessary -> If it is true and configuration is absent, raise an error.

### Returns
* dict -> A dictionary that contains all the key-value pairs.
"""
function inp_toml(f::String, necessary::Bool)
    if isfile(f)
        dict = TOML.parsefile(f)
        return dict
    else
        if necessary
            error("Please make sure that the file $f really exists")
        else
            nothing
        end
    end
end

"""
    fil_dict(cfg::Dict{String,Any})

Transfer configurations from dict `cfg` to internal dict (`PTEST`). In
other words, all the relevant internal dicts should be filled / updated
in this function.

### Arguments
* cfg -> A dict that contains all the configurations (from act.toml).

### Returns
N/A
"""
function fil_dict(cfg::Dict{String,Any})
    TEST = cfg["Test"]
    for key in keys(TEST)
        if haskey(PTEST, key)
            PTEST[key][1] = TEST[key]
        else
            error("Sorry, $key is not supported currently")
        end
    end
end

"""
    see_dict()

Display all of the relevant configuration parameters to the terminal.

### Arguments
N/A

### Returns
N/A

See also: [`fil_dict`](@ref).
"""
function see_dict()
    println("[ Test ]")
    #
    println("solver  : ", get_t("solver") )
    println("ptype   : ", get_t("ptype")  )
    println("ktype   : ", get_t("ktype")  )
    println("grid    : ", get_t("grid")   )
    println("mesh    : ", get_t("mesh")   )
    println("ngrid   : ", get_t("ngrid")  )
    println("nmesh   : ", get_t("nmesh")  )
    println("ntest   : ", get_t("ntest")  )
    println("wmax    : ", get_t("wmax")   )
    println("wmin    : ", get_t("wmin")   )
    println("pmax    : ", get_t("pmax")   )
    println("pmin    : ", get_t("pmin")   )
    println("beta    : ", get_t("beta")   )
    println("noise   : ", get_t("noise")  )
    println("lcorr   : ", get_t("lcorr")  )
    println("tcorr   : ", get_t("tcorr")  )
    println("offdiag : ", get_t("offdiag"))
    println("lpeak   : ", get_t("lpeak")  )
    println("pmesh   : ", get_t("pmesh")  )
    #
    println()
end

"""
    rev_dict(TEST::Dict{String,Any})

Setup the configuration dictionary: `PTEST`.

### Arguments
* TEST -> A dict that contains configurations from the [Test] block.

### Returns
N/A

See also: [`PTEST`](@ref).
"""
function rev_dict(TEST::Dict{String,Any})
    for key in keys(TEST)
        if haskey(PTEST, key)
            PTEST[key][1] = TEST[key]
        else
            error("Sorry, $key is not supported currently")
        end
    end
    foreach(x -> _v(x.first, x.second), PTEST)
end

"""
    rev_dict(TEST::Dict{String,Vector{Any}})

Setup the configuration dictionary: `PTEST`.

### Arguments
* TEST -> A dict that contains configurations from the [Test] block.

### Returns
N/A

See also: [`PTEST`](@ref).
"""
function rev_dict(TEST::Dict{String,Vector{Any}})
    for key in keys(TEST)
        if haskey(PTEST, key)
            PTEST[key][1] = TEST[key][1]
        else
            error("Sorry, $key is not supported currently")
        end
    end
    foreach(x -> _v(x.first, x.second), PTEST)
end

"""
    chk_dict()

Validate the correctness and consistency of configurations.

### Arguments
N/A

### Returns
N/A

See also: [`fil_dict`](@ref), [`_v`](@ref).
"""
function chk_dict()
    #
    # The MiniPole solver implements the minimal pole method (MPM), which
    # is provided by the MiniPole code. Please see the following pages
    # and references for more details:
    #
    # https://github.com/Green-Phys/MiniPole
    # Phys. Rev. B 110, 235131 (2024)
    # Phys. Rev. B 110, 035154 (2024)
    #
    # The other solvers are provided by the ACFlow toolkit.
    #
    @assert get_t("solver") in (
        "MaxEnt",
        "BarRat", "NevanAC",
        "StochAC", "StochSK", "StochOM", "StochPX",
        "MiniPole"
    )
    @assert get_t("ptype") in (
        "gauss", "lorentz", "delta", "rectangle", "risedecay",
        "random1", "random2"
    )
    @assert get_t("ktype") in ("fermi", "boson", "bsymm")
    @assert get_t("grid") in ("ftime", "btime", "ffreq", "bfreq")
    @assert get_t("mesh") in ("linear", "tangent", "lorentz", "halflorentz")
    @assert get_t("ngrid") ≥ 1
    @assert get_t("nmesh") ≥ 1
    @assert get_t("ntest") ≥ 1
    @assert get_t("wmax") > get_t("wmin")
    @assert get_t("pmax") > get_t("pmin")
    @assert get_t("beta") ≥ 0.0
    @assert get_t("noise") ≥ 0.0
    @assert get_t("lcorr") > 0.0
    @assert length(get_t("lpeak")) ≥ 1
    @assert all(x -> x > 0, get_t("lpeak"))

    foreach(x -> _v(x.first, x.second), PTEST)
end

"""
    _v(key::String, val::Array{Any,1})

Verify the value array. Called by chk_dict() function only.

### Arguments
* key -> Key of parameter.
* val -> Value of parameter.

### Returns
N/A

See also: [`chk_dict`](@ref).
"""
@inline function _v(key::String, val::Array{Any,1})
    # To check if the value is updated
    if isa(val[1], Missing) && val[2] > 0
        error("Sorry, key ($key) shoule be set")
    end

    # To check if the type of value is correct
    if !isa(val[1], Missing) && !isa(val[1], eval(val[3]))
        error("Sorry, type of key ($key) is wrong")
    end
end

"""
    get_t(key::String)

Extract configurations from dict: PTEST.

### Arguments
* key -> Key of parameter.

### Returns
* value -> Value of parameter.

See also: [`PTEST`](@ref).
"""
@inline function get_t(key::String)
    if haskey(PTEST, key)
        PTEST[key][1]
    else
        error("Sorry, PTEST does not contain key: $key")
    end
end
