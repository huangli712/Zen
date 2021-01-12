#
# project : pansy
# source  : config.jl
# author  : Li Huang (lihuang.dmft@gmail.com)
# status  : unstable
# comment :
#
# last modified: 2021/01/12
#

"""
    parse_toml(f::String, key::String, necessary::Bool)

Parse the configuration file (toml format)
"""
function parse_toml(f::String, key::String, necessary::Bool)
    if isfile(f)
        dict = TOML.parsefile(f)

        if haskey(dict, key)
            dict[key]
        else
            error("Do not have this key: $key in file: $f")
        end
    else
        if necessary
            error("Please make sure that the file $f really exists")
        else
            nothing
        end
    end
end

"""
    parse_toml(f::String, necessary::Bool)

Parse the configuration file (toml format)
"""
function parse_toml(f::String, necessary::Bool)
    if isfile(f)
        dict = TOML.parsefile(f)
    else
        if necessary
            error("Please make sure that the file $f really exists")
        else
            nothing
        end
    end
end

"""
    renew_config(cfg::Dict{String,Any})

Copy configurations from cfg to PCASE, PDFT, PDMFT, PIMP, and PSOLVER
"""
function renew_config(cfg::Dict{String,Any})
    # for case block
    PCASE["case"][1] = cfg["case"]

    # for dft block
    dft = cfg["dft"]
    for key in keys(dft)
        if haskey(PDFT, key)
            PDFT[key][1] = dft[key]
        else
            error("Sorry, $key is not supported currently")
        end
    end

    # for dmft block
    dmft = cfg["dmft"]
    for key in keys(dmft)
        if haskey(PDMFT, key)
            PDMFT[key][1] = dmft[key]
        else
            error("Sorry, $key is not supported currently")
        end
    end

    # for impurity block
    impurity = cfg["impurity"]
    for key in keys(impurity)
        if haskey(PIMP, key)
            PIMP[key][1] = impurity[key]
        else
            error("Sorry, $key is not supported currently")
        end
    end

    # for solver block
    solver = cfg["solver"]
    for key in keys(solver)
        if haskey(PSOLVER, key)
            PSOLVER[key][1] = solver[key]
        else
            error("Sorry, $key is not supported currently")
        end
    end
end

"""
    check_config()

Validate the correctness and consistency of configurations
"""
function check_config()
    # 1. check types and existences
    #
    # check case block
    for key in keys(PCASE)
        val = PCASE[key]
        _v(val)
    end
    #
    # check dft block
    for key in keys(PDFT)
        val = PDFT[key]
        _v(val)
    end
    #
    # check dmft block
    for key in keys(PDMFT)
        val = PDMFT[key]
        _v(val)
    end
    #
    # check impurity block
    for key in keys(PIMP)
        val = PIMP[key]
        _v(val)
    end
    #
    # check solver block
    for key in keys(PSOLVER)
        val = PSOLVER[key]
        _v(val)
    end

    # 2. check rationalities
    #
    # check dft block
    @assert _d("engine") in ("vasp", "wannier")
    @assert _d("smear") in ("m-p", "gauss", "tetra", missing)
    @assert _d("kmesh") in ("accurate", "medium", "coarse", "file", missing)
    #
    # check dmft block
    @assert _m("mode") in (1, 2)
    @assert _m("axis") in (1, 2)
    @assert _m("niter") > 0
    @assert _m("dcount") in ("fll1", "fll2", "amf")
    @assert _m("beta") >= 0.0
    #
    # check solver block
    @assert _s("engine") in ("ct_hub1", "ct_hub2", "hub1", "norg")
    #
    # please add more assertion statements here

    # 3. check self-consistency
    #
    # check dft block
    if _d("lspinorb")
        @assert _d("lspins")
    end
    if _d("lproj")
        @assert !_d("lsymm") && !isa(_d("sproj"), Missing)
    end
    #
    # check solver block
    if _s("engine") in ("ct_hub1", "ct_hub2", "hub1")
        @assert _m("axis") === 1 # imaginary axis
    elseif _s("engine") in ("norg")
        @assert _m("axis") === 2 # real axis
    end
    #
    # please add more assertion statements here
end

"""
    _v(val::Array{Any,1})

Verify the value array
"""
@inline function _v(val::Array{Any,1})
    # to check if the value is updated
    if isa(val[1], Missing) && val[2] > 0
        error("Sorry, key shoule be set")
    end

    # to check if the type of value is correct
    if !isa(val[1], Missing) && !isa(val[1], eval(val[3]))
        error("Sorry, type of key is wrong")
    end
end

"""
    _c(key::String)

Extract configurations from dict: PCASE
"""
@inline function _c(key::String)
    if haskey(PCASE, key)
        PCASE[key][1]
    else
        error("Sorry, PCASE does not contain key: $key")
    end
end

"""
    _d(key::String)

Extract configurations from dict: PDFT
"""
@inline function _d(key::String)
    if haskey(PDFT, key)
        PDFT[key][1]
    else
        error("Sorry, PDFT does not contain key: $key")
    end
end

"""
    str_d(key::String)

Extract configurations from dict: PDFT, convert them into strings
"""
@inline function str_d(key::String)
    if haskey(PDFT, key)
        if PDFT[key][3] === :Array 
            join(PDFT[key][1], "; ")
        else
            string(PDFT[key][1])
        end
    else
        error("Sorry, PDFT does not contain key: $key")
    end
end

"""
    _m(key::String)

Extract configurations from dict: PDMFT
"""
@inline function _m(key::String)
    if haskey(PDMFT, key)
        PDMFT[key][1]
    else
        error("Sorry, PDMFT does not contain key: $key")
    end
end

"""
    str_m(key::String)

Extract configurations from dict: PDMFT, convert them into strings
"""
@inline function str_m(key::String)
    if haskey(PDMFT, key)
        if PDMFT[key][3] === :Array
            join(PDMFT[key][1], "; ")
        else
            string(PDMFT[key][1])
        end
    else
        error("Sorry, PDMFT does not contain key: $key")
    end
end

"""
    _i(key::String)

Extract configurations from dict: PIMP
"""
@inline function _i(key::String)
    if haskey(PIMP, key)
        PIMP[key][1]
    else
        error("Sorry, PIMP does not contain key: $key")
    end
end

"""
    str_i(key::String)

Extract configurations from dict: PIMP, convert them into strings
"""
@inline function str_i(key::String)
    if haskey(PIMP, key)
        if PIMP[key][3] === :Array
            join(PIMP[key][1], "; ")
        else
            string(PIMP[key][1])
        end
    else
        error("Sorry, PIMP does not contain key: $key")
    end
end

"""
    _s(key::String)

Extract configurations from dict: PSOLVER
"""
@inline function _s(key::String)
    if haskey(PSOLVER, key)
        PSOLVER[key][1]
    else
        error("Sorry, PSOLVER does not contain key: $key")
    end
end

@inline function str_s(key::String)
    if haskey(PSOLVER, key)
        if PSOLVER[key][3] === :Array
            join(PSOLVER[key][1], "; ")
        else
            string(PSOLVER[key][1])
        end
    else
        error("Sorry, PSOLVER does not contain key: $key")
    end
end
