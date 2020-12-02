"""
    parse_toml(f::AbstractString, key::AbstractString, necessary::Bool)

Parse the configuration file (toml format)
"""
function parse_toml(f::AbstractString, key::AbstractString, necessary::Bool)
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
    parse_toml(f::AbstractString, necessary::Bool)

Parse the configuration file (toml format)
"""
function parse_toml(f::AbstractString, necessary::Bool)
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
    renew_params(cfg::Dict{String,Any})

Copy parameters from cfg to PCASE, PDFT, PDMFT, PIMP, and PSOLVER
"""
function renew_params(cfg::Dict{String,Any})
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
    check_params()

Validate the correctness and consistency of parameters
"""
function check_params()
    # check case block
    for key in keys(PCASE)
        if isa(PCASE[key][1], Missing) && PCASE[key][2] > 0
            error("Sorry, $key shoule be set")
        end

        if !isa(PCASE[key][1], PCASE[key][3])
            error("Sorry, type of $key is wrong")
        end
    end

    # check dft block
    for key in keys(PDFT)
        if isa(PDFT[key][1], Missing) && PDFT[key][2] > 0
            error("Sorry, $key shoule be set")
        end

        if !isa(PDFT[key][1], PDFT[key][3])
            error("Sorry, type of $key is wrong")
        end
    end

    # check dmft block
    for key in keys(PDMFT)
        if isa(PDMFT[key][1], Missing) && PDMFT[key][2] > 0
            error("Sorry, $key shoule be set")
        end

        if !isa(PDMFT[key][1], PDMFT[key][3])
            error("Sorry, type of $key is wrong")
        end
    end

    # check impurity block
    for key in keys(PIMP)
        if isa(PIMP[key][1], Missing) && PIMP[key][2] > 0
            error("Sorry, $key shoule be set")
        end

        if !isa(PIMP[key][1], PIMP[key][3])
            error("Sorry, type of $key is wrong")
        end
    end

    # check solver block
    for key in keys(PSOLVER)
        if isa(PSOLVER[key][1], Missing) && PSOLVER[key][2] > 0
            error("Sorry, $key shoule be set")
        end

        if !isa(PSOLVER[key][1], PSOLVER[key][3])
            error("Sorry, type of $key is wrong")
        end
    end
end

"""
    Param(dict::Dict{String,Any}, key::String)

Extract parameters from the given dictionary
"""
function Param(dict::Dict{String,Any}, key::String)
    dict[key][1]
end

@inline function _c(key::String)
    if haskey(PCASE, key)
        PCASE[key][1]
    else
        error("Sorry, PCASE does not contain key: $key")
    end
end

@inline function _d(key::String)
    if haskey(PDFT, key)
        PDFT[key][1]
    else
        error("Sorry, PDFT does not contain key: $key")
    end
end

@inline function _m(key::String)
    if haskey(PDMFT, key)
        PDMFT[key][1]
    else
        error("Sorry, PDMFT does not contain key: $key")
    end
end

@inline function _i(key::String)
    if haskey(PIMP, key)
        PIMP[key][1]
    else
        error("Sorry, PIMP does not contain key: $key")
    end
end

@inline function _s(key::String)
    if haskey(PSOLVER, key)
        PSOLVER[key][1]
    else
        error("Sorry, PSOLVER does not contain key: $key")
    end
end
