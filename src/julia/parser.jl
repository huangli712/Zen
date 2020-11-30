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
    parse_dict(cfg::Dict{String,Any})

Copy parameters from cfg to PCASE, PDFT, PDMFT, PIMP, and PSOLVER
"""
function parse_dict(cfg::Dict{String,Any})
    case = cfg["case"]
    PCASE["case"][1] = case

    dft = cfg["dft"]
    for key in keys(dft)
        if haskey(PDFT, key)
            PDFT[key][1] = dft[key]
        else
            error("Sorry, $key is not supported currently")
        end
    end

    dmft = cfg["dmft"]
    for key in keys(dmft)
        if haskey(PDMFT, key)
            PDMFT[key][1] = dmft[key]
        else
            error("Sorry, $key is not supported currently")
        end
    end

    impurity = cfg["impurity"]
    for key in keys(impurity)
        if haskey(PIMP, key)
            PIMP[key][1] = impurity[key]
        else
            error("Sorry, $key is not supported currently")
        end
    end

    solver = cfg["solver"]
    for key in keys(solver)
        if haskey(PSOLVER, key)
            PSOLVER[key][1] = solver[key]
        else
            error("Sorry, $key is not supported currently")
        end
    end
end

function validate_params()
    for key in keys(PCASE)
        if isa(PCASE[key][1], Missing)
            error("Sorry, $key shoule be set")
        end
    end
end
