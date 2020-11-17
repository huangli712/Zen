using TOML

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

function parse_poscar()
end

function parse_projcar()
end

function parse_locproj()
end

function parse_ibzkpt()
end

function parse_eigenval()
end
