using TOML

"""
    parse_config(f::AbstractString, key::AbstractString)

Parse the configuration file (toml format)
"""
function parse_config(f::AbstractString, key::AbstractString)
    dict = TOML.parsefile(f)

    if haskey(dict, key)
        dict[key]
    else
        error("Do not have this key: $key")
    end
end

function parse_mpi(key::AbstractString)
    f = "MPI.toml"
    if isfile(f)
        dict = TOML.parsefile("MPI.toml")

        if haskey(dict, key)
            dict[key]
        else
            error("Do not have this key: $key")
        end
    else
        ""
    end
end
