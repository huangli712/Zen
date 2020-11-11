using TOML

"""
    parse_config(f::AbstractString)

Parse the configuration file (toml format)
"""
function parse_config(f::AbstractString)
    dict = TOML.parsefile(f)

    case = dict["case"]
    dft = dict["dft"]
    dmft = dict["dmft"]
    dft_dmft = dict["dft_dmft"]

    return case, dft, dmft, dft_dmft
end
