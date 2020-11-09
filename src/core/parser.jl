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
    solver = dict["solver"]
    impurity = dict["impurity"]
    dft_dmft = dict["dft_dmft"]

    return case, dft, dmft, solver, impurity, dft_dmft
end
