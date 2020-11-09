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
    adaptor = dict["adaptor"]
    impurity = dict["impurity"]
    @show case
    @show dft
    @show dmft
    @show solver
    @show adaptor
    @show impurity
end
