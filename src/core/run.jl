include("zen.jl")
using .Zen

welcome()

check_version()

toml = check_toml()

case, dft, dmft, solver, impurity, dft_dmft = parse_config(toml)

@show case
@show dft
@show dmft
@show solver
@show impurity
@show dft_dmft

make_trees(impurity)

dft_driver(dft)

goodbye()
