include("zen.jl")
using .Zen

welcome()

check_version()

toml = check_toml()

case, dft, dmft, solver, impurity, dft_dmft = parse_config(toml)

make_trees(impurity)

it = IterInfo(0, 0, 0)

dft_driver(it, dft)
@show it.dft_iter

goodbye()
