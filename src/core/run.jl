include("zen.jl")
using .Zen

welcome()

check_version()

check_toml()
exit(-1)

case, dft, dmft, solver, impurity, dft_dmft = parse_config("case.toml")

@show case
@show dft
@show dmft
@show solver
@show impurity
@show dft_dmft

goodbye()
