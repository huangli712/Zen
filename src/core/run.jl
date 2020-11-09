include("zen.jl")
using .Zen

welcome()

check_version()

home = check_home()

case, dft, dmft, solver, impurity, dft_dmft = parse_config("case.toml")

@show case
@show dft
@show dmft

goodbye()
