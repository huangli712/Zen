include("zen.jl")
using .Zen

check_version()
home = check_home()
@show home


case, dft, dmft, solver, adaptor, impurity = parse_config("case.toml")

@show case
@show dft
@show dmft
