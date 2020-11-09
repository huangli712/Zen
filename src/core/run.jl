
include("zen.jl")
using .Zen

check_version()
case, dft, dmft, solver, adaptor, impurity = parse_config("case.toml")
