#!/usr/bin/env julia

#=

=#

include("zen.jl")
using .Zen

# check the version of julia runtime environment
check_version()

# check whether the environment variable ZEN_HOME is on  
ZEN_HOME = check_home()

# check the file case.toml, which contains the configuration 
toml = check_toml()

welcome()



case, dft, dmft, solver, impurity, dft_dmft = parse_config(toml)

make_trees(impurity)

it = IterInfo(0, 0, 0)

dft_driver(it, dft)
@show it.dft_iter

goodbye()
