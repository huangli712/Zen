#!/usr/bin/env julia

#=

=#

include("zen.jl")
using .Zen

# check the version of julia runtime environment
check_version()

# check whether the environment variable ZEN_HOME is on  
check_home()

println("here")
exit(-1)


welcome()


toml = check_toml()

case, dft, dmft, solver, impurity, dft_dmft = parse_config(toml)

make_trees(impurity)

it = IterInfo(0, 0, 0)

dft_driver(it, dft)
@show it.dft_iter

goodbye()
