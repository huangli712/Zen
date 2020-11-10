#!/usr/bin/env julia

#=

=#

include("zen.jl")
using .Zen

# check the version of julia runtime environment
check_version()

# check the file case.toml, which contains the configuration 
toml = check_toml()

# print the welcome message 
welcome()

# parse the file case.toml to extract parameters
case, dft, dmft, solver, impurity, dft_dmft = parse_config(toml)

# validate the parameters
# plug your codes here

# print the job's summary 
# plug your codes here

# check whether the environment variable ZEN_HOME is on  
ZEN_HOME = check_home()

# check the dft engine
DFT_HOME = check_dft(dft)

# prepare the working directories
make_trees(impurity)

it = IterInfo(0, 0, 0)

dft_driver(it, dft)
@show it.dft_iter

goodbye()
