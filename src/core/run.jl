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
case, dft, dmft, dft_dmft = parse_config(toml)

# validate the parameters
# plug your codes here

# check the zen engine 
ZEN_HOME = check_home()

# check the dft engine
DFT_HOME = check_dft(dft)

# check the input files (which are essential for the dft engine)
check_inputs(dft)

# prepare the working directories
make_trees(dmft["impurity"])

# create a IterInfo object
it = IterInfo(0, 0, 0, 0)

if dft_dmft["mode"] == 1

    dft_init(it, case, dft)
    dft_run(it, dft, DFT_HOME)
    for iter in 1:dft_dmft["niter"]
        println("iter: $iter")
    end

else
    sorry()
end

goodbye()
