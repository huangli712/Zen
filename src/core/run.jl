#!/usr/bin/env julia

#=

=#

include("zen.jl")
using .Zen

# check the version of julia runtime environment
check_version()

# print the welcome message 
welcome()

# parse the file case.toml to extract parameters
case = parse_config(check_toml(), "case") 
dft = parse_config(check_toml(), "dft")
dmft = parse_config(check_toml(), "dmft")
dft_dmft = parse_config(check_toml(), "dft_dmft")

# validate the parameters
# plug your codes here

# check the input files (which are essential for the dft engine)
check_inputs(dft)

# prepare the working directories
make_trees(dmft["impurity"])

# create a IterInfo object
it = IterInfo(0, 0, 0, 0)

if dft_dmft["mode"] == 1

    dft_init(it, case, dft)
    dft_run(it, dft)
    dft_save(it, dft)
    for iter in 1:dft_dmft["niter"]
        println("iter: $iter")
    end

else
    sorry()
end

goodbye()
