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
case = parse_config(toml, "case") 
dft = parse_config(toml, "dft")
dmft = parse_config(toml, "dmft")
dft_dmft = parse_config(toml, "dft_dmft")
@show dft_dmft
exit(-1)

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
    for iter in 1:dft_dmft["niter"]
        println("iter: $iter")
    end

else
    sorry()
end

goodbye()
