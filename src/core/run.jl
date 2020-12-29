#!/usr/bin/env julia

#
# project : pansy
# source  : run.jl
# author  : Li Huang (lihuang.dmft@gmail.com)
# status  : unstable
# comment :
#
# last modified: 2020/12/29
#

include("Zen.jl")
using .Zen

# S-1: check the version of julia runtime environment
require()

# S00: print the welcome message
welcome()

# S01: parse the configuration file, get job's description
message("ZEN", "Parsing Job")
#
# S01.1: parse the case.toml file to extract configuration parameters
cfg = parse_toml(query_args(), true)
println("Extract configuration parameters from $(query_args())")
#
# S01.2: build the configuration dictionaries
renew_config(cfg)
println("Update configuration dictionaries")
#
# S01.3: validate the configuration dictionaries
check_config()
println("Check configuration dictionaries\n")

# S02: write the configuration parameters to stdout
message("ZEN", "Viewing Job")
#
# S02.1: show dict PCASE
view_case()
#
# S02.2: show dict PDFT
view_dft()
#
# S02.3: show dict PDMFT
view_dmft()
#
# S02.4: show dict PIMP
view_impurity()
#
# S02.5: show dict PSOLVER
view_solver()

# check the input files (which are essential for the calculation)
message("ZEN", "Preparing Job")
query_inps()

# prepare the working directories
make_trees()

# create a IterInfo object
it = IterInfo()

exit(-1)

if _m("mode") === 1

    message("zen", "enter one-shot mode")
    message("zen", "begin < dft block >")
    message("zen", "dft -> init")
    dft_init(it)
    message("zen", "dft -> run")
    dft_run(it)
    message("zen", "dft -> save")
    message("zen", "e_n_d < dft block >")
    dft_save(it)
    for iter = 1:_m("niter")
        message("zen", "dmft_cycle -> 0  dmft1_iter -> $iter dmft2_iter -> 0")
    end

else
    sorry()
end

goodbye()
