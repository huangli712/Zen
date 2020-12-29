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

# S01: print the overview of zen
message("ZEN", "Overview")
overview()

# S02: parse the configuration file, get job's description
message("ZEN", "Parsing Job")
#
# S02.1: parse the case.toml file to extract configuration parameters
cfg = parse_toml(query_args(), true)
println("Extract configuration parameters from $(query_args())")
#
# S02.2: build the configuration dictionaries
renew_config(cfg)
println("Update configuration dictionaries")
#
# S02.3: validate the configuration dictionaries
check_config()
println("Verify configuration dictionaries\n")

# S03: write the configuration parameters to stdout
message("ZEN", "Viewing Job")
#
# S03.1: show dict PCASE
view_case()
#
# S03.2: show dict PDFT
view_dft()
#
# S03.3: show dict PDMFT
view_dmft()
#
# S03.4: show dict PIMP
view_impurity()
#
# S03.5: show dict PSOLVER
view_solver()

# S04: initialize the job 
message("ZEN", "Initializing")
#
# S04.1: check the input files (which are essential for the calculation)
query_inps()
println("Check essential input files")
#
# S04.2: prepare the working directories
make_trees()
println("Create working directories")
#
# S04.3: create a IterInfo object
it = IterInfo()
println("Make self-consistent iterator\n")

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
