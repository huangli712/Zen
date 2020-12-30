#!/usr/bin/env julia

#
# project : pansy
# source  : run.jl
# author  : Li Huang (lihuang.dmft@gmail.com)
# status  : unstable
# comment :
#
# last modified: 2020/12/30
#

include("Zen.jl")
using .Zen

# S-1: check the version of julia runtime environment
require()

# S00: print the welcome message
welcome()

# S01: print the overview of zen
message("ZEN", "OVERVIEW")
overview()

# S02: parse the configuration file, get job's description
message("ZEN", "PARSER")
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
message("ZEN", "PARAMETER")
#
# S03.1: show dict PCASE
list_case()
#
# S03.2: show dict PDFT
list_dft()
#
# S03.3: show dict PDMFT
list_dmft()
#
# S03.4: show dict PIMP
list_impurity()
#
# S03.5: show dict PSOLVER
list_solver()

# S04: initialize the job 
message("ZEN", "PREPARATION")
#
# S04.1: check the input files (which are essential for the calculation)
query_inps()
println("Check essential input files")
#
# S04.2: prepare the working directories
make_trees()
println("Create working directories")
#
# S04.3: create an IterInfo object
it = IterInfo()
println("Create self-consistent iterator")
#
# S04.4: check self-consistent mode
mode = _m("mode") === 1 ? "one-shot" : "fully self-consistent"
println("Check self-consistent mode: $mode\n")

message("ZEN", "START")

# well, we choose the one-shot mode
if _m("mode") === 1

#
# remarks:
# 
# we would like to perform two successive dft runs if _d("loptim") is true.
# the purpose of the first dft run is to evaluate the fermi level. then a
# energy window is determined. we will use this window to generate optimal
# projectors in the second dft run. On the other hand, if _d("loptim") is
# false, only the first dft run is enough.
#

    # S05: perform dft calculation (for the first time).
    message("ZEN", "DFT")
    #
    # S05.1: prepare and check essential files for the dft engine 
    dft_init(it)
    println("Initialize everything needed by the dft engine")
    #
    # S05.2: perform a self-consitent calculation at the dft level
    dft_run(it)
    println("Launch the dft engine")
    #
    # S05.3: backup the output files of the dft engine
    dft_save(it)
    println("Save the output files\n")

    # we want better optimal projectors
    # in the previous dft run, initial fermi level = 0 -> wrong energy
    # window -> wrong optimial projectors. but at this point, the fermi
    # level is updated, so we have to generate the optimal projectors
    # again within this window
    if _d("loptim")

        # S06: perform dft calculation (for the second time).
        message("ZEN", "DFT")
        #
        # S06.1: prepare and check essential files for the dft engine 
        dft_init(it)
        println("Initialize everything needed by the dft engine")
        #
        # S06.2: perform a self-consitent calculation at the dft level
        dft_run(it)
        println("Launch the dft engine")
        #
        # S06.3: backup the output files of the dft engine
        dft_save(it)
        println("Save the output files\n")

    end

    # S07:
    message("ZEN", "ADAPTOR")
    adaptor_init(it)
    adaptor_run(it)
    adaptor_save(it)

else
    sorry()
end

goodbye()
