#!/usr/bin/env julia

#
# project : pansy
# source  : run.jl
# author  : Li Huang (lihuang.dmft@gmail.com)
# status  : unstable
# comment :
#
# last modified: 2021/01/19
#

include("Zen.jl")
using .Zen

# S-1: check the version of julia runtime environment
require()

# S00: print the welcome message
welcome()

# S01: print the overview for zen
prompt("ZEN", "OVERVIEW")
overview()

# S02: parse the configuration file, get job's description
prompt("ZEN", "PARSER")
setup()

# S03: write the configuration parameters to stdout
prompt("ZEN", "VIEWER")
exhibit()

# S04: initialize the job
prompt("ZEN", "CREATOR")
ready()

prompt("ZEN", "START")
go()
exit(-1)

goodbye()
