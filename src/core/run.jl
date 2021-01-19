#!/usr/bin/env julia

#
# Project : Pansy
# Source  : run.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
# Comment :
#
# Last modified: 2021/01/19
#

include("Zen.jl")
using .Zen

# S-1: Check the version of julia runtime environment
require()

# S00: Print the welcome message
welcome()

# S01: Print the overview for Zen
prompt("ZEN", "OVERVIEW")
overview()

# S02: Setup the configuration parameters
prompt("ZEN", "PARSER")
setup()

# S03: Print the configuration parameters to stdout
prompt("ZEN", "VIEWER")
exhibit()

# S04: Initialize the task
prompt("ZEN", "CREATOR")
ready()

# S05: Do the calculations
prompt("ZEN", "START")
go()

# S06: Finalize the calculations
final()

# S07: Say good bye
goodbye()
