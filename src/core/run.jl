#!/usr/bin/env julia

#
# Project : Pansy
# Source  : run.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
# Comment :
#
# Last modified: 2021/01/20
#

include("Zen.jl")
using .Zen

# S-1: Check the version of julia runtime environment
require()

# S00: Print the welcome message
welcome()

# S01: Print the overview for Zen
overview()

# S02: Setup the configuration parameters
setup()

# S03: Print the configuration parameters to stdout
exhibit()

# S04: Initialize the task
ready()

# S05: Do the calculations
prompt("ZEN", "START")
go()

# S06: Finalize the calculations
final()

# S07: Say good bye
goodbye()
