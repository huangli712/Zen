#!/usr/bin/env julia

#
# Project : Begonia
# Source  : run.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/02/28
#

# Update LOAD_PATH
push!(LOAD_PATH, ENV["ZEN_CORE"])

# Use the ZEN Framework
using Zen

# R-1: Check the version of julia runtime environment
require()

# R00: Print the welcome message
welcome()

# R01: Print the overview for Zen
overview()

# R02: Setup the configuration parameters
setup()

# R03: Print the configuration parameters
exhibit()

# R04: Initialize the task
ready()

# R05: Carry out the task
go()

# R06: Finalize the task
final()

# R07: Say good bye
goodbye()
