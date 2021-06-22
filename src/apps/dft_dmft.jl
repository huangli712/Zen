#!/usr/bin/env julia

#
# Project : Daisy
# Source  : dft_dmft.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Testing
#
# Last modified: 2021/06/22
#

#=
*Remarks 1*:

Since `ZenCore` has not been registered as a regular julia package. You
can not use the Package Manager to install it. In order to solve this
problem, you have to setup the environment variable `ZEN_CORE`, which
should point to the directory where the `ZenCore` package is downloaded
and uncompressed (in other words, `ZEN_CORE` is the directory that the
`ZenCore.jl` file is available).

*Remarks 2*:

Basically, you can execute the following codes in `REPL`. They work. 
=#

# Update LOAD_PATH
push!(LOAD_PATH, ENV["ZEN_CORE"])

# Use the ZEN Framework
using ZenCore

# R-1: Check the version of julia runtime environment
require()

# R00: Print the welcome message
welcome()

# R01: Print the overview for Zen
overview()

# R02: Setup the configuration parameters
setup()

# R03: Initialize the task
ready()

# R04: Carry out the task
go()

# R05: Finalize the task
final()

# R06: Say good bye
goodbye()
