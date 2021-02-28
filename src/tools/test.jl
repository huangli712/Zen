#!/usr/bin/env julia

#
# Project : Begonia
# Source  : test.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
# Comment :
#
# Last modified: 2021/02/28
#

# Update LOAD_PATH
push!(LOAD_PATH, joinpath(ENV["ZEN_HOME"], "/src/core"))

# Use the ZEN Framework
include("Zen.jl")
using Zen

# Put your codes here
setup()
exhibit()
it = IterInfo()
adaptor_run(it)
