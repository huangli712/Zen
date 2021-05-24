#!/usr/bin/env julia

#
# Project : Begonia
# Source  : test.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/05/24
#

# Update LOAD_PATH
push!(LOAD_PATH, ENV["ZEN_CORE"])

# Use the ZEN Framework
using ZenCore

# Put your codes here
setup()
exhibit()
it = IterInfo()
lr = Logger(query_case())
adaptor_run(it, lr)
