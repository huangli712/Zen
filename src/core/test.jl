#!/usr/bin/env julia

#
# Project : Pansy
# Source  : test.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
# Comment :
#
# Last modified: 2021/02/03
#

include("Zen.jl")
using .Zen
using Profile
using BenchmarkTools

setup()
exhibit()
it = IterInfo()
@timev adaptor_run(it)
