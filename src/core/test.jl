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
using ProfileView

setup()
exhibit()
it = IterInfo()
@profile adaptor_run(it)
ProfileView.view()
