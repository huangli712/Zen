#!/usr/bin/env julia

#
# Project : Pansy
# Source  : test.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
# Comment :
#
# Last modified: 2021/02/05
#

include("Zen.jl")
using .Zen

setup()
exhibit()
it = IterInfo()
adaptor_run(it)
