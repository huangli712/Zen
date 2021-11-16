#!/usr/bin/env julia

#
# Project : Begonia
# Source  : ircheck.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Testing
#
# Last modified: 2021/11/16
#

#=
*Remarks*:

Are the projectors correct? It is an important question. Here we will
try to calculate some physical quantitites, which will be written to
external files or terminal for reference.

These physical quantities include density matrix, overlap matrix, local
hamiltonian, full hamiltonian, and partial density of states. Of course,
it is time-comsuming to do these things.

This script is only for debug purpose.
=#

# Update LOAD_PATH
push!(LOAD_PATH, ENV["ZEN_CORE"])

# Using standard library
using Printf

# Using the ZenCore library
using ZenCore

# Read the Kohn-Sham dataset
D = ir_read("dft")

# Check the validity of the `D` dict
key_list = [:MAP, :PG, :PW,
            :latt, :kmesh, :weight,
            :enk, :occupy, :Fchipsi, :fermi, :chipsi]
for k in key_list
    @assert haskey(D, k)
end

# Check the projectors
plo_check(D)
