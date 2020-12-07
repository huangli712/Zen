#!/usr/bin/env julia

include("Zen.jl")
using .Zen

#kmesh, weight = from_ibzkpt(pwd() * "/dft")
#
#for i in eachindex(weight)
#    println(i, " ", kmesh[i,:], " ", weight[i])
#end

#from_poscar(pwd() * "/dft")
from_eigenval(pwd() * "/dft")
