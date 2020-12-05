#!/usr/bin/env julia

include("Zen.jl")
using .Zen

kmesh, weight = from_ibzkpt(pwd() * "/dft")

for i in length(kmesh)
    println(i, " ", kmesh[i,:], " ", weight[i])
end
