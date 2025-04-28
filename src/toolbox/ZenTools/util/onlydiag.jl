#!/usr/bin/env julia

#
# Project : Begonia
# Source  : onlydiag.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Testing
#
# Last modified: 2025/04/03
#

#=
*Remarks*:

This script is used to parse the matrix functions, which are essential
output of the dmft engine (`Dyson`), and extract the diagonal elements
only. Only for debug purpose.
=#

# Update LOAD_PATH
haskey(ENV,"ZEN_CORE") && pushfirst!(LOAD_PATH, ENV["ZEN_CORE"])

# Using standard library
using Printf

# Using the ZenCore library
using ZenCore

# Please setup what you want to do
#
# Total configuration file
fcfg = "SrVO3.toml"
#
# File that contains original matrix functions
fmat = "dmft1/dmft.delta"
#
# File that contains the diagonal elements
fdia = fmat * ".diag"

# Parse the configuration
cfg = inp_toml(fcfg, true)
rev_dict(cfg)
ai = GetImpurity()

# Read the matrix functions
fmesh, Delta = read_delta(ai, fmat)
_, qdim, nmesh, nspin, nsite = size(Delta)

# Write the diagonal elements
for t = 1:nsite
    for s = 1:nspin
        open(fdia * ".$s.$t", "w") do fout
            for m = 1:nmesh
                @printf(fout, "%6i%16.8f", m, fmesh[m])
                for q = 1:ai[t].nband
                    z = Delta[q,q,m,s,t]
                    @printf(fout, "%16.8f%16.8f", real(z), imag(z))
                end # END OF Q LOOP
                println(fout)
            end # END OF M LOOP
        end # END OF IOSTREAM
    end # END OF S LOOP
end # END OF T LOOP
