#!/usr/bin/env julia

#
# Project : Begonia
# Source  : ploband.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Testing
#
# Last modified: 2021/11/08
#

#=
*Remarks*:

This script is used to generate the DFT band structures via the
projected local orbitals scheme. Only for debug purpose.
=#

# Update LOAD_PATH
push!(LOAD_PATH, ENV["ZEN_CORE"])

# Using standard library
using Printf

# Using the ZenCore library
using ZenCore

# Read the Kohn-Sham dataset
println("Parse the Kohn-Sham dataset")
D = ir_read("dft")

# Check the validity of the `D` dict
key_list = [:MAP, :PG, :PW,
            :latt, :kmesh, :weight,
            :enk, :occupy, :Fchipsi, :fermi]
for k in key_list
    @assert haskey(D, k)
end

# Generate 𝑟-points in Wigner-Seitz cell
println("Generate the Wigner-Seitz cell")
rdeg, rvec = w90_make_rcell(D[:latt])

# Build high-symmetry 𝑘-path
println("Generate the high-symmetry 𝑘-path in the Brillouin zone")
#
# Number of 𝑘-points per direction
ndiv = 100
#
# Please modify the following 𝑘-points to define high-symmetry 𝑘-paths
kstart = [0.0 0.0 0.0; # Γ
          0.5 0.0 0.0; # X
          0.5 0.5 0.0; # M
          0.0 0.0 0.0] # Γ
kend   = [0.5 0.0 0.0; # X
          0.5 0.5 0.0; # M
          0.0 0.0 0.0; # Γ
          0.5 0.5 0.5] # R
#
# Generate 𝑘-list
kpath, xpath = w90_make_kpath(ndiv, kstart, kend)

# Build the hamiltonian in an uniform 𝑘-mesh
println("Restore the Kohn-Sham Hamiltonian")
hamk = calc_hamk(D[:PW], D[:Fchipsi], D[:enk])

# Get nspin
_, _, nspin = size(D[:enk])

# Go through the groups and spins
for p in eachindex(hamk)
    for s = 1:nspin
        println("Group: [$p] Spin: [$s]")

        # Calculate H(𝑟)
        println("  > Generate H(R)")
        HR = w90_make_hamr(D[:kmesh], rvec, hamk[p][:,:,:,s])

        # Build H(𝑘) along high-symmetry directions
        println("  > Generate H(K)")
        HK = w90_make_hamk(kpath, rdeg, rvec, HR)

        # Calculate the band structures
        println("  > Diagonalize H(K)")
        eigs, evec = w90_diag_hamk(HK)
        nband, nkpt = size(eigs)

        # Dump the band structures
        println("  > Dump band structures into band.dat.h$p.s$s")
        open("band.dat.h$p.s$s", "w") do fout
            for b = 1:nband
                for k = 1:nkpt
                    println(fout, xpath[k], " ", eigs[b,k])
                end
                println(fout)
            end
        end # END OF IOSTREAM
    end # END OF S LOOP
end # END OF P LOOP

# Dump the 𝑘-list
nkpt, _ = size(kpath)
open("kpath.dat", "w") do fout
    for k = 1:nkpt
        @printf(fout, "%12.6f %8.6f %8.6f %6.4f\n", kpath[k,:]..., 1.00)
    end
end
