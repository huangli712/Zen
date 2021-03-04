#!/usr/bin/env julia

#
# Project : Begonia
# Source  : analyze.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
# Comment :
#
# Last modified: 2021/03/04
#

# Update LOAD_PATH
push!(LOAD_PATH, ENV["ZEN_CORE"])

# Use the ZEN Framework
using Zen
using Printf

# Define orbital labels
orb_labels = ["s",
              "py", "pz", "px",
              "dxy", "dyz", "dz2", "dxz", "dx2-y2",
              "fz3", "fxz2", "fyz2", "fz(x2-y2)", "fxyz", "fx(x2-3y2)", "fy(3x2-y2)"]

# Parse the PROCAR file if it is available
print("Please specify the folder that contains the PROCAR file: ")
path = readline(stdin)
oab, enk, occ = vaspio_procar(path)

# Extract key parameters
norbs, natom, nband, nspin = size(oab)
_, nkpt, _ = size(enk)

# The data are ready.
# Then this function will interact with the users.
# Print essential information
println()
println("Number of spins: $nspin")
println("Number of k-points: $nkpt")
println("Number of bands: $nband")
println("Number of atoms: $natom")
println("Number of atomic orbitals: $norbs")

# Enter a infinite loop until the users enter `q`
while true
    # Get spin index
    print("Please input spin index (integer, from 1 to $nspin): ")
    spin_index = parse(I64, readline(stdin))

    # Get atom index
    print("Please input atom index (integer, from 1 to $natom): ")
    atom_index = parse(I64, readline(stdin))

    # Get atomic orbital index
    println("Atomic orbitals: ", orb_labels)
    print("Please input atomic orbital index (integer, from 1 to $norbs): ")
    orbital_index = parse(I64, readline(stdin))

    # How many bands would you like to see?
    print("How many bands would you like to see (integer, 1, 3, 5, or 7): ")
    nview = parse(I64, readline(stdin))
    @assert nview < nband

    # Output the gathered information
    println("Selected spin index: $spin_index")
    println("Selected atom index: $atom_index")
    println("Selected atomic orbital index: $orbital_index")
    println("Selected atomic orbital label: $(orb_labels[orbital_index])")

    # Sort, find out the most relevant orbitals
    v = sortperm(oab[orbital_index, atom_index, :, spin_index], rev = true)

    # Output the band indices and weights
    print("Band index :")
    foreach(x -> @printf("%12i", x), v[1:nview])
    println()
    print("Band weight:")
    foreach(x -> @printf("%12.7f", x), oab[orbital_index, atom_index, v[1:nview], spin_index])
    println()

    # Prompt whether the users want to continue or quit
    println("If you want to continue, please enter `c` key, or else press `q` key")
    q = readline(stdin)

    # Quit the loop
    if q === "q"
        break
    end
end

println("Haha")
