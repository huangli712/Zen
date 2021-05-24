#!/usr/bin/env julia

#
# Project : Begonia
# Source  : analyze.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/05/24
#

#
# Remarks:
#
# This script is used to parse the PROCAR file. The users can use it to
# figure out which bands are the most relevant and the corresponding
# energy window for these bands. These information is quite useful. Then
# we can determine the `window` entry in the `case.toml` file.  
#

# Update LOAD_PATH
push!(LOAD_PATH, ENV["ZEN_CORE"])

# Use the ZEN Framework
using ZenCore

# Define orbital labels
orb_labels = ["1:s",
              "2:py", "3:pz", "4:px",
              "5:dxy", "6:dyz", "7:dz2", "8:dxz", "9:dx2-y2",
              "10:fz3", "11:fxz2", "12:fyz2", "13:fz(x2-y2)", "14:fxyz", "15:fx(x2-3y2)", "16:fy(3x2-y2)"]

# Parse the PROCAR file if it is available
print("Please specify the folder that contains the PROCAR file: ")
path = readline(stdin)
oab, enk, occ = vaspio_procar(path)

# Extract key parameters
norbs, natom, nband, nspin = size(oab)
_, nkpt, _ = size(enk)

# Now the data are ready.
# Then this function will interact with the users.
# Print essential information
println()
println("Number of spins: $nspin")
println("Number of k-points: $nkpt")
println("Number of bands: $nband")
println("Number of atoms: $natom")
println("Number of atomic orbitals: $norbs")
println()

# To tell us which bands are relevant
# Enter an infinite loop until the users enter `q`
println("Next we will enter an infinite loop to find out which bands are relevant")
while true
    # Get spin index
    print("Please input spin index (integer, from 1 to $nspin): ")
    spin_index = parse(I64, readline(stdin))

    # Get atom index
    print("Please input atom index (integer, from 1 to $natom): ")
    atom_index = parse(I64, readline(stdin))

    # Get atomic orbital index
    println("Atomic orbitals: ", orb_labels[1:norbs])
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
    println("The following bands are relevant:")
    print("band index :")
    foreach(x -> @printf("%12i", x), v[1:nview])
    println()
    print("band weight:")
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

println()
println("Next we will try to determine the energy window for the selected bands")
while true
    # Get spin index
    print("Please input spin index (integer, from 1 to $nspin): ")
    spin_index = parse(I64, readline(stdin))

    # Get fermi level
    print("Please input fermi level (float): ")
    fermi = parse(F64, readline(stdin))

    # Get band index
    print("Please input index for the low-lying band (integer, from 1 to $nband):")
    low_band_index = parse(I64, readline(stdin))

    # Get band index
    print("Please input index for the high-lying band (integer, from 1 to $nband):")
    high_band_index = parse(I64, readline(stdin))

    # Sanity check
    @assert low_band_index >= 1 && low_band_index <= nband
    @assert high_band_index >= 1 && high_band_index <= nband
    @assert low_band_index <= high_band_index

    # Evaluate the energy window for the selected band window
    min_ene = minimum(enk[low_band_index:high_band_index, :, spin_index])
    max_ene = maximum(enk[low_band_index:high_band_index, :, spin_index])
    println("The energy window for the selected band window is:")
    println("band window: $low_band_index -> $high_band_index")
    println("energy window: $min_ene -> $max_ene")
    println("energy window: $(min_ene - fermi) -> $(max_ene - fermi) (adjusted)")

    # Prompt whether the users want to continue or quit
    println("If you want to continue, please enter `c` key, or else press `q` key")
    q = readline(stdin)

    # Quit the loop
    if q === "q"
        break
    end
end
