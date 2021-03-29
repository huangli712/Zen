#
# Project : Pansy
# Source  : dmft.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/03/29
#

"""
    dmft_init(it::IterInfo, task::I64)

See also: [`dmft_exec`](@ref), [`dmft_save`](@ref).
"""
function dmft_init(it::IterInfo, task::I64)
    # Well, determine which files are necessary.
    #
    # Self-energy functions
    fsig = ("sigma.bare", "sigma.dc")
    #
    # Kohn-Sham data (including projectors) in IR format
    fir  = ("params.ir", "lattice.ir")
    #
    # Configuration file for DMFT engine
    fdmft = ("dmft.in")

    # Extract key parameters
    axis = get_m("axis")
    beta = get_m("beta")

    # Prepare essential input files, i.e., dmft.in.
    # If the `dmft.in` file exists already, it will be overwritten.
    open("dmft.in", "w") do fout
        println(fout, "task = $task")
        println(fout, "axis = $axis")
        println(fout, "beta = $beta")
    end

    # Check essential input files
    flist = (fdmft, fsig..., fir...)
    for i in eachindex(flist)
        filename = flist[i]
        if !isfile(filename)
            error("Please make sure the file $filename is available")
        end
    end
    exit(-1)
end

"""
    dmft_exec(it::IterInfo)

See also: [`dmft_init`](@ref), [`dmft_save`](@ref).
"""
function dmft_exec(it::IterInfo)
end

"""
    dmft_save(it::IterInfo)

See also: [`dmft_init`](@ref), [`dmft_exec`](@ref).
"""
function dmft_save(it::IterInfo)
end
