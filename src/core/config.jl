#
# Project : Pansy
# Source  : config.jl
# Author  : Li Huang (lihuang.dmft@gmail.com)
# Status  : Unstable
#
# Last modified: 2021/04/14
#

#
# Driver Functions
#

"""
    setup()

Read parameters from configuration file, and then setup the related dicts.

See also: [`inp_toml`](@ref), [`rev_dict`](@ref), [`chk_dict`](@ref).
"""
function setup()
    # S1: Parse the case.toml file to extract configuration parameters
    cfg = inp_toml(query_args(), true)

    # S2: Revise the configuration dictionaries
    rev_dict(cfg)

    # S3: Validate the configuration parameters
    chk_dict()
end

"""
    exhibit()

Display the configuration parameters for reference.

See also: [`setup`](@ref).
"""
function exhibit()
    # E0: Print the header
    prompt("ZEN", "Configuration")

    # E1: Show dict PCASE
    println("Parameters : General")
    cat_c()

    # E2: Show dict PDFT
    println("Parameters : Density Functional Theory")
    cat_d()

    # E3: Show dict PDMFT
    println("Parameters : Dynamical Meal-Field Theory")
    cat_m()

    # E4: Show dict PIMP
    println("Parameters : Quantum Impurity Atoms")
    cat_i()

    # E5: Show dict PSOLVER
    println("Parameters : Quantum Impurity Solvers")
    cat_s()
end

#
# Service Functions (Group A)
#

"""
    inp_toml(f::String, key::String, necessary::Bool)

Parse the configuration file (toml format). It reads only parts of the
configuration file, which depends on the value of `key`.

See also: [`setup`](@ref).
"""
function inp_toml(f::String, key::String, necessary::Bool)
    if isfile(f)
        dict = TOML.parsefile(f)

        if haskey(dict, key)
            dict[key]
        else
            error("Do not have this key: $key in file: $f")
        end
    else
        if necessary
            error("Please make sure that the file $f really exists")
        else
            nothing
        end
    end
end

"""
    inp_toml(f::String, necessary::Bool)

Parse the configuration file (toml format). It reads the whole file.

See also: [`setup`](@ref).
"""
function inp_toml(f::String, necessary::Bool)
    if isfile(f)
        dict = TOML.parsefile(f)
    else
        if necessary
            error("Please make sure that the file $f really exists")
        else
            nothing
        end
    end
end

"""
    rev_dict(cfg::Dict{String,Any})

Transfer configurations from dict `cfg` to internal dicts (including
`PCASE`, `PDFT`, `PDMFT`, `PIMP`, and `PSOLVER`).

See also: [`chk_dict`](@ref).
"""
function rev_dict(cfg::Dict{String,Any})
    # For case block
    #
    # Pay attention to that the case block only includes one child element
    PCASE["case"][1] = cfg["case"]

    # For dft block
    dft = cfg["dft"]
    for key in keys(dft)
        if haskey(PDFT, key)
            PDFT[key][1] = dft[key]
        else
            error("Sorry, $key is not supported currently")
        end
    end

    # For dmft block
    dmft = cfg["dmft"]
    for key in keys(dmft)
        if haskey(PDMFT, key)
            PDMFT[key][1] = dmft[key]
        else
            error("Sorry, $key is not supported currently")
        end
    end

    # For impurity block
    impurity = cfg["impurity"]
    for key in keys(impurity)
        if haskey(PIMP, key)
            PIMP[key][1] = impurity[key]
        else
            error("Sorry, $key is not supported currently")
        end
    end

    # For solver block
    solver = cfg["solver"]
    for key in keys(solver)
        if haskey(PSOLVER, key)
            PSOLVER[key][1] = solver[key]
        else
            error("Sorry, $key is not supported currently")
        end
    end
end

"""
    chk_dict()

Validate the correctness and consistency of configurations.

See also: [`rev_dict`](@ref), [`_v`](@ref).
"""
function chk_dict()

#
# Remarks:
#
# This function is far away from completeness. We should add more
# constraints here, both physically and numerically.
#

    # C1. Check types and existences
    #
    # Check all blocks one by one
    for P in (PCASE, PDFT, PDMFT, PIMP, PSOLVER)
        foreach(x -> _v(x.first, x.second), P)
    end

    # C2. Check rationalities
    #
    # Check dft block
    @assert get_d("engine") in ("vasp",)
    @assert get_d("projtype") in ("plo", "wannier")
    @assert get_d("smear") in ("mp2", "mp1", "gauss", "tetra")
    @assert get_d("kmesh") in ("accurate", "medium", "coarse", "file")
    @assert iseven( length( get_d("window") ) )
    #
    # Check dmft block
    @assert get_m("mode") in (0, 1, 2)
    @assert get_m("axis") in (1, 2)
    @assert get_m("niter") > 0
    @assert get_m("nmesh") > 0
    @assert get_m("dcount") in ("fll1", "fll2", "amf", "exact")
    @assert get_m("beta") >= 0.0
    @assert get_m("mixer") >= 0.0 && get_m("mixer") <= 1.0
    #
    # Check impurity block
    @assert get_i("nsite") >= 1 && get_i("nsite") <= 99
    #
    # Check solver block
    @assert get_s("engine") in ("ct_hyb1", "ct_hyb2", "hub1", "norg")
    #
    # Please add more assertion statements here

    # C3. Check self-consistency
    #
    # Check dft block
    @assert get_d("projtype") === "plo"
    if get_d("lspinorb")
        @assert !get_d("lspins")
    end
    if get_d("lproj")
        @assert !get_d("lsymm") && !isa(get_d("sproj"), Missing)
    end
    #
    # Check solver block
    if get_s("engine") in ("ct_hyb1", "ct_hyb2", "hub1")
        @assert get_m("axis") === 1 # Imaginary axis
    elseif get_s("engine") in ("norg",)
        @assert get_m("axis") === 2 # Real axis
    end
    #
    # Please add more assertion statements here
end

#
# Service Functions (Group B)
#

"""
    _v(key::String, val::Array{Any,1})

Verify the value array. Called by chk_dict() function only.

See also: [`chk_dict`](@ref).
"""
@inline function _v(key::String, val::Array{Any,1})
    # To check if the value is updated
    if isa(val[1], Missing) && val[2] > 0
        error("Sorry, key ($key) shoule be set")
    end

    # To check if the type of value is correct
    if !isa(val[1], Missing) && !isa(val[1], eval(val[3]))
        error("Sorry, type of key ($key) is wrong")
    end
end

#
# Service Functions (Group C: CAT Functions)
#

"""
    cat_c()

Print the configuration parameters to stdout: for PCASE dict.

See also: [`get_c`](@ref), [`str_c`](@ref).
"""
function cat_c()
    # See comments in cat_d()
    println("  case     | case     -> ", str_c("case"))
    println()
end

"""
    cat_d()

Print the configuration parameters to stdout: for PDFT dict.

See also: [`get_d`](@ref), [`str_d`](@ref).
"""
function cat_d()

#
# Remarks:
#
# The get_d("sproj") is actually an Array{String,1}. It would be quite
# low efficiency if we print it directly. So we have to convert it into
# a string by using the join() function at first.
#
# The get_d("window") is actually an Array{F64,1}. It would be quite
# low efficiency if we print it directly. So we have to convert it into
# a string by using the join() function at first.
#
# See config.jl/str_d() function for more details.
#

    println("  dft      | engine   -> ", str_d("engine"))
    println("  dft      | projtype -> ", str_d("projtype"))
    println("  dft      | smear    -> ", str_d("smear"))
    println("  dft      | kmesh    -> ", str_d("kmesh"))
    println("  dft      | magmom   -> ", str_d("magmom"))
    println("  dft      | lsymm    -> ", str_d("lsymm"))
    println("  dft      | lspins   -> ", str_d("lspins"))
    println("  dft      | lspinorb -> ", str_d("lspinorb"))
    println("  dft      | loptim   -> ", str_d("loptim"))
    println("  dft      | lproj    -> ", str_d("lproj"))
    println("  dft      | sproj    -> ", str_d("sproj"))
    println("  dft      | window   -> ", str_d("window"))
    println()
end

"""
    cat_m()

Print the configuration parameters to stdout: for PDMFT dict.

See also: [`get_m`](@ref), [`str_m`](@ref).
"""
function cat_m()
    # See comments in cat_d()
    println("  dmft     | mode     -> ", str_m("mode"))
    println("  dmft     | axis     -> ", str_m("axis"))
    println("  dmft     | niter    -> ", str_m("niter"))
    println("  dmft     | nmesh    -> ", str_m("nmesh"))
    println("  dmft     | dcount   -> ", str_m("dcount"))
    println("  dmft     | beta     -> ", str_m("beta"))
    println("  dmft     | mixer    -> ", str_m("mixer"))
    println("  dmft     | mc       -> ", str_m("mc"))
    println("  dmft     | cc       -> ", str_m("cc"))
    println("  dmft     | ec       -> ", str_m("ec"))
    println("  dmft     | fc       -> ", str_m("fc"))
    println("  dmft     | lfermi   -> ", str_m("lfermi"))
    println("  dmft     | lcharge  -> ", str_m("lcharge"))
    println("  dmft     | lenergy  -> ", str_m("lenergy"))
    println("  dmft     | lforce   -> ", str_m("lforce"))
    println()
end

"""
    cat_i()

Print the configuration parameters to stdout: for PIMP dict.

See also: [`get_i`](@ref), [`str_i`](@ref).
"""
function cat_i()
    # See comments in cat_d()
    println("  impurity | nsite    -> ", str_i("nsite"))
    println("  impurity | atoms    -> ", str_i("atoms"))
    println("  impurity | equiv    -> ", str_i("equiv"))
    println("  impurity | shell    -> ", str_i("shell"))
    println("  impurity | ising    -> ", str_i("ising"))
    println("  impurity | occup    -> ", str_i("occup"))
    println("  impurity | upara    -> ", str_i("upara"))
    println("  impurity | jpara    -> ", str_i("jpara"))
    println("  impurity | lpara    -> ", str_i("lpara"))
    println()
end

"""
    cat_s()

Print the configuration parameters to stdout: for PSOLVER dict.

See also: [`get_s`](@ref), [`str_s`](@ref).
"""
function cat_s()
    # See comments in cat_d()
    println("  solver   | engine   -> ", str_s("engine"))
    println("  solver   | params   -> ", str_s("params"))
    println()
end

#
# Service Functions (Group D: GET Functions)
#

"""
    get_c(key::String)

Extract configurations from dict: PCASE.

See also: [`cat_c`](@ref), [`str_c`](@ref).
"""
@inline function get_c(key::String)
    if haskey(PCASE, key)
        PCASE[key][1]
    else
        error("Sorry, PCASE does not contain key: $key")
    end
end

"""
    get_d(key::String)

Extract configurations from dict: PDFT.

See also: [`cat_d`](@ref), [`str_d`](@ref).
"""
@inline function get_d(key::String)
    if haskey(PDFT, key)
        PDFT[key][1]
    else
        error("Sorry, PDFT does not contain key: $key")
    end
end

"""
    get_m(key::String)

Extract configurations from dict: PDMFT.

See also: [`cat_m`](@ref), [`str_m`](@ref).
"""
@inline function get_m(key::String)
    if haskey(PDMFT, key)
        PDMFT[key][1]
    else
        error("Sorry, PDMFT does not contain key: $key")
    end
end

"""
    get_i(key::String)

Extract configurations from dict: PIMP.

See also: [`cat_i`](@ref), [`str_i`](@ref).
"""
@inline function get_i(key::String)
    if haskey(PIMP, key)
        PIMP[key][1]
    else
        error("Sorry, PIMP does not contain key: $key")
    end
end

"""
    get_s(key::String)

Extract configurations from dict: PSOLVER.

See also: [`cat_s`](@ref), [`str_s`](@ref).
"""
@inline function get_s(key::String)
    if haskey(PSOLVER, key)
        PSOLVER[key][1]
    else
        error("Sorry, PSOLVER does not contain key: $key")
    end
end

#
# Service Functions (Group E: STR Functions)
#

"""
    str_c(key::String)

Extract configurations from dict: PCASE, convert them into strings.

See also: [`cat_c`](@ref), [`get_c`](@ref).
"""
@inline function str_c(key::String)
    if haskey(PCASE, key)
        c = PCASE[key][1]
        if PCASE[key][3] === :Array && !isa(c, Missing)
            join(c, "; ")
        else
            c isa String ? c : string(c)
        end
    else
        error("Sorry, PCASE does not contain key: $key")
    end
end

"""
    str_d(key::String)

Extract configurations from dict: PDFT, convert them into strings.

See also: [`cat_d`](@ref), [`get_d`](@ref).
"""
@inline function str_d(key::String)
    if haskey(PDFT, key)
        d = PDFT[key][1]
        if PDFT[key][3] === :Array && !isa(d, Missing)
            join(d, "; ")
        else
            d isa String ? d : string(d)
        end
    else
        error("Sorry, PDFT does not contain key: $key")
    end
end

"""
    str_m(key::String)

Extract configurations from dict: PDMFT, convert them into strings.

See also: [`cat_m`](@ref), [`get_m`](@ref).
"""
@inline function str_m(key::String)
    if haskey(PDMFT, key)
        m = PDMFT[key][1]
        if PDMFT[key][3] === :Array && !isa(m, Missing)
            join(m, "; ")
        else
            m isa String ? m : string(m)
        end
    else
        error("Sorry, PDMFT does not contain key: $key")
    end
end

"""
    str_i(key::String)

Extract configurations from dict: PIMP, convert them into strings.

See also: [`cat_i`](@ref), [`get_i`](@ref).
"""
@inline function str_i(key::String)
    if haskey(PIMP, key)
        i = PIMP[key][1]
        if PIMP[key][3] === :Array && !isa(i, Missing)
            join(i, "; ")
        else
            i isa String ? i : string(i)
        end
    else
        error("Sorry, PIMP does not contain key: $key")
    end
end

"""
    str_s(key::String)

Extract configurations from dict: PSOLVER, convert them into strings.

See also: [`cat_s`](@ref), [`get_s`](@ref).
"""
@inline function str_s(key::String)
    if haskey(PSOLVER, key)
        s = PSOLVER[key][1]
        if PSOLVER[key][3] === :Array && !isa(s, Missing)
            join(PSOLVER[key][1], "; ")
        else
            s isa String ? s : string(s)
        end
    else
        error("Sorry, PSOLVER does not contain key: $key")
    end
end
