#
# project : pansy
# source  : config.jl
# author  : Li Huang (lihuang.dmft@gmail.com)
# status  : unstable
# comment :
#
# last modified: 2021/01/15
#

"""
    setup()
"""
function setup() end

"""
    exhibit()
"""
function exhibit() end

#
# internal functions
#

"""
    inp_toml(f::String, key::String, necessary::Bool)

Parse the configuration file (toml format)
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

Parse the configuration file (toml format)
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
    new_dict(cfg::Dict{String,Any})

Copy configurations from cfg to PCASE, PDFT, PDMFT, PIMP, and PSOLVER
"""
function new_dict(cfg::Dict{String,Any})
    # for case block
    PCASE["case"][1] = cfg["case"]

    # for dft block
    dft = cfg["dft"]
    for key in keys(dft)
        if haskey(PDFT, key)
            PDFT[key][1] = dft[key]
        else
            error("Sorry, $key is not supported currently")
        end
    end

    # for dmft block
    dmft = cfg["dmft"]
    for key in keys(dmft)
        if haskey(PDMFT, key)
            PDMFT[key][1] = dmft[key]
        else
            error("Sorry, $key is not supported currently")
        end
    end

    # for impurity block
    impurity = cfg["impurity"]
    for key in keys(impurity)
        if haskey(PIMP, key)
            PIMP[key][1] = impurity[key]
        else
            error("Sorry, $key is not supported currently")
        end
    end

    # for solver block
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

Validate the correctness and consistency of configurations
"""
function chk_dict()
    # C1. check types and existences
    #
    # check all blocks
    for P in (PCASE, PDFT, PDMFT, PIMP, PSOLVER)
        foreach(x -> _v(x.second), P)
    end

    # C2. check rationalities
    #
    # check dft block
    @assert get_d("engine") in ("vasp", "wannier")
    @assert get_d("smear") in ("m-p", "gauss", "tetra", missing)
    @assert get_d("kmesh") in ("accurate", "medium", "coarse", "file", missing)
    #
    # check dmft block
    @assert get_m("mode") in (1, 2)
    @assert get_m("axis") in (1, 2)
    @assert get_m("niter") > 0
    @assert get_m("dcount") in ("fll1", "fll2", "amf")
    @assert get_m("beta") >= 0.0
    #
    # check solver block
    @assert get_s("engine") in ("ct_hub1", "ct_hub2", "hub1", "norg")
    #
    # please add more assertion statements here

    # C3. check self-consistency
    #
    # check dft block
    if _d("lspinorb")
        @assert _d("lspins")
    end
    if _d("lproj")
        @assert !_d("lsymm") && !isa(_d("sproj"), Missing)
    end
    #
    # check solver block
    if _s("engine") in ("ct_hub1", "ct_hub2", "hub1")
        @assert _m("axis") === 1 # imaginary axis
    elseif _s("engine") in ("norg")
        @assert _m("axis") === 2 # real axis
    end
    #
    # please add more assertion statements here
end

"""
    _v(val::Array{Any,1})

Verify the value array
"""
@inline function _v(val::Array{Any,1})
    # to check if the value is updated
    if isa(val[1], Missing) && val[2] > 0
        error("Sorry, key shoule be set")
    end

    # to check if the type of value is correct
    if !isa(val[1], Missing) && !isa(val[1], eval(val[3]))
        error("Sorry, type of key is wrong")
    end
end

#
# cat_x() functions
#

"""
    cat_c()

Print the configuration parameters to stdout: for PCASE dict
"""
function cat_c()
    # see comments in cat_d()
    println("< Parameters: case >")
    println("  case     -> ", str_c("case"))
    println()
end

"""
    cat_d()

Print the configuration parameters to stdout: for PDFT dict
"""
function cat_d()
    #
    # remarks:
    #
    # get_d("sproj") is actually an Array{String,1}. it would be quite low
    # efficiency if we print it directly. so we convert it into a string
    # by using the join() function at first.
    #
    # get_d("ewidth") is actually a real number. it would be quite low
    # efficiency if we print it directly. so we convert it into a string
    # by using the string() function at first.
    #
    # see config.jl/str_d() function for more details
    #
    println("< Parameters: dft engine >")
    println("  engine   -> ", str_d("engine"))
    println("  smear    -> ", str_d("smear"))
    println("  kmesh    -> ", str_d("kmesh"))
    println("  magmom   -> ", str_d("magmom"))
    println("  lsymm    -> ", str_d("lsymm"))
    println("  lspins   -> ", str_d("lspins"))
    println("  lspinorb -> ", str_d("lspinorb"))
    println("  loptim   -> ", str_d("loptim"))
    println("  lproj    -> ", str_d("lproj"))
    println("  ewidth   -> ", str_d("ewidth"))
    println("  sproj    -> ", str_d("sproj"))
    println()
end

"""
    cat_m()

Print the configuration parameters to stdout: for PDMFT dict
"""
function cat_m()
    # see comments in cat_d()
    println("< Parameters: dmft engine >")
    println("  mode     -> ", str_m("mode"))
    println("  axis     -> ", str_m("axis"))
    println("  niter    -> ", str_m("niter"))
    println("  dcount   -> ", str_m("dcount"))
    println("  beta     -> ", str_m("beta"))
    println("  mixer    -> ", str_m("mixer"))
    println("  cc       -> ", str_m("cc"))
    println("  ec       -> ", str_m("ec"))
    println("  fc       -> ", str_m("fc"))
    println("  lcharge  -> ", str_m("lcharge"))
    println("  lenergy  -> ", str_m("lenergy"))
    println("  lforce   -> ", str_m("lforce"))
    println()
end

"""
    cat_i()

Print the configuration parameters to stdout: for PIMP dict
"""
function cat_i()
    # see comments in cat_d()
    println("< Parameters: quantum impurity atoms >")
    println("  nsite    -> ", str_i("nsite"))
    println("  atoms    -> ", str_i("atoms"))
    println("  equiv    -> ", str_i("equiv"))
    println("  shell    -> ", str_i("shell"))
    println("  ising    -> ", str_i("ising"))
    println("  occup    -> ", str_i("occup"))
    println("  upara    -> ", str_i("upara"))
    println("  jpara    -> ", str_i("jpara"))
    println("  lpara    -> ", str_i("lpara"))
    println()
end

"""
    cat_s()

Print the configuration parameters to stdout: for PSOLVER dict
"""
function cat_s()
    # see comments in cat_d()
    println("< Parameters: quantum impurity solvers >")
    println("  engine   -> ", str_s("engine"))
    println("  params   -> ", str_s("params"))
    println()
end

#
# get_x() functions
#

"""
    get_c(key::String)

Extract configurations from dict: PCASE
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

Extract configurations from dict: PDFT
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

Extract configurations from dict: PDMFT
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

Extract configurations from dict: PIMP
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

Extract configurations from dict: PSOLVER
"""
@inline function get_s(key::String)
    if haskey(PSOLVER, key)
        PSOLVER[key][1]
    else
        error("Sorry, PSOLVER does not contain key: $key")
    end
end

#
# str_x() functions
#

"""
    str_c(key::String)

Extract configurations from dict: PCASE, convert them into strings
"""
@inline function str_c(key::String)
    if haskey(PCASE, key)
        if PCASE[key][3] === :Array
            join(PCASE[key][1], "; ")
        else
            ( c = PCASE[key][1] ) isa String ? c : string(c)
        end
    else
        error("Sorry, PCASE does not contain key: $key")
    end
end

"""
    str_d(key::String)

Extract configurations from dict: PDFT, convert them into strings
"""
@inline function str_d(key::String)
    if haskey(PDFT, key)
        if PDFT[key][3] === :Array
            join(PDFT[key][1], "; ")
        else
            ( d = PDFT[key][1] ) isa String ? d : string(d)
        end
    else
        error("Sorry, PDFT does not contain key: $key")
    end
end

"""
    str_m(key::String)

Extract configurations from dict: PDMFT, convert them into strings
"""
@inline function str_m(key::String)
    if haskey(PDMFT, key)
        if PDMFT[key][3] === :Array
            join(PDMFT[key][1], "; ")
        else
            ( m = PDMFT[key][1] ) isa String ? m : string(m)
        end
    else
        error("Sorry, PDMFT does not contain key: $key")
    end
end

"""
    str_i(key::String)

Extract configurations from dict: PIMP, convert them into strings
"""
@inline function str_i(key::String)
    if haskey(PIMP, key)
        if PIMP[key][3] === :Array
            join(PIMP[key][1], "; ")
        else
            ( i = PIMP[key][1] ) isa String ? i : string(i)
        end
    else
        error("Sorry, PIMP does not contain key: $key")
    end
end

"""
    str_s(key::String)

Extract configurations from dict: PSOLVER, convert them into strings
"""
@inline function str_s(key::String)
    if haskey(PSOLVER, key)
        if PSOLVER[key][3] === :Array
            join(PSOLVER[key][1], "; ")
        else
            ( s = PSOLVER[key][1] ) isa String ? s : string(s)
        end
    else
        error("Sorry, PSOLVER does not contain key: $key")
    end
end
