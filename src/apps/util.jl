#
# project : pansy
# source  : util.jl
# author  : Li Huang (lihuang.dmft@gmail.com)
# status  : unstable
# comment :
#
# last modified: 2020/12/15
#

"""
    require()

Check the version of julia runtime environment
"""
function require()
    if VERSION < v"1.5-"
        error("Please upgrade your julia to v1.5.0 or higher")
    end
end

"""
    query_args()

Check whether the configuration file (case.toml) is provided
"""
function query_args()
    nargs = length(ARGS)
    if nargs < 1
        error("Please specify the configuration file at least")
    else
        ARGS[1]
    end
end

"""
    query_cars()

Check whether the essential input files exist
"""
function query_cars()
    if _d("engine") === "vasp"
        if !isfile("POSCAR") || !isfile("POTCAR")
            error("Please provide both POSCAR and POTCAR files")
        end
    else
        sorry()
    end
end

"""
    query_zen()

Query the home directory for zen
"""
function query_zen()
    if haskey(ENV, "ZEN_HOME")
        ENV["ZEN_HOME"]
    else
        error("ZEN_HOME is undefined")
    end
end

"""
    query_dft()

Query the home directory of the dft engine
"""
function query_dft()
    if _d("engine") === "vasp"
        if haskey(ENV, "VASP_HOME")
            ENV["VASP_HOME"]
        else
            error("VASP_HOME must be defined")
        end
    else
        sorry()
    end
end

"""
    view_case()

Print the configuration parameters to stdout: for case
"""
function view_case()
    println()
    message("case summary")
    println("case -> ", _c("case"))
    println()
end

"""
    view_dft()

Print the configuration parameters to stdout: for dft
"""
function view_dft()
    message("dft parameters")
    println("dft  -> engine   -> ", _d("engine"))
    println("dft  -> smear    -> ", _d("smear"))
    println("dft  -> kmesh    -> ", _d("kmesh"))
    println("dft  -> magmom   -> ", _d("magmom"))
    println("dft  -> lsymm    -> ", _d("lsymm"))
    println("dft  -> lspins   -> ", _d("lspins"))
    println("dft  -> lspinorb -> ", _d("lspinorb"))
    println("dft  -> window   -> ", _d("window"))
    println("dft  -> loptim   -> ", _d("loptim"))
    println("dft  -> lproj    -> ", _d("lproj"))
    println("dft  -> nproj    -> ", _d("nproj"))
    println("dft  -> sproj    -> ", _d("sproj"))
    println()
end

"""
    view_dmft()

Print the configuration parameters to stdout: for dmft
"""
function view_dmft()
    message("dmft parameters")
    println("dmft -> mode     -> ", _m("mode"))
    println("dmft -> axis     -> ", _m("axis"))
    println("dmft -> beta     -> ", _m("beta"))
    println("dmft -> niter    -> ", _m("niter"))
    println("dmft -> mixer    -> ", _m("mixer"))
    println("dmft -> dcount   -> ", _m("dcount"))
    println("dmft -> cc       -> ", _m("cc"))
    println("dmft -> ec       -> ", _m("ec"))
    println("dmft -> fc       -> ", _m("fc"))
    println("dmft -> lcharge  -> ", _m("lcharge"))
    println("dmft -> lenergy  -> ", _m("lenergy"))
    println("dmft -> lforce   -> ", _m("lforce"))
    println()
end

"""
    view_impurity()

Print the configuration parameters to stdout: for impurity
"""
function view_impurity()
    message("impurity parameters")
    println("impurity -> nsite  -> ", _i("nsite"))
    println("impurity -> atoms  -> ", _i("atoms"))
    println("impurity -> equiv  -> ", _i("equiv"))
    println("impurity -> shell  -> ", _i("shell"))
    println("impurity -> ising  -> ", _i("ising"))
    println("impurity -> occup  -> ", _i("occup"))
    println("impurity -> upara  -> ", _i("upara"))
    println("impurity -> jpara  -> ", _i("jpara"))
    println("impurity -> lpara  -> ", _i("lpara"))
    println()
end

"""
    view_solver()

Print the configuration parameters to stdout: for solver
"""
function view_solver()
    message("solver parameters")
    println("solver   -> engine -> ", _s("engine"))
    println("solver   -> params -> ", _s("params"))
    println()
end

"""
    welcome()

Print out the welcome messages to the screen
"""
function welcome()
    printstyled("                                   |\n", color = :green)
    printstyled(
        "========== ========== ===     ===  | A Modern DFT + DMFT Computation Framework\n",
        color = :green,
    )
    printstyled("        // ||         || n     ||  |\n", color = :green)
    printstyled("       //  ||         ||  n    ||  |\n", color = :green)
    printstyled("  //zz//   ||eeeeeeee ||   n   ||  |\n", color = :green)
    printstyled(
        " //        ||         ||    n  ||  | Version: $__version__\n",
        color = :green,
    )
    printstyled(
        "//         ||         ||     n ||  | Release: $__release__\n",
        color = :green,
    )
    printstyled(
        "========== ========== ===     ===  | Powered by the julia programming language\n",
        color = :green,
    )
    printstyled("                                   |\n", color = :green)
    println()
end

"""
    goodbye()

Print the goodbye messages to the screen
"""
function goodbye()
    println("See you later")
end

function sorry()
    error("Sorry, this feature has not been implemented")
end

function message(from::String, msg::String)
    printstyled("[" * from * "]: ", color = :green)
    println(msg)
end

function message(from::String)
    message(from, "")
end

"""
    line_to_array(io::IOStream)

Convert a line (reading from an iostream) to a string array
"""
@inline function line_to_array(io::IOStream)
    split(readline(io), " ", keepempty = false)
end

"""
    line_to_array(str::AbstractString)

Convert a string (AbstractString) to a string array
"""
@inline function line_to_array(str::AbstractString)
    split(str, " ", keepempty = false)
end
