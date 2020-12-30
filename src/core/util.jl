#
# project : pansy
# source  : util.jl
# author  : Li Huang (lihuang.dmft@gmail.com)
# status  : unstable
# comment :
#
# last modified: 2020/12/30
#

"""
    @cswitch

Provides C-like switch statement with the ``falling through'' behavior. this
implement is borrowed from the following github repp.:
    https://github.com/Gnimuc/CSyntax.jl

# Examples
```julia
engine = _d("engine")
@cswitch engine begin
    @case "vasp"
        vasp_init(it)
        break

    @default
        sorry()
        break
end
```
"""
macro cswitch(constexpr, body)
    case2label = Dict{Any,Symbol}()
    flow = Expr(:block)
    end_label = gensym("end")
    default_label = end_label

    for arg in body.args
        if Meta.isexpr(arg, :macrocall) && arg.args[1] == Symbol("@case")
            label = gensym("case")
            case2label[arg.args[3]] = label
            labelexpr = Expr(:symboliclabel, label)
            push!(flow.args, labelexpr)
        elseif Meta.isexpr(arg, :macrocall) && arg.args[1] == Symbol("@default")
            default_label = gensym("default")
            labelexpr = Expr(:symboliclabel, default_label)
            push!(flow.args, labelexpr)
        elseif arg == Expr(:break)
            labelexpr = Expr(:symbolicgoto, end_label)
            push!(flow.args, labelexpr)
        else
            push!(flow.args, arg)
        end
    end
    push!(flow.args, Expr(:symboliclabel, end_label))

    jumptable = Expr(:block)
    for (case, label) in case2label
        condition = Expr(:call, :(==), constexpr, case)
        push!(jumptable.args, Expr(:if, condition, Expr(:symbolicgoto, label)))
    end
    push!(jumptable.args[end].args, Expr(:symbolicgoto, default_label))

    return esc(Expr(:block, jumptable, flow))
end

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
        error("Please specify the configuration file")
    else
        ARGS[1]
    end
end

"""
    query_inps()

Check whether the essential input files exist
"""
function query_inps()
    # if the dft engine is vasp, we have to ensure that the required input
    # files (POSCAR and POTCAR) are present
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
    # we have to setup environment variable ZEN_HOME
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
    # we have to setup environment variable VASP_HOME
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
    list_case()

Print the configuration parameters to stdout: for PCASE dict
"""
function list_case()
    println("< Parameters: case >")
    println("  case     -> ", _c("case"))
    println()
end

"""
    list_dft()

Print the configuration parameters to stdout: for PDFT dict
"""
function list_dft()
    println("< Parameters: dft engine >")
    println("  engine   -> ", _d("engine"))
    println("  smear    -> ", _d("smear"))
    println("  kmesh    -> ", _d("kmesh"))
    println("  magmom   -> ", _d("magmom"))
    println("  lsymm    -> ", _d("lsymm"))
    println("  lspins   -> ", _d("lspins"))
    println("  lspinorb -> ", _d("lspinorb"))
    println("  loptim   -> ", _d("loptim"))
    println("  window   -> ", _d("window"))
    println("  lproj    -> ", _d("lproj"))
    println("  sproj    -> ", _d("sproj"))
    println()
end

"""
    list_dmft()

Print the configuration parameters to stdout: for PDMFT dict
"""
function list_dmft()
    println("< Parameters: dmft engine >")
    println("  mode     -> ", _m("mode"))
    println("  axis     -> ", _m("axis"))
    println("  beta     -> ", _m("beta"))
    println("  niter    -> ", _m("niter"))
    println("  mixer    -> ", _m("mixer"))
    println("  dcount   -> ", _m("dcount"))
    println("  cc       -> ", _m("cc"))
    println("  ec       -> ", _m("ec"))
    println("  fc       -> ", _m("fc"))
    println("  lcharge  -> ", _m("lcharge"))
    println("  lenergy  -> ", _m("lenergy"))
    println("  lforce   -> ", _m("lforce"))
    println()
end

"""
    list_impurity()

Print the configuration parameters to stdout: for PIMP dict
"""
function list_impurity()
    println("< Parameters: quantum impurity atoms >")
    println("  nsite    -> ", _i("nsite"))
    println("  atoms    -> ", _i("atoms"))
    println("  equiv    -> ", _i("equiv"))
    println("  shell    -> ", _i("shell"))
    println("  ising    -> ", _i("ising"))
    println("  occup    -> ", _i("occup"))
    println("  upara    -> ", _i("upara"))
    println("  jpara    -> ", _i("jpara"))
    println("  lpara    -> ", _i("lpara"))
    println()
end

"""
    list_solver()

Print the configuration parameters to stdout: for PSOLVER dict
"""
function list_solver()
    println("< Parameters: quantum impurity solvers >")
    println("  engine   -> ", _s("engine"))
    println("  params   -> ", _s("params"))
    println()
end

"""
    welcome()

Print out the welcome messages to the screen
"""
function welcome()
    printstyled("                                   |\n", color = :green)
    printstyled("========== ========== ===     ===  | ", color = :green)
    printstyled("A Modern DFT + DMFT Computation Framework\n", color = :magenta)
    printstyled("        //            || n     ||  |\n", color = :green)
    printstyled("       //             ||  n    ||  |\n", color = :green)
    printstyled("  zzzzzz    eeeeeeee  ||   n   ||  |\n", color = :green)
    printstyled(" //                   ||    n  ||  | ", color = :green)
    printstyled("Version: $__version__\n", color = :magenta)
    printstyled("//                    ||     n ||  | ", color = :green)
    printstyled("Release: $__release__\n", color = :magenta)
    printstyled("========== ========== ===     ===  | ", color = :green)
    printstyled("Powered by the julia programming language\n", color = :magenta)
    printstyled("                                   |\n", color = :green)
    println()
end

"""
    overview()

Print out the overview of zen to the screen
"""
function overview()
    println("Starting time: ", Dates.format(now(), "yyyy-mm-dd / HH:MM:SS"))
    println("Parallel execution: using ", nprocs(), nprocs() == 1 ? " processor" : " processors")
    println("Current working directory: ", pwd())
    println("Job description file: ", query_args())
    println()
end

"""
    goodbye()

Print the goodbye messages to the screen
"""
function goodbye()
    println("See you later")
end

"""
    sorry()

Print an error message to the screen
"""
function sorry()
    error("Sorry, this feature has not been implemented")
end

"""
    message(from::String, msg::String)

Print an standard zen message to the screen
"""
function message(from::String, msg::String)
    printstyled(from * " > ", color = :green)
    printstyled(msg, color = :magenta)
    println()
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
