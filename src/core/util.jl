#
# project : pansy
# source  : util.jl
# author  : Li Huang (lihuang.dmft@gmail.com)
# status  : unstable
# comment :
#
# last modified: 2021/01/12
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
            error("VASP_HOME is undefined")
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
    #
    # remarks:
    #
    # _d("sproj") is actually an Array{String,1}. it would be quite low
    # efficiency if we print it directly. so we convert it into a string
    # by using the join() function at first.
    #
    # _d("ewidth") is actually a real number. it would be quite low
    # efficiency if we print it directly. so we convert it into a string
    # by using the string() function at first.
    #
    # see config.jl/str_d() for more details
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
    list_dmft()

Print the configuration parameters to stdout: for PDMFT dict
"""
function list_dmft()
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
    list_impurity()

Print the configuration parameters to stdout: for PIMP dict
"""
function list_impurity()
    # see comments in list_dft()
    println("< Parameters: quantum impurity atoms >")
    println("  nsite    -> ", _i("nsite"))
    println("  atoms    -> ", join(_i("atoms"), "; "))
    println("  equiv    -> ", join(_i("equiv"), "; "))
    println("  shell    -> ", join(_i("shell"), "; "))
    println("  ising    -> ", join(_i("ising"), "; "))
    println("  occup    -> ", join(_i("occup"), "; "))
    println("  upara    -> ", join(_i("upara"), "; "))
    println("  jpara    -> ", join(_i("jpara"), "; "))
    println("  lpara    -> ", join(_i("lpara"), "; "))
    println()
end

"""
    list_solver()

Print the configuration parameters to stdout: for PSOLVER dict
"""
function list_solver()
    # see comments in list_solver()
    println("< Parameters: quantum impurity solvers >")
    println("  engine   -> ", _s("engine"))
    println("  params   -> ", join(_s("params"), "; "))
    println()
end

"""
    welcome()

Print out the welcome messages to the screen
"""
function welcome()
    printstyled("                                        |\n", color = :green)
    printstyled("ZZZZZZZZZZZZ EEEEEEEEEEEE NNNNNNNNNNNN  | ", color = :green)
    printstyled("A Modern DFT + DMFT Computation Framework\n", color = :magenta)
    printstyled("          Z               N          N  |\n", color = :green)
    printstyled("         Z                N          N  |\n", color = :green)
    printstyled("   ZZZZZZ    EEEEEEEEEEEE N          N  |\n", color = :green)
    printstyled("  Z                       N          N  | ", color = :green)
    printstyled("Version: $__VERSION__\n", color = :magenta)
    printstyled(" Z                        N          N  | ", color = :green)
    printstyled("Release: $__RELEASE__\n", color = :magenta)
    printstyled("ZZZZZZZZZZZZ EEEEEEEEEEEE N          N  | ", color = :green)
    printstyled("Powered by the julia programming language\n", color = :magenta)
    printstyled("                                        |\n", color = :green)
    println()
end

"""
    overview()

Print out the overview of zen to the screen
"""
function overview()
    # build strings
    str1 = nprocs() === 1 ? " processor " : " processors "
    str2 = "(myid = $(myid()))"

    # write the information
    println("Starting time: ", Dates.format(now(), "yyyy-mm-dd / HH:MM:SS"))
    println("Parallel execution: using ", nprocs(), str1, str2)
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
