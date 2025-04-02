#
# Project : Pansy
# Source  : util.jl
# Author  : Li Huang (huangli@caep.cn)
# Status  : Unstable
#
# Last modified: 2025/03/26
#

#=
### *Basic Macros*
=#

"""
    @cswitch(constexpr, body)

Provides a C-like switch statement with the *falling through* behavior.
This implementation was borrowed from the following github repository:

* https://github.com/Gnimuc/CSyntax.jl

### Examples
```julia
engine = get_d("engine")
@cswitch engine begin
    @case "vasp"
        just_do_it()
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
    @time_call(ex)

Evaluate a function call (`ex`), and then print the elapsed time (number
of seconds) it took to execute.

This macro is a variation of the standard `@elapsed` macro.
"""
macro time_call(ex)
    quote
        while false; end
        local t₀ = time_ns()
        $(esc(ex))
        δt = (time_ns() - t₀) / 1e9
        println("Report: Total elapsed time $(δt) s\n")
        flush(stdout)
    end
end

"""
    @pcs(x...)

Try to print colorful strings. Here `x` is a combination of strings and
colors. Its format likes `string1 color1 string2 color2 (repeat)`. For
the supported colors, please check the global dict `COLORS`.

### Examples
```julia-repl
julia> @pcs "Hello world!" blue
julia> @pcs "Hello " red "world!" green
```

See also: [`COLORS`](@ref), [`welcome`](@ref).
"""
macro pcs(x...)
    ex = quote
        # The `args` is actually a Tuple
        args = $x

        # We have to make sure the strings and colors are paired.
        @assert iseven(length(args))

        for i = 1:2:length(args)
            # Construct and check string
            # Sometimes args[i] contains interpolated variables, its
            # type is `Expr`. At this time, we have to evaluate this
            # `Expr` at first to convert it to a format `String`.
            str   = eval(args[i])
            @assert str isa AbstractString
            #
            # Construct and check color
            color = args[i+1]
            @assert color isa Symbol

            # Generate expression
            print(eval(color)(str))
        end
    end

    return :( $(esc(ex)) )
end

#=
### *Query Runtime Environment*
=#

"""
    require()

Check the version of julia runtime environment. It should be higher
than v1.6.x. One of the most important philosophies of the `ZenCore`
package is minimizing the dependence on the third-party libraries as
far as possible. Note that the `ZenCore` package relys on the `TOML`
package to parse the *.toml file. Only in v1.6.0 and higher versions,
julia includes the `TOML` package in its standard library.
"""
function require()
    if VERSION < v"1.6-"
        error("Please upgrade your julia to v1.6.0 or higher")
    end
end

"""
    setup_args(x::Vararg{String})

Setup `ARGS` manually. This function is used only in `REPL` environment.
We can use this function to update `ARGS`, so that the `query_args()`
and the other related functions can work correctly.

### Examples
```julia-repl
julia> setup_args("SrVO3.toml")
1-element Array{String,1}:
 "SrVO3.toml"
```

See also: [`query_args`](@ref).
"""
function setup_args(x::Vararg{String})
    # Make sure it is the REPL
    @assert isinteractive()

    # Clean `ARGS`
    empty!(ARGS)

    # Convert the arguments to an array of strings
    X = collect(x)

    # Push them into `ARGS` one by one
    for i in eachindex(X)
        push!(ARGS, X[i])
    end

    # Return `ARGS`, only for debug.
    ARGS
end

"""
    query_args()

Check whether the configuration file (`case.toml`) is provided.

See also: [`setup_args`](@ref).
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
    query_case()

Return `case`, in other words, the job's name.

See also: [`query_stop`](@ref).
"""
function query_case()
    basename( splitext(query_args())[1] )
end

#=
*Remarks 1* :

For `vasp`, the essential input files include:

* POSCAR
* POTCAR

As for the `INCAR` and `KPOINTS` files, they should be generated
automatically. Do not worry about them.

*Remarks 2* :

For `quantum espresso` (`pwscf`), the essential input files include:

* QE.INP
* Pseudopotential files

As for the other files, they should be generated automatically. Do
not worry about them.
=#

"""
    query_inps(::NULLEngine)
    query_inps(::VASPEngine)
    query_inps(::QEEngine)
    query_inps(::WANNIEREngine)

Check whether the essential input files exist. It acts as a dispatcher.
This function is designed for the DFT engine only. The input files for
the DMFT engine, quantum impurity solver, and Kohn-Sham adaptor will be
generated automatically by default. The `ZenCore` package will take care
of them. Do not worry about that.

See also: [`query_case`](@ref).
"""
function query_inps(::NULLEngine)
    sorry()
end
#
function query_inps(::VASPEngine)
    # For vasp code.
    if !isfile("POSCAR") || !isfile("POTCAR")
        error("Please provide both POSCAR and POTCAR files")
    end
end
#
function query_inps(::QEEngine)
    # For quantum espresso code.
    if !isfile("QE.INP")
        error("Please provide the QE.INP file")
    end
end
#
function query_inps(::WANNIEREngine)
    # For wannier90 code.
    sorry()
end

"""
    query_stop()

Query whether the `case.stop` file exists.

See also: [`query_case`](@ref).
"""
function query_stop()
    isfile(query_case() * ".stop")
end

"""
    query_test()

Query whether the `case.test` file exists.

See also: [`query_case`](@ref).
"""
function query_test()
    isfile(query_case() * ".test")
end

#=
*Remarks* :

In the `ZenCore` package, the following environment variables matter:

* ZEN_HOME
* ZEN_CORE
* ZEN_DMFT
* ZEN_SOLVER
* VASP_HOME (If we are using `vasp` as our DFT backend)
* QE_HOME (If we are using `quantum espresso` as our DFT backend)
* WAN90_HOME (If we are using `wannier90` to generate projectors)

Please setup them in your `.bashrc` (Lniux) or `.profile` (macOS) files.
=#

"""
    query_home()

Query the home directory of `Zen`. Actually, the `ZEN_HOME` means the
directory that the `Zen` Framework is installed.

See also: [`query_core`](@ref).
"""
function query_home()
    # We have to setup the environment variable ZEN_HOME
    if haskey(ENV, "ZEN_HOME")
        ENV["ZEN_HOME"]
    else
        error("ZEN_HOME is undefined")
    end
end

"""
    query_core()

Query the src/core directory of `Zen`. Actually, the `ZEN_CORE` denotes
the directory that contains the ZenCore.jl file. Be careful, `ZEN_CORE`
must be included in `LOAD_PATH`.

See also: [`query_home`](@ref).
"""
function query_core()
    # We have to setup the environment variable ZEN_CORE
    if haskey(ENV, "ZEN_CORE")
        ENV["ZEN_CORE"]
    else
        error("ZEN_CORE is undefined")
    end
end

"""
    query_dft(ae::AbstractEngine)

Query the home directory of the chosen DFT backend. It supports vasp,
quantum espresso, and wannier90 by now.

See also: [`query_dmft`](@ref), [`query_solver`](@ref).
"""
function query_dft(ae::AbstractEngine)
    # Build valid name for environment variable
    var = uppercase(nameof(ae)) * "_HOME"

    # Query the environment variable
    if haskey(ENV, var)
        return ENV[var]
    else
        error("$var is undefined")
    end
end

"""
    query_dmft()

Query the home directory of the DMFT engine.

See also: [`query_dft`](@ref), [`query_solver`](@ref).
"""
function query_dmft()
    # We have to setup the environment variable ZEN_DMFT
    if haskey(ENV, "ZEN_DMFT")
        ENV["ZEN_DMFT"]
    # For develop stage only
    else
        joinpath(query_home(), "src/dmft")
    end
end

"""
    query_solver(as::AbstractSolver)

Query the home directories of various quantum impurity solvers. Now it
supports CTHYB₁, CTHYB₂, HIA, and NORG.

See also: [`query_dft`](@ref), [`query_dmft`](@ref).
"""
function query_solver(as::AbstractSolver)
    # We have to setup the environment variable ZEN_SOLVER
    if haskey(ENV, "ZEN_SOLVER")
        ENV["ZEN_SOLVER"]
    # For develop stage only
    else
        joinpath(query_home(), "src/solver/" * nameof(as))
    end
end

"""
    is_vasp()

Test whether the DFT backend is the `vasp` code.

See also: [`is_qe`](@ref).
"""
function is_vasp()
    get_d("engine") == "vasp"
end

"""
    is_qe()

Test whether the DFT backend is the `quantum espresso` (`pwscf`) code.

See also: [`is_vasp`](@ref).
"""
function is_qe()
    get_d("engine") == "qe"
end

"""
    is_plo()

Test whether the projector is the projected local orbitals.

See also: [`is_wannier`](@ref).
"""
function is_plo()
    get_d("projtype") == "plo"
end

"""
    is_wannier()

Test whether the projector is the wannier functions.

See also: [`is_plo`](@ref).
"""
function is_wannier()
    get_d("projtype") == "wannier"
end

#=
### *Colorful Outputs*
=#

"""
    welcome()

Print out the welcome messages to the screen.
"""
function welcome()
    @pcs "                                       |\n" green
    @pcs "ZZZZZZZZZZZZ EEEEEEEEEEEE NNNNNNNNNNNN | "  green "A Modern DFT + DMFT Computation Framework\n" magenta
    @pcs "          Z               N          N |\n" green
    @pcs "         Z                N          N |\n" green
    @pcs "   ZZZZZZ    EEEEEEEEEEEE N          N | "  green "Package: $__LIBNAME__\n" magenta
    @pcs "  Z                       N          N | "  green "Version: $__VERSION__\n" magenta
    @pcs " Z                        N          N | "  green "Release: $__RELEASE__\n" magenta
    @pcs "ZZZZZZZZZZZZ EEEEEEEEEEEE N          N | "  green "Powered by the julia programming language\n" magenta
    @pcs "                                       |\n" green
    println()
    #
    flush(stdout)
end

"""
    overview()

Print out the overview of Zen to the screen.
"""
function overview()
    # Build strings
    str1 = nprocs() == 1 ? " processor " : " processors "
    str2 = "(myid = $(myid()))"

    # Write the information
    prompt("Overview")
    println("Time : ", Dates.format(now(), "yyyy-mm-dd / HH:MM:SS"))
    println("Para : Using ", nprocs(), str1, str2)
    println("Dirs : ", pwd())
    println("Task : ", query_args())
    println()
    #
    flush(stdout)
end

"""
    goodbye()

Print the goodbye messages to the screen.
"""
function goodbye()
    println(  red("╔═╗┌─┐┌┐┌"), magenta("╔═╗┌─┐┬─┐┌─┐"))
    println(green("╔═╝├┤ │││"), magenta("║  │ │├┬┘├┤ "))
    println( blue("╚═╝└─┘┘└┘"), magenta("╚═╝└─┘┴└─└─┘"))
    #
    flush(stdout)
end

"""
    sorry()

Print an error message to the screen.
"""
function sorry()
    error("Sorry, this feature has not been implemented")
end

"""
    prompt(msg::String)

Print a stylized Zen message to the screen.
"""
function prompt(msg::String)
    print(green("ZEN > "))
    print(magenta(msg))
    println()
    #
    flush(stdout)
end

"""
    prompt(msg1::String, msg2::String)

Print a stylized Zen message to the screen.
"""
function prompt(msg1::String, msg2::String)
    print(blue("Task -> "))
    print(light_red(msg1))
    print(magenta(msg2))
    println()
    #
    flush(stdout)
end

"""
    prompt(io::IOStream, msg::String)

Print a stylized Zen message to the given IOStream. This function is used
to log the key events during DFT + DMFT iterations.
"""
function prompt(io::IOStream, msg::String)
    date = Dates.format(now(), "yyyy-mm-dd / HH:MM:SS")
    println(io, "$date  $msg")
    #
    flush(io)
end

#=
### *I/O Operations*
=#

"""
    line_to_array(io::IOStream)

Convert a line (reading from an IOStream) to a string array.
"""
@inline function line_to_array(io::IOStream)
    split(readline(io), " ", keepempty = false)
end

"""
    line_to_array(str::AbstractString)

Convert a string (AbstractString) to a string array.
"""
@inline function line_to_array(str::AbstractString)
    split(str, " ", keepempty = false)
end

"""
    line_to_cmplx(io::IOStream)

Convert a line (reading from an IOStream) to a cmplx number. It is used
to parse the `LOCPROJ` file only.

See also: [`vaspio_projs`](@ref).
"""
@inline function line_to_cmplx(io::IOStream)
    str = readline(io)
    _re = str[10:28]
    _im = str[30:end]
    return parse(F64, _re) + parse(F64, _im) * im
end

#=
### *Mathematical Functions*
=#

#=
*Remarks 1* :

The definition of Gauss error function is as follows:

```math
erf(x) = \frac{2}{\sqrt{\pi}}\int^{x}_{0} e^{-\eta^2} d\eta.
```

*Remarks 2* :

We call the `erf()` function defined in the mathematical library `libm`
or `openlibm` directly, instead of implementing it again by ourselves.
For more details about `libm` and `openlibm`, please visit the following
websites:

* https://openlibm.org
* https://github.com/JuliaMath/openlibm
* https://sourceware.org/newlib/libm.html

*Remarks 3* :

This below implementation is taken from the `SpecialFunctions.jl`. See:

* https://github.com/JuliaMath/SpecialFunctions.jl

*Remarks 4* :

`Base.Math.libm` is actually a string. It denotes `libopenlibm`.
=#

"""
    erf(x::F64)

Calculate the Gauss error function.

See also: [`gauss_weight`](@ref).
"""
function erf(x::F64)
    ccall(("erf", libm), F64, (F64,), x)
end

#=
### *Utility Functions*
=#

"""
    subscript(num::I64)

Convert a number (it must be in [0,9]) to subscript.
"""
function subscript(num::I64)
    @assert 0 ≤ num ≤ 9
    SUB = ["\u2080" "\u2081" "\u2082" "\u2083" "\u2084" "\u2085" "\u2086" "\u2087" "\u2088" "\u2089"]
    return SUB[num + 1]
end

"""
    str_to_struct(str::AbstractString, postfix::AbstractString)

Convert a string (`str`) to an instance of struct. Here `postfix` could
be `Engine`, `Solver`, `Adaptor`, `Mode`, and `Mixer`.

### Examples
```julia-repl
julia> str_to_struct("vasp", "Engine")
ZenCore.VASPEngine()

julia> str_to_struct("hia", "Solver")
ZenCore.HIASolver()
```
"""
@inline function str_to_struct(str::AbstractString, postfix::AbstractString)
    # Assemble the name of the struct
    fstr = replace(uppercase(str), "_" => "") * postfix

    # Perhaps `fstr` contains some numbers, such as CTHYB1Solver.
    # Replace number with its subscript's form (CTHYB₁Solver)
    for i in 0:9
        if contains(fstr, string(i))
            fstr = replace(fstr, string(i) => subscript(i))
        end
    end

    # Generate an instance
    sym = Symbol(fstr)
    @eval machine = ($sym)()

    # Return the desired value
    return machine
end

#=
### *Color Tools*
=#

#=
*Remarks* :

The purpose of the following codes is to provide some convenient tools
to output colorful and stylized texts in the terminal. Actually, these
codes are inspried by this repository:

* https://github.com/Aerlinger/AnsiColor.jl

For more information about the ANSI color escape sequences, please check
the following websites further:

* https://stackoverflow.com/questions/4842424/
* https://en.wikipedia.org/wiki/ANSI_escape_code

Note that the macro `@pcs` and functions `prompt()` rely on these codes.
=#

"""
    COLORS

A global dict, which is used to specify the system colors.
"""
const COLORS = Dict{String,I64}(
    "black"          => 0,
    "red"            => 1,
    "green"          => 2,
    "yellow"         => 3,
    "blue"           => 4,
    "magenta"        => 5,
    "cyan"           => 6,
    "white"          => 7,
    "default"        => 9,
    "light_black"    => 60,
    "light_red"      => 61,
    "light_green"    => 62,
    "light_yellow"   => 63,
    "light_blue"     => 64,
    "light_magenta"  => 65,
    "light_cyan"     => 66,
    "light_white"    => 67
)

"""
    MODES

A global dict, which is used to specify the mode for output characters.
"""
const MODES = Dict{String,I64}(
    "default"        => 0,
    "bold"           => 1,
    "underline"      => 4,
    "blink"          => 5,
    "swap"           => 7,
    "hide"           => 8
)

"""
    colorize(c::String, s::String; bg::String = "default", m::String="default")

Return some escape sequences, which will be displayed as colorized texts
in the terminal.
"""
function colorize(c::String, s::String; bg::String = "default", m::String="default")
    C_OFFSET = 30
    B_OFFSET = 40
    "\033[$(MODES[m]);$(C_OFFSET + COLORS[c]);$(B_OFFSET + COLORS[bg])m$(s)\033[0m"
end

"""
    colorize(c::String, s::String; bg::String = "default", m::String="default")

Return some escape sequences, which will be displayed as colorized texts
in the terminal.
"""
function colorize(c::Symbol, s::String; bg::String = "default", m::String="default")
    colorize(string(c), s; bg=bg, m=m)
end

#=
*Remarks* :

The following codes will generate and export dynamically some color
functions, including:

```julia
# For standard colors
black(str::String)
red(str::String)
green(str::String)
yellow(str::String)
blue(str::String)
magenta(str::String)
cyan(str::String)
white(str::String)
```

and their light color versions

```julia
# For light colors
light_black(str::String)
light_red(str::String)
light_green(str::String)
light_yellow(str::String)
light_blue(str::String)
light_magenta(str::String)
light_cyan(str::String)
light_white(str::String)
```

These functions provide some shortcuts to create texts decorated by
special escape sequences. These texts will be show as colorized texts
in the terminal.

### Examples
```julia-repl
julia> println(red("hello world!"))
```
=#

export COLORS
export MODES
export colorize

for k in keys(COLORS)
    f = Symbol(k)
    k == "default" && continue
    @eval ($f)(str::String) = colorize(Symbol($f), str)
    @eval export ($f)
end
