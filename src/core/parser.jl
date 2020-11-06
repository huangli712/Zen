using TOML

println("in parser.jl")

dict = TOML.parsefile("case.toml")
@show typeof(dict)
@show dict
