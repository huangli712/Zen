
include("util.jl")
include("parser.jl")

using .ZenUtil
using .ZenParser


check_version()
parse_zen_config("case.toml")
