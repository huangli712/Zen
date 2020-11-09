module Zen

include("util.jl")
export check_version
export check_home
export welcome
export goodbye

include("parser.jl")
export parse_config

end
