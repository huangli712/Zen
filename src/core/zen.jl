module Zen

const version = v"0.0.1"
const authors = "Li Huang" 

include("util.jl")
export check_version
export check_home
export welcome
export goodbye

include("parser.jl")
export parse_config

end
