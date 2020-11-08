module ZenParser

using TOML

export parse_zen_config

function parse_zen_config(f::AbstractString)
    dict = TOML.parsefile(f)
    @show dict
end

end
