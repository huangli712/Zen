using TOML

function parse_zen_config(f::AbstractString)
    dict = TOML.parsefile(f)
    @show dict
end
