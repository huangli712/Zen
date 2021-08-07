push!(LOAD_PATH, ENV["ZEN_CORE"])

using Documenter, ZenCore

makedocs(
    sitename="Zen",
    clean = false,
    authors = "Li Huang",
    format = Documenter.HTML(
        prettyurls = false,
        ansicolor = true,
    ),
    pages = [
        "Home" => "index.md",
#        "Introduction" => "intro/readme.md",
#        "Getting started" => "start/readme.md",
#        "Tutorials" => "tutor/readme.md",
#        "Guide" => "guide/readme.md",
        "Internals" => Any[
            "README" => "internals/README.md",
            "ZenCore APIs" => Any[
                "ZenCore" => "internals/apis/zencore.md",
                "Global" => "internals/apis/global.md",
                "Util" => "internals/apis/util.md",
                "Tetra" => "internals/apis/tetra.md",
            ],
        ],
#        "Theory" => "theory/readme.md",
    ],
    modules = [ZenCore],
)
