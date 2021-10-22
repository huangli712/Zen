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
    modules = [ZenCore],
    pages = [
        "Home" => "index.md",
        "Introduction" => Any[],
        "Getting started" => Any[],
        "Tutorials" => Any[],
        "Guide" => Any[],
        "Internals" => Any[
            "README" => "internals/README.md",
            "ZenCore APIs" => Any[
                "ZenCore" => "internals/apis/zencore.md",
                "Global" => "internals/apis/global.md",
                "Util" => "internals/apis/util.md",
                "Tetra" => "internals/apis/tetra.md",
                "Types" => "internals/apis/types.md",
                "Config" => "internals/apis/config.md",
                "Base" => "internals/apis/base.md",
                "VASP" => "internals/apis/vasp.md",
                "QE" => "internals/apis/qe.md",
                "PLO" => "internals/apis/plo.md",
                "WANNIER" => "internals/apis/wannier.md",
                "IR" => "internals/apis/ir.md",
                "DMFT" => "internals/apis/dmft.md",
                "Solver" => "internals/apis/solver.md",
                "Sigma" => "internals/apis/sigma.md",
                "Mixer" => "internals/apis/mixer.md",
            ],
        ],
        "Theory" => Any[],
    ],
)
