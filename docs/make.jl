push!(LOAD_PATH, ENV["ZEN_CORE"])

using Documenter, ZenCore

makedocs(
    sitename = "Zen",
    clean = false,
    authors = "Li Huang",
    format = Documenter.HTML(
        prettyurls = false,
        ansicolor = true,
    ),
    modules = [ZenCore],
    pages = [
        "Home" => "index.md",
        "Introduction" => "intro.md",
        "Getting started" => Any[
            "README" => "start/README.md",
            "Download Zen" => Any[],
            "Compile Zen" => Any[],
            "Configure Zen" => Any[],
        ],
        "Tutorials" => Any[
            "README" => "tutor/README.md",
            "SrVO3" => Any[],
            "FeSe" => Any[],
            "NiO" => Any[],
            "Ce" => Any[],
        ],
        "Guide" => Any[
            "README" => "guide/README.md",
            "Mode" => Any[
                "Interactive mode" => Any[],
                "Batch Mode" => Any[],    
            ],
            "Core applications" => Any[
                "Input parameters" => Any[
                    "PCASE block" => "guide/para/case.md",
                    "PDFT block" => "guide/para/dft.md",
                    "PDMFT block" => "guide/para/dmft.md",
                    "PIMP block" => "guide/para/impurity.md",
                    "PSOLVER block" => "guide/para/solver.md",
                ],    
            ],
            "Auxiliary tools" => Any[],
            "Components" => Any[
                "Density functional theory" => Any[],
                "Dynamical mean-field theory" => Any[],
                "Quantum impurity solver" => Any[],
                "Adaptor" => Any[],    
            ],
            "Files" => Any[
                "Standard output" => Any[],
                "case.cycle" => Any[],
                "case.log" => Any[],
                "case.stop" => Any[],
                "case.test" => Any[],    
            ],
        ],
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
                "Wannier" => "internals/apis/wannier.md",
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
