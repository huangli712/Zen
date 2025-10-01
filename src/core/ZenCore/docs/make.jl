haskey(ENV,"ZEN_CORE") && pushfirst!(LOAD_PATH, ENV["ZEN_CORE"])

using Documenter
using ZenCore

makedocs(
    sitename = "ZenCore: The User Guide",
    clean = true,
    authors = "Li Huang <huangli@caep.cn> and contributors",
    format = Documenter.HTML(
        prettyurls = false,
        ansicolor = true,
        repolink = "https://github.com/huangli712/ZenCore",
        size_threshold = 409600, # 400kb
        assets = ["assets/zencore.css"],
        collapselevel = 1,
    ),
    #format = Documenter.LaTeX(platform = "none"),
    remotes = nothing,
    modules = [ZenCore],
    pages = [
        "Welcome" => "index.md",
        "Introduction" => "intro.md",
        "Installation" => "install.md",
        "Inside The Library" => Any[
            "Outline" => "library/outline.md",
            "The ZenCore Module" => "library/zencore.md",
            "Global Constants" => "library/global.md",
            "Auxiliary Utilities" => "library/util.md",
            "Analytical Tetrahedron Method" => "library/tetra.md",
            "Fundamental Types" => "library/types.md",
            "Configuration" => "library/config.md",
            "Basic Workflow" => "library/base.md",
            "Vienna Ab Initio Simulation Package" => "library/vasp.md",
            "Quantum ESPRESSO" => "library/qe.md",
            "Projected Local Orbitals" => "library/plo.md",
            "Wannier Functions" => "library/wannier.md",
            "Intermediate Representation" => "library/ir.md",
            "Dynamical Mean-Field Theory" => "library/dmft.md",
            "Quantum Impurity Solvers" => "library/solver.md",
            "Self-Energy Functions" => "library/sigma.md",
            "Mixing Algorithms" => "library/mixing.md",
        ],
    ],
)
