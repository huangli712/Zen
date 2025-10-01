haskey(ENV,"ACTEST_HOME") && pushfirst!(LOAD_PATH, ENV["ACTEST_HOME"])

using Documenter
using Random
using ACTest

makedocs(
    sitename = "ACTest: The User Guide",
    clean = true,
    authors = "Li Huang <huangli@caep.cn> and contributors",
    format = Documenter.HTML(
        prettyurls = false,
        ansicolor = true,
        repolink = "https://github.com/huangli712/ACTest",
        size_threshold = 409600, # 400kb
        assets = ["assets/actest.css"],
        collapselevel = 1,
    ),
    #format = Documenter.LaTeX(platform = "none"),
    remotes = nothing,
    modules = [ACTest],
    pages = [
        "Welcome" => "index.md",
        "Introduction" => Any[
            "Motivation" => "intro/motivation.md",
            "Acknowledgements" => "intro/ack.md",
            "Citation" => "intro/cite.md",
        ],
        "Manual" => Any[
            "Main Features" => "man/feature.md",
            "Installation" => "man/install.md",
            "Scripts" => "man/script.md",
            "Inputs" => "man/input.md",
            "Outputs" => "man/output.md",
            "Parameters" => "man/param.md",
            "Built-in Testing Dataset" => "man/act100.md",
            "Interface To Analytic Continuation Toolkits" => "man/interface.md",
            "Graphic User Interface" => "man/gui.md",
        ],
        "Theory" => Any[
            "Grids" => "theory/grid.md",
            "Meshes" => "theory/mesh.md",
            "Peaks" => "theory/peak.md",
            "Kernels" => "theory/kernel.md",
            "Noise" => "theory/noise.md",
        ],
        "Examples" => Any[
            "Generating Spectra And Correlators" => "examples/generate.md",
            "Analytic Continuation Simulations" => "examples/acflow.md",
            "Visualizations" => "examples/plot.md",
        ],
        "Library" => Any[
            "Outline" => "library/outline.md",
            "ACTest" => "library/actest.md",
            "Constants" => "library/global.md",
            "Types" => "library/type.md",
            "Core" => "library/base.md",
            "Peaks" => "library/peak.md",
            "Spectra" => "library/spectrum.md",
            "Standard Dataset" => "library/dataset.md",
            "Grids" => "library/grid.md",
            "Meshes" => "library/mesh.md",
            "Kernels" => "library/kernel.md",
            "Configuration" => "library/config.md",
            "Input And Output" => "library/inout.md",
            "Math" => "library/math.md",
            "Utilities" => "library/util.md",
        ],
    ],
)
