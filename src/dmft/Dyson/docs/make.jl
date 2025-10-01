using Documenter

makedocs(
    sitename="Dyson: The User Guide",
    clean = true,
    authors = "Li Huang <huangli@caep.cn> and contributors",
    format = Documenter.HTML(
        prettyurls = false,
        ansicolor = true,
        repolink = "https://github.com/huangli712/Dyson",
        size_threshold = 409600, # 400kb
        assets = ["assets/dyson.css"],
        collapselevel = 1,
    ),
    #format = Documenter.LaTeX(platform = "none"),
    remotes = nothing,
    modules = Module[],
    pages = [
        "Welcome" => "index.md",
        "Introduction" => "intro.md",
        "Installation" => "install.md",
        "Usage" => "usage.md",
        "Input Files" => "input.md",
        "Output Files" => "output.md",
        "Parameters" => "param.md", 
    ],
)
