haskey(ENV,"ACFLOW_HOME") && pushfirst!(LOAD_PATH, ENV["ACFLOW_HOME"])
haskey(ENV,"ACGUI_HOME") && pushfirst!(LOAD_PATH, ENV["ACGUI_HOME"])

using Documenter
using Random
using ACGui

makedocs(
    sitename = "ACGui: The User Guide",
    clean = true,
    authors = "Li Huang <huangli@caep.cn> and contributors",
    format = Documenter.HTML(
        prettyurls = false,
        ansicolor = true,
        repolink = "https://github.com/huangli712/ACGui",
        size_threshold = 409600, # 400kb
        assets = ["assets/acgui.css"],
        collapselevel = 1,
    ),
    #format = Documenter.LaTeX(platform = "none"),
    remotes = nothing,
    modules = [ACGui],
    pages = [
        "Welcome" => "index.md",
        "Introduction" => "intro.md",
        "Installation" => "install.md",
        "Usage" => "usage.md",
        "Library" => "library.md",
    ],
)
