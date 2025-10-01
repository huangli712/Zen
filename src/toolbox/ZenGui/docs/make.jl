haskey(ENV,"ZEN_GUI") && pushfirst!(LOAD_PATH, ENV["ZEN_GUI"])

using Documenter
using ZenGui

makedocs(
    sitename = "ZenGui: The User Guide",
    clean = true,
    authors = "Li Huang <huangli@caep.cn> and contributors",
    format = Documenter.HTML(
        prettyurls = false,
        ansicolor = true,
        repolink = "https://github.com/huangli712/ZenGui",
        size_threshold = 409600, # 400kb
        assets = ["assets/zengui.css"],
        collapselevel = 1,
    ),
    #format = Documenter.LaTeX(platform = "none"),
    remotes = nothing,
    modules = [ZenGui],
    pages = [
        "Welcome" => "index.md",
        "Introduction" => "intro.md",
        "Installation" => "install.md",
        "Usage" => "usage.md",
        "Library" => "library.md",
    ],
)
