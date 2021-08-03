using Documenter

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
        "Introduction" => "intro/readme.md",
        "Getting started" => "start/readme.md",
        "Tutorials" => "tutor/readme.md",
        "Guide" => "guide/readme.md",
        "Internals" => "internals/readme.md",
        "Theory" => "theory/readme.md",
    ],
)
