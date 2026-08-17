using GEMB
using Documenter
using DocumenterVitepress

DocMeta.setdocmeta!(GEMB, :DocTestSetup, :(using GEMB); recursive=true)

makedocs(;
    modules=[GEMB],
    authors="Alex Gardner <alex.s.gardner@jpl.nasa.gov> and contributors",
    repo="https://github.com/alex-s-gardner/GEMB.jl",
    sitename="GEMB.jl",
    format=DocumenterVitepress.MarkdownVitepress(;
        repo="https://github.com/alex-s-gardner/GEMB.jl",
        devbranch="main",
        devurl="dev",
    ),
    pages=[
        "Home" => "index.md",
        "Model Architecture" => "architecture.md",
        "Variable Reference" => "variables.md",
        "API Reference" => "api.md",
        "Internals" => [
            "Overview" => "internals.md",
            "Physics" => "internals_physics.md",
            "Grid and Column" => "internals_grid.md",
            "Support" => "internals_support.md",
        ],
    ],
)

deploydocs(;
    repo="github.com/alex-s-gardner/GEMB.jl",
    target="build",   # DocumenterVitepress writes the rendered site here
    branch="gh-pages",
    devbranch="main",
    push_preview=true,
)
