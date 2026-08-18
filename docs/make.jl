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
        "Model Parameters" => "parameters.md",
        "Physics Notes" => "physics_notes.md",
        "Variable Reference" => "variables.md",
        "API Reference" => "api.md",
        "Model Comparisons" => [
            "Community Firn Model" => "cfm_comparison.md",
            "IMAU-FDM" => "imau_fdm_comparison.md",
        ],
        "Internals" => [
            "Overview" => "internals.md",
            "Physics" => "internals_physics.md",
            "Grid and Column" => "internals_grid.md",
            "Support" => "internals_support.md",
        ],
    ],
)

# DocumenterVitepress.deploydocs, *not* Documenter.deploydocs: since DV 0.2 the two are
# incompatible, because one build can emit several version folders (v1.2.3, v1.2, v1,
# stable) that VitePress renders into numbered subdirectories of `target`. Documenter's
# own deploydocs uploads `target` verbatim, which publishes those intermediates
# (`build/1/`, `build/.documenter/`) instead of the site itself, and every page 404s.
DocumenterVitepress.deploydocs(;
    repo="github.com/alex-s-gardner/GEMB.jl",
    target=joinpath(@__DIR__, "build"),
    branch="gh-pages",
    devbranch="main",
    push_preview=true,
)
