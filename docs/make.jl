using GEMB
using Documenter

DocMeta.setdocmeta!(GEMB, :DocTestSetup, :(using GEMB); recursive=true)

makedocs(;
    modules=[GEMB],
    authors="Alex Gardner <alex.s.gardner@jpl.nasa.gov> and contributors",
    sitename="GEMB.jl",
    format=Documenter.HTML(;
        canonical="https://alex-s-gardner.github.io/GEMB.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Variable Reference" => "variables.md",
        "API Reference" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/alex-s-gardner/GEMB.jl",
    devbranch="main",
)
