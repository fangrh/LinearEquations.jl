using LinearEquations
using Documenter

DocMeta.setdocmeta!(LinearEquations, :DocTestSetup, :(using LinearEquations); recursive=true)

makedocs(;
    modules=[LinearEquations],
    authors="Ruihuan Fang",
    repo="https://github.com/fangrh/LinearEquations.jl/blob/{commit}{path}#{line}",
    sitename="LinearEquations.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://fangrh.github.io/LinearEquations.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/fangrh/LinearEquations.jl",
    devbranch="main",
)
