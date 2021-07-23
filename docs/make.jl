using QuantumAlgebra
using Documenter

DocMeta.setdocmeta!(QuantumAlgebra, :DocTestSetup, :(using QuantumAlgebra); recursive=true)

makedocs(;
    modules=[QuantumAlgebra],
    authors="Johannes Feist <johannes.feist@gmail.com> and contributors",
    repo="https://github.com/jfeist/QuantumAlgebra.jl/blob/{commit}{path}#{line}",
    sitename="QuantumAlgebra.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jfeist.github.io/QuantumAlgebra.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Interface" => "interface.md",
    ],
)

deploydocs(;
    repo="github.com/jfeist/QuantumAlgebra.jl",
    devbranch = "main",
)
