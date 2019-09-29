using Documenter, QuantumAlgebra

makedocs(;
    modules=[QuantumAlgebra],
    format=Documenter.HTML(assets=String[]),
    pages=[
        "Home" => "index.md",
        "Interface" => "interface.md",
    ],
    repo="https://github.com/jfeist/QuantumAlgebra.jl/blob/{commit}{path}#L{line}",
    sitename="QuantumAlgebra.jl",
    authors="Johannes Feist"
)

deploydocs(;
    repo="github.com/jfeist/QuantumAlgebra.jl",
)
