import Pkg
Pkg.activate("..")

using SnoopCompile

inf_timing = @snoopi begin
    using QuantumAlgebra
    include("./snoop_statements.jl")
    include("../test/runtests.jl")
end

pc = SnoopCompile.parcel(inf_timing)

SnoopCompile.write("/tmp/precompile", pc)
