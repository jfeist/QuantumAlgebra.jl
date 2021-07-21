import Pkg
Pkg.activate("..")

using SnoopCompile

inf_timing = @snoopi tmin=0.0 begin
    using QuantumAlgebra
    include("./snoop_statements.jl")
    include("../test/runtests.jl")
end

pc = SnoopCompile.parcel(inf_timing)

SnoopCompile.write("/tmp/precompile", pc)
