module QuantumAlgebra

include("operator_defs.jl")
include("operator_iteration.jl")
include("operator_baseops.jl")
include("index_handling.jl")
include("output.jl")
include("correlations.jl")
include("vacuum_expvals.jl")

include("precompile.jl")
_precompile_()

end # module
