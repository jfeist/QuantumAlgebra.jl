module QuantumAlgebra

# helper function used in some definitons, copied from ColorTypes.jl
Base.@pure basetype(T::Type) = Base.typename(T).wrapper
Base.@pure basetype(A) = basetype(typeof(A))

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
