module QuantumAlgebra

include("tools.jl")
include("operator_defs.jl")
include("index_handling.jl")
include("operator_baseops.jl")
include("output.jl")
include("correlations.jl")
include("vacuum_expvals.jl")
include("convert_to_expression.jl")

include("precompile.jl")
_precompile_()

function __init__()
    # WARNING: This does not work since using @eval in __init__ means the module
    # cannot be used inside any other module. The error upon precompilation of
    # the other module is:
    # ERROR: LoadError: InitError: Evaluation into the closed module
    # `QuantumAlgebra` breaks incremental compilation because the side effects
    # will not be permanent. This is likely due to some other module mutating
    # `QuantumAlgebra` with `eval` during precompilation - don't do this.

    # default_ops = parse(Bool,get(ENV,"QUANTUMALGEBRA_DEFAULT_OPS","true"))
    # if default_ops
    #     @eval begin
    #         @boson_ops a
    #         @fermion_ops f
    #         @tlspm_ops σ
    #         @tlsxyz_ops σ
    #         export σx, σy, σz, σp, σm
    #         export a, adag, f, fdag
    #     end
    # end

    auto_normal = parse(Bool,get(ENV,"QUANTUMALGEBRA_AUTO_NORMAL_FORM","false"))
    auto_normal_form(auto_normal)
end

end # module
