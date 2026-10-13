"""
    QuantumAlgebra

Quantum operator algebra: define operators (bosons, fermions, two-level systems, qudits, spins and
generators of so(n), see [`@boson`](@ref) etc.) and parameters ([`@qparam`](@ref)), compute normal
forms ([`normal_form`](@ref)), commutators ([`comm`](@ref)), vacuum expectation values
([`vacExpVal`](@ref)), Heisenberg equations of motion ([`heisenberg_eom`](@ref),
[`heisenberg_eom_system`](@ref)) and cumulant expansions ([`expval_as_corrs`](@ref)), and convert
the results to LaTeX ([`latex`](@ref)) or Julia code ([`julia_expression`](@ref)).
"""
module QuantumAlgebra

include("tools.jl")
include("operator_defs.jl")
include("index_handling.jl")
include("operator_rules.jl")
include("operator_baseops.jl")
include("sum_indices.jl")
include("alias.jl")
include("output.jl")
include("correlations.jl")
include("vacuum_expvals.jl")
include("convert_to_expression.jl")
include("eqsofmotion.jl")
include("eqsys.jl")

include("precompile.jl")

function __init__()
    auto_normal = @load_preference("auto_normal_form", false)
    if haskey(ENV,"QUANTUMALGEBRA_AUTO_NORMAL_FORM")
        val = parse(Bool,ENV["QUANTUMALGEBRA_AUTO_NORMAL_FORM"])
        @warn "Environment variable QUANTUMALGEBRA_AUTO_NORMAL_FORM is set to $val. This is deprecated, please instead set it through Preferences.jl by calling `QuantumAlgebra.auto_normal_form($val; set_preference=true) or `Preferences.set_preferences!(QuantumAlgebra,\"auto_normal_form\"=>$val)."
        if !@has_preference("auto_normal_form")
            auto_normal = val
        end
    end
    auto_normal_form(auto_normal)

    use_σpm(@load_preference("use_σpm", false))
    use_bracket_display(@load_preference("bracket_display", false))
end

end # module
