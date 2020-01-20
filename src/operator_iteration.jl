import Base.iterate

struct OpProdIter{swallow_∑} A::Operator end
Base.eltype(::OpProdIter) = Operator
proditer(A::Operator,swallow_∑::Bool) = OpProdIter{swallow_∑}(A)
proditer(A::OpSum,swallow_∑::Bool) = throw(ArgumentError("cannot get proditer for OpSum!"))
# if iter.A is an OpProd, default function argument will call more specific function below
iterate(iter::OpProdIter, state::Operator=iter.A) = (state,nothing)
iterate(iter::OpProdIter, state::OpProd) = (@assert !isa(state.A,OpProd); (state.A,state.B))
iterate(iter::OpProdIter{true},  state::OpSumAnalytic) = iterate(iter, state.A)
iterate(iter::OpProdIter{false}, state::OpSumAnalytic) = (state,nothing)
iterate(iter::OpProdIter, state::Nothing) = state

prodtuple(A::Operator) = (A,)
prodtuple(A::OpSum) = throw(ArgumentError("cannot get prodtuple for OpSum!"))
prodtuple(A::OpProd) = (prodtuple(A.A)...,prodtuple(A.B)...)

# make a tuple from a product, containing only either prefactor types, expectation value types, or operators
for (name,types) in [(:pref,(scal,param,δ)),(:exp,(ExpVal,Corr)),(:op,(BosonDestroy,BosonCreate,FermionDestroy,FermionCreate,σ,σplus,σminus))]
    types = Union{types...}
    @eval $(Symbol(name,:iter))(A::Operator,swallow_∑::Bool) = Iterators.filter(x->x isa $types, proditer(A,swallow_∑))
    name = Symbol(name,:tuple)
    @eval $name(A::$types) = (A,)
    @eval $name(A::Operator) = ()
    @eval $name(A::OpSum) = throw(ArgumentError("cannot get $($name) for OpSum!"))
    @eval $name(A::OpSumAnalytic) = $name(A.A)
    @eval $name(A::OpProd) = ($name(A.A)...,$name(A.B)...)
end
prodtuples(A::Operator) = (preftuple(A), exptuple(A), optuple(A))

sumiter(A::OpSum) = (s*t for (t,s) in A.terms)
sumiter(A::Operator) = (A,)

sumtuple(A::Operator) = (sumiter(A)...,)

isscalar(A::Scalar) = true
isscalar(A::BaseOperator) = false
isscalar(A::OpProd) = all(map(isscalar,proditer(A)))
isscalar(A::OpSum) = all(map(isscalar,sumiter(A)))
isscalar(A::OpSumAnalytic) = isscalar(A.A)
