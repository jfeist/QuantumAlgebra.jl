replace_index(A::scal,iold,inew) = A
replace_index(A::param,iold,inew) = param(A.name,A.state,(n-> n==iold ? inew : n).(A.inds))
replace_index(A::ExpVal,iold,inew) = ExpVal(replace_index(A.A,iold,inew))
replace_index(A::Corr,iold,inew) = Corr(replace_index(A.A,iold,inew))
for op in (δ,a,adag,f,fdag,σminus,σplus)
    @eval replace_index(A::$op,iold,inew) = $op((n-> n==iold ? inew : n).(A.inds))
end
replace_index(A::σ,iold,inew) = σ(A.a,(n-> n==iold ? inew : n).(A.inds))
replace_index(A::OpProd,iold,inew) = replace_index(A.A,iold,inew)*replace_index(A.B,iold,inew)
replace_index(A::OpSum,iold,inew) = replace_index(A.A,iold,inew) + replace_index(A.B,iold,inew)
replace_index(A::OpSumAnalytic,iold,inew) = begin
    (A.ind==iold || A.ind==inew) && throw(ArgumentError("replace_index in OpSumAnalytic cannot have iold ($iold) or inew ($inew) be the same as the sum index ($(A.ind))!"))
    tmp = replace_index(A.A,iold,inew)
    # we have to be careful here - when replacing an index inside the expression, there might be reorderings that
    # have be done with the same index explicitly to make sure that the commutator shows up
    # so take out the term with the same index explicitly, and redo that one term
    OpSumAnalytic(A.ind,tmp) - replace_index(tmp,A.ind,inew) + replace_index(replace_index(A.A,A.ind,iold),iold,inew)
end
replace_index(A::Operator,reps) = begin
    Anew = A
    for (iold,inew) in reps
        Anew = replace_index(Anew,iold,inew)
    end
    Anew
end

function exchange_inds(A::Operator,i1,i2)
    ss = gensym()
    A = replace_index(A,i1,ss)
    A = replace_index(A,i2,i1)
    A = replace_index(A,ss,i2)
    A
end

distribute_indices!(inds,A::scal) = A
distribute_indices!(inds,A::param) = param(A.name,A.state,(popfirst!(inds) for _ in A.inds)...)
distribute_indices!(inds,A::ExpVal) = ExpVal(distribute_indices!(inds,A.A))
distribute_indices!(inds,A::Corr) = Corr(distribute_indices!(inds,A.A))
for op in (a,adag,f,fdag,σminus,σplus)
    @eval distribute_indices!(inds,A::$op) = $op((popfirst!(inds) for _ in A.inds)...)
end
distribute_indices!(inds,A::σ) = σ(A.a,(popfirst!(inds) for _ in A.inds)...)
distribute_indices!(inds,A::OpProd) = distribute_indices!(inds,A.A)*distribute_indices!(inds,A.B)
distribute_indices!(inds,A::OpSum) = distribute_indices!(inds,A.A) + distribute_indices!(inds,A.B)
# on purpose do not define this for OpSumAnalytic or δ

indextuple(A::scal)::OpIndices = ()
indextuple(A::Union{param,δ,a,adag,f,fdag,σ,σminus,σplus})::OpIndices = A.inds
indextuple(A::Union{OpProd,OpSum})::OpIndices = (indextuple(A.A)...,indextuple(A.B)...)
indextuple(A::Union{ExpVal,Corr})::OpIndices = indextuple(A.A)
indextuple(A::OpSumAnalytic)::OpIndices = (A.ind,indextuple(A.A)...)
indexset(A) = Set{OpIndex}(indextuple(A))
sumindextuple(A::Operator)::OpIndices = ()
sumindextuple(A::Union{OpProd,OpSum})::OpIndices = (sumindextuple(A.A)...,sumindextuple(A.B)...)
sumindextuple(A::Union{ExpVal,Corr})::OpIndices = sumindextuple(A.A)
sumindextuple(A::OpSumAnalytic)::OpIndices = (A.ind,sumindextuple(A.A)...)
sumindexset(A) = Set{OpIndex}(sumindextuple(A))

"`extindices(A::Operator)` return externally visible indices of an expression"
extindices(A::Operator) = [ind for ind in indextuple(A) if !in(ind,sumindextuple(A))]

"`symmetric_index_nums(A::Operator)` return sequence of numbers of exchange-symmetric indices"
function symmetric_index_nums(A::Operator)
    inds = extindices(A)
    Nsyms = [1]
    for ii=2:length(inds)
        if A == exchange_inds(A,inds[ii-1],inds[ii])
            Nsyms[end] += 1
        else
            push!(Nsyms,1)
        end
    end
    Nsyms
end

"sum indices have no semantic meaning, so rename them in case they happen to occur in the other expression"
function ensure_compatible_sumind(S::OpSumAnalytic,A::Operator)
    Ainds = indextuple(A)
    if S.ind in Ainds
        oldinds = Set{OpIndex}((Ainds...,indextuple(S)...))
        m = match(r"(.*)_([0-9]+)",string(S.ind))
        indstem, ii = (m === nothing) ? (string(S.ind), 1) : (m.captures[1], 1+parse(Int,m.captures[2]))
        while (newind = Symbol(indstem,:_,ii)) in oldinds
            ii += 1
        end
        OpSumAnalytic(newind,replace_index(S.A,S.ind,newind))
    else
        S
    end
end
