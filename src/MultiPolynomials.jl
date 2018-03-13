module MultiPolynomials

import 
Base: (*), /, ^, +, -, <, >, rem, divrem, diff, isless, lexless, string, show, write, convert, zero, one

export Field, Index, Term, Polynomial, PolyAsArray, @variables
export terms, polynomial, coeff, deg, index, divides, evaluate
export extract_univariate_polynomial, convert_base_field
export grlexless, grevlexless
export set_order_lex, set_order_grlex, set_order_grevlex
export grevlexelimless, set_order_grevlexelimless
export LT, LC, LM, multidegree
export LCM, S_polynomial, groebner 
export groebner_fgb, call_fgb
export fgb_gbasis, fgb_gbasis_elim
export call_giac, giac_gbasis, giac_greduce

const Field = Number

Index{N} = NTuple{N,UInt16}

immutable Term{K<:Field,N}
    index::Index{N}
    coefficient::K
end

Polynomial{K<:Field,N} = Dict{Index{N},K}

PolyAsArray{K<:Field,N} = Array{Term{K,N},1}

function polynomial{K<:Field, N}(terms::PolyAsArray{K, N})
    Polynomial{K, N}([(t.index, t.coefficient) for t in terms])
end

function terms{K<:Field,N}(f::Polynomial{K, N})
    PolyAsArray{K,N}(Term{K,N}[Term{K,N}(key, value) for (key,value) in f])
end

convert{K<:Field,N}(::Type{Polynomial{K,N}}, t::Term{K,N}) = 
    Polynomial{K,N}([(t.index, t.coefficient)])

zero{K,N}(::Type{Polynomial{K,N}}) = Polynomial{K,N}([])
one{K,N}(::Type{Polynomial{K,N}}) = Polynomial{K,N}([((zeros(K,N)...), 1)])

deg(i::Index) = Int(sum(i))
deg(t::Term) = deg(t.index)

function deg(f::Polynomial)
    d = 0
    for key in keys(f)
        d = max(d, deg(key))
    end
    d
end

function deg{K<:Field, N}(f::PolyAsArray{K, N})
    d = 0
    for t in f
        d = max(d, deg(t))
    end
    d
end

index(t::Term) = map(x->convert(Int,x),t.index)

function coeff{K<:Field,N}(t::Term{K,N}, f::Polynomial{K,N})
    get(f, t.index, zero(K))/t.coefficient
end

function divides{K<:Field,N}(t1::Term{K,N}, t2::Term{K,N})
    for i=1:N
        if t1.index[i]>t2.index[i]
            return false
        end
    end
    return true
end    

function divides{K<:Field,N}(t::Term{K,N}, f::Polynomial{K,N})
    for (i,c) in f
        if !divides(t, Term{K,N}(i,c))
            return false
        end    
    end
    return true
end

function divides{K<:Field,N}(t::Term{K,N}, tt::PolyAsArray{K,N})
    for s in tt
        if !divides(t, s)
            return false
        end    
    end
    return true
end

+{K<:Field, N}(f::Polynomial{K, N}) = copy(f)
+{K<:Field, N}(t::Term{K, N}) = t
-{K<:Field, N}(t::Term{K, N}) = Term(t.index, -t.coefficient)

function -{K<:Field, N}(f::Polynomial{K, N})
    Polynomial{K,N}([(i,-c) for (i,c) in f])
end        
        
function +{K<:Field, N}(f::Polynomial{K, N}, g::Polynomial{K, N})
    r = copy(f)
    for (key, val) in g
        if haskey(r, key)
            r[key] += val
        else
            r[key] = val
        end    
        if r[key]==zero(K)
            delete!(r, key) 
        end     
    end
    r
end

function -{K<:Field, N}(f::Polynomial{K, N}, g::Polynomial{K, N})
    r = copy(f)
    for (key, val) in g
        if haskey(r, key)
            r[key] -= val
        else
            r[key] = -val
        end    
        if r[key]==zero(K)
            delete!(r, key) 
        end     
    end
    r
end

function +{K<:Field, N}(f::Polynomial{K, N}, t::Term{K, N})
    r = copy(f)
    key, val = (t.index, t.coefficient)
    if haskey(r, key)
        r[key] += val
    else
        r[key] = val
    end    
    if r[key]==zero(K)
        delete!(r, key) 
    end     
    r
end

function -{K<:Field, N}(f::Polynomial{K, N}, t::Term{K, N})
    r = copy(f)
    key, val = (t.index, t.coefficient)
    if haskey(r, key)
        r[key] -= val
    else
        r[key] = -val
    end    
    if r[key]==zero(K)
        delete!(r, key) 
    end     
    r
end

function -{K<:Field, N}(t::Term{K, N}, f::Polynomial{K, N},)
    r = -f
    key, val = (t.index, t.coefficient)
    if haskey(r, key)
        r[key] += val
    else
        r[key] = val
    end    
    if r[key]==zero(K)
        delete!(r, key) 
    end     
    r
end

+{K<:Field, N}(t::Term{K, N}, f::Polynomial{K, N}) = +(f,t)

+{K<:Field, N}(t::Term{K, N}, s::Term{K, N}) = +(Polynomial{K,N}(t.index=>t.coefficient), s)

+{K<:Field, N}(f::Polynomial{K, N}, a::Field) = +(f, Term{K,N}((zeros(K,N)...), convert(K, a)))

+{K<:Field, N}(a::Field, f::Polynomial{K, N}) = +(f, a)

+{K<:Field, N}(t::Term{K, N}, a::Field) = +(Polynomial{K,N}(t.index=>t.coefficient), convert(K, a))

+{K<:Field, N}(a::Field, t::Term{K, N}) = +(t, a)



-{K<:Field, N}(t::Term{K, N}, s::Term{K, N}) = -(Polynomial{K,N}(t.index=>t.coefficient), s)

-{K<:Field, N}(f::Polynomial{K, N}, a::Field) = -(f, Term{K,N}((zeros(K,N)...), convert(K, a)))

-{K<:Field, N}(a::Field, f::Polynomial{K, N}) = -(Term{K,N}((zeros(K,N)...), convert(K, a)), f)

-{K<:Field, N}(t::Term{K, N}, a::Field) = -(Polynomial{K,N}(t.index=>t.coefficient), convert(K,a))

-{K<:Field, N}(a::Field, t::Term{K, N}) = -(convert(K,a), Polynomial{K,N}(t.index=>t.coefficient))


function _mul{N}(i1::Index{N}, i2::Index{N})
    Index{N}(([i1[k]+i2[k] for k=1:N]...))
end

function *{K<:Field, N}(t::Term{K, N}, a::Field)
    Term{K, N}(t.index, t.coefficient*convert(K, a))
end

*{K<:Field, N}(a::Field, t::Term{K, N}) = *(t, a)

function *{K<:Field, N}(t::Term{K, N}, s::Term{K, N})
    Term{K, N}( _mul(t.index, s.index), t.coefficient*s.coefficient)
end

function *{K<:Field, N}(f::Polynomial{K, N}, a::Field)
    a1 = convert(K, a)
    if a1==zero(K)
        return Polynomial{K, N}([])
    end
    r = copy(f)
    for key in keys(r)
        r[key]*=a1
    end
    r 
end

*{K<:Field, N}(a::Field, f::Polynomial{K, N}) = *(f, a)

function *{K<:Field, N}(f::Polynomial{K, N}, t::Term{K, N})
    Polynomial{K, N}([(_mul(i, t.index), c*t.coefficient) for (i, c) in f])
end

*{K<:Field, N}(t::Term{K, N}, f::Polynomial{K, N}) = *(f, t)

function *{K<:Field, N}(f::Polynomial{K, N}, g::Polynomial{K, N})
    tt = terms(g)
    if length(tt) == 0
        return copy(g)
    end
    r = f*tt[1]
    for k=2:length(tt)
        r += f*tt[k]
    end
    r
end

function _div{N}(i1::Index{N}, i2::Index{N})
    Index{N}(([ (i1[k]>=i2[k] ? i1[k]-i2[k] : -1) for k=1:N]...))
    #deliberately throws error if not divisible
end

function /{K<:Field, N}(t::Term{K, N}, s::Term{K, N})
    Term{K, N}( _div(t.index, s.index), t.coefficient/s.coefficient)
end

function /{K<:Field, N}(f::Polynomial{K, N}, t::Term{K, N})
    Polynomial{K, N}([(_div(i, t.index), c/t.coefficient) for (i, c) in f])
end

function /{K<:Field, N}(t::Term{K, N}, a::Field)
    Term{K, N}(t.index, t.coefficient/convert(K, a))
end

function /{K<:Field, N}(f::Polynomial{K, N}, a::Field)
    a1 = convert(K, a)
    r = copy(f)
    for key in keys(r)
        r[key]/=a1
    end
    r 
end



function _pow{N}(i::Index{N}, e::Integer)
    Index{N}(([e*i[k] for k=1:N]...))
end

function ^{K<:Field, N}(t::Term{K, N}, e::Integer)
    Term{K, N}(_pow(t.index, e), t.coefficient^e)
end

function ^{K<:Field, N}(f::Polynomial{K, N}, e::Integer)
    if e<0
        error("Nonnegative exponent expected")
    elseif e==0
        return Term{K,N}((zeros(K,N)...), 1)
    elseif e==1
        return f
    else
        s = bits(e)
        s = reverse(s[first(search(s,"1")):end])
        r = 0
        q = 0
        isfirst = true
        for k=1:length(s)
            if k==1
                q = f
            else
                q = q*q
            end    
            if s[k:k]=="1"
                if isfirst
                    r = q
                    isfirst = false
                else
                    r = r*q
                end
            end    
        end        
        return r
    end    
end

_varnames = []


macro variables(K, x...)
    N = length(x)
    T = parse(string("Term{", K, ",", N, "}")) 
    q = Expr(:block)    
    if length(x) == 1 && isa(x[1],Expr)
        @assert x[1].head === :tuple "@variables expected a list of symbols"
        x = x[1].args
    end
    i = 1
    for s in x
        index = zeros(UInt16, N)
        index[i] = 1
        @assert isa(s,Symbol) "@variables expected a list of symbols"
        push!(q.args, Expr(:(=), s, Expr(:call, T, Expr(:tuple, index...), 1) ) )
        i = i + 1
    end    
    push!(q.args, Expr(:tuple, x...))
    global _varnames = map(string, x)
    eval(Main, q)
end


function _str(i::Index; latex::Bool=false)
    s = ""
    f = false
    for k = 1:length(i)
        if i[k]>0
            n = ""
            if k<=length(_varnames)
                n = _varnames[k]
            else
                if latex
                    n = string("X_{", k, "}")
                else
                    n = string("X", k)
                end                    
            end                
            if latex
                s = string(s, n, i[k]>1?string("^{", i[k], "}"):"")
            else
                s = string(s, f?"*":"", n, i[k]>1?string("^", i[k]):"")                
            end    
            f = true
        end
    end
    s
end

_str(x; latex::Bool=false) = string(x)

function _str(t::Term; latex::Bool=false)
    return string(sum(t.index)!=0 && abs(t.coefficient)==1?(t.coefficient<0?"-":""):string(_str(t.coefficient, latex=latex),(sum(t.index)==0||latex)?"":"*"), _str(t.index, latex=latex)) 
end

#function _str{K<:Field,N}(tt::PolyAsArray{K,N}; latex::Bool=false)    
function _str{K<:Field,N}(tt::Array{Term{K,N},1}; latex::Bool=false)
#function _str(tt::PolyAsArray; latex::Bool=false)    
    if length(tt)==0
        return "0"
    end
    s = ""
    f = false
    for t in tt
        if f && t.coefficient>0
            s = string(s, "+")
        end
        f = true
        s = string(s, _str(t, latex=latex))        
    end
    s
end

function _str{K<:Field,N}(f::Polynomial{K,N}; latex::Bool=false)
#function _str(f::Polynomial; latex::Bool=false)
#    _str(terms(f), latex=latex)
    tt = terms(f)
    if length(tt)==0
        return "0"
    end
    s = ""
    f = false
    for t in tt
        if f && t.coefficient>0
            s = string(s, "+")
        end
        f = true
        s = string(s, _str(t, latex=latex))        
    end
    s    
end

function _str{T}(q::Rational{T}; latex::Bool=false)
    if latex
        if denominator(q)==1
            return string(numerator(q))
        end
        s = ""
        if numerator(q)<0
            s = "-"
        end
        return string(s, "\\frac{", abs(numerator(q)), "}{", denominator(q),"}")
    else 
        return string(q)
    end    
end

string(t::Term) = _str(t)
show(io::IO, t::Term) = print(io, _str(t))
#show(io::IO, ::MIME"application/x-latex", t::Term) = write(io, "\$", _str(t, latex=true), "\$")
#show(io::IO, ::MIME"text/latex",  t::Term) = write(io, "\$", _str(t, latex=true), "\$")

string{K<:Field,N}(tt::PolyAsArray{K,N}) = _str(tt)
show{K<:Field,N}(io::IO, tt::PolyAsArray{K,N}) = print(io, _str(tt))
#show{K<:Field,N}(io::IO, ::MIME"application/x-latex", tt::PolyAsArray{K,N}) = write(io, "\$", _str(tt, latex=true), "\$")
#show{K<:Field,N}(io::IO, ::MIME"text/latex",  tt::PolyAsArray{K,N}) = write(io, "\$", _str(tt, latex=true), "\$")

#string(f::Polynomial) = _str(f)
#show(io::IO, f::Polynomial) = print(io, _str(f))
string{K<:Field,N}(f::Polynomial{K,N}) = _str(f)
show{K<:Field,N}(io::IO, f::Polynomial{K,N}) = print(io, _str(f))
#show(io::IO, ::MIME"application/x-latex", f::Polynomial) = write(io, "\$", _str(f, latex=true), "\$")
#show(io::IO, ::MIME"text/latex",  f::Polynomial) = write(io, "\$", _str(f, latex=true), "\$")

#show(io::IO, ::MIME"application/x-latex", q::Rational) = write(io, "\$", _str(q, latex=true), "\$")
#show(io::IO, ::MIME"text/latex",  q::Rational) = write(io, "\$", _str(q, latex=true), "\$")

evaluate{K<:Field,N}(t::Term{K,N}, x) = t.coefficient * prod([x[k]^t.index[k] for k=1:N])
evaluate{K<:Field,N}(f::Polynomial{K,N}, x) = sum([evaluate(t, x) for t in terms(f)])

function diff{K<:Field,N}(t::Term{K,N}, k::Int)
    if t.index[k] > 0
        i = [t.index...]
        i[k] -= 1
        return Term{K,N}((i...), t.index[k]*t.coefficient)
    else
        return zero(K)
    end
end

function diff{K<:Field,N}(t::Term{K,N}, v::Term{K,N}) 
    @assert sum(v.index)==1 
    diff(t, findfirst(v.index,1))
end

function diff{K<:Field,N}(f::Polynomial{K,N}, k::Int)
    r = Polynomial{K,N}()
    for (i,c) in f
        r = r + diff(Term{K,N}(i,c), k)
    end
    r	
end

function diff{K<:Field,N}(f::Polynomial{K,N}, v::Term{K,N})
    @assert sum(v.index)==1 
    diff(f, findfirst(v.index,1))
end


function convert_base_field{K<:Field, N}(K_new::DataType, p::Polynomial{K,N})
    Polynomial{K_new,N}([(i, convert(K_new, c))  for (i,c) in p])
end


function extract_univariate_polynomial{K<:Field, N}(p::Polynomial{K,N})
    t = terms(p)
    d = deg(t)
    up = zeros(K, d+1)
    v = 0
    for i = 1:length(t)
        if v==0 && deg(t[i])>0
            v = findfirst(t[i].index)            
        end
        if v>0
            for j=1:N
                if j!=v && t[i].index[j]!=0
                    error("not a univariate polynomial")
                end
            end
            up[t[i].index[v]+1] = t[i].coefficient
        else
            up[1] = t[i].coefficient
        end
    end
    up
end

include("groebner_et_al.jl")
include("interfaces.jl")


end #MultiPolynomials
