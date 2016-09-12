## Reference implementations of algorithms from
## D.Cox, J.Little, D.O'Shea: Ideals, Varieties, and Algorithms (3rd ed)  


lexless{K<:Field,N}(t1::Term{K,N}, t2::Term{K,N}) = lexless(t1.index, t2.index)

function grlexless{N}(i1::Index{N}, i2::Index{N}) 
    d1 = deg(i1)
    d2 = deg(i2)
    (d1<d2) || ((d1==d2) && lexless(i1, i2))
end
grlexless{K<:Field,N}(t1::Term{K,N}, t2::Term{K,N}) = grlexless(t1.index, t2.index)

function grevlexless{N}(i1::Index{N}, i2::Index{N}) 
    d1 = deg(i1)
    d2 = deg(i2)
    (d1<d2) || ((d1==d2) && lexless(reverse(i2), reverse(i1)))
end
grevlexless{K<:Field,N}(t1::Term{K,N}, t2::Term{K,N}) = grevlex(t1.index, t2.index)

_order = lexless


function set_order_lex()
    global _order = lexless
end

function set_order_grlex()
    global _order = grlexless
end

function set_order_grevlex()
    global _order = grevlexless
end

_grevlexelim_l = 0
function grevlexelimless{N}(i1::Index{N}, i2::Index{N}) 
    d1 = sum(i1[1:_grevlexelim_l])
    d2 = sum(i2[1:_grevlexelim_l])
    (d1<d2) || ((d1==d2) && grevlexless(i1, i2))
end
grevlexelimless{K<:Field,N}(t1::Term{K,N}, t2::Term{K,N}) = grevlexelimless(t1.index, t2.index)

function set_order_grevlexelimless(l)
    global _order = grevlexelimless
    global _grevlexelim_l = l
end

isless{K<:Field,N}(t1::Term{K,N}, t2::Term{K,N}) = _order(t1.index, t2.index)
<{K<:Field,N}(t1::Term{K,N}, t2::Term{K,N}) = _order(t1.index, t2.index)
>{K<:Field,N}(t1::Term{K,N}, t2::Term{K,N}) = _order(t2.index, t1.index)

LT{K<:Field,N}(f::Polynomial{K,N}) = maximum(terms(f))
LT{K<:Field,N}(tt::PolyAsArray{K,N}) = maximum(tt)
LC{K<:Field,N}(f::Polynomial{K,N}) = LT(f).coefficient
LC{K<:Field,N}(tt::PolyAsArray{K,N}) = LT(tt).coefficient
LM{K<:Field,N}(f::Polynomial{K,N}) = Term(LT(f).index, one(K))
LM{K<:Field,N}(tt::PolyAsArray{K,N}) = Term(LT(tt).index, one(K))
multidegree{K<:Field,N}(f::Polynomial{K,N}) = index(LT(f))
multidegree{K<:Field,N}(tt::PolyAsArray{K,N}) = index(LT(tt))




function divrem{K<:Field,N}(f::Polynomial{K,N}, ff::Polynomial{K,N} ...)
    # Algorithm from Theorem 2.3 in 
    # D.Cox, J.Little, D.O'Shea: Ideals, Varieties, and Algorithms (3rd ed)  
    s = length(ff)
    zero = Polynomial{K,N}()
    r = copy(zero)
    a = [copy(zero) for i=1:s]
    p = copy(f)
    while p!=zero 
        i = 1
        divisionoccurred = false
        while i<=s && !divisionoccurred
            ltf = LT(ff[i])
            ltp = LT(p)
            if divides(ltf, ltp)
                a[i] = a[i] + ltp/ltf
                p = p - (ltp/ltf)*ff[i]
                divisionoccurred = true
            else    
                i = i + 1
            end
        end    
        if !divisionoccurred
            ltp = LT(p)
            r = r + ltp
            p = p - ltp
        end    
    end    
    a, r
end

function rem{K<:Field,N}(f::Polynomial{K,N}, ff::Polynomial{K,N} ...)
    #The same as a divrem without computing a
    s = length(ff)
    zero = Polynomial{K,N}()
    r = copy(zero)
    p = copy(f)
    while p!=zero 
        i = 1
        divisionoccurred = false
        while i<=s && !divisionoccurred
            ltf = LT(ff[i])
            ltp = LT(p)
            if divides(ltf, ltp)
                p = p - (ltp/ltf)*ff[i]
                divisionoccurred = true
            else    
                i = i + 1
            end
        end    
        if !divisionoccurred
            ltp = LT(p)
            r = r + ltp
            p = p - ltp
        end    
    end    
    r
end

LCM{K<:Field,N}(t1::Term{K,N}, t2::Term{K,N}) = Term{K,N}(Index{N}(([max(t1.index[k], t2.index[k]) for k=1:N]...)), one(K))

function S_polynomial{K<:Field,N}(f::Polynomial{K,N}, g::Polynomial{K,N})
    ltf = LT(f)
    ltg = LT(g)
    lcm = LCM(ltf, ltg)
    (lcm/ltf)*f - (lcm/ltg)*g
end

function groebner{K<:Field, N}(ff::Polynomial{K,N} ...)
    # Algorithm from Theorem 2.11 in 
    # D.Cox, J.Little, D.O'Shea: Ideals, Varieties, and Algorithms (3rd ed)
    zero = Polynomial{K,N}()
    s = length(ff)
    B = Set(vcat([[(i,j)  for i=1:j-1] for j=1:s]...))
    G = [ff[k]/LC(ff[k]) for k=1:s]
    sort(G, lt=(x,y)->LT(x)<LT(y)) # Sort G[i] such that their leading terms are increasing
    t = s 
    lt = [LT(G[k]) for k=1:s]
    while !isempty(B)
        first = true
        (i,j) = (0,0)
        lcm = 0
        for (i1,j1) in B # Select (i,j)∈B of such that LCM(lt[i], lt[j]) is as small as possible.
            if first
                (i,j) = (i1,j1)
                lcm = LCM(lt[i1], lt[j1])
                first = false
            else 
                lcm1 = LCM(lt[i1], lt[j1])
                if lcm1<lcm
                    (i,j) = (i1,j1)
                    lcm = lcm1
                end
            end
        end        
        #println((i,j), " selected")
        if lcm.index != (lt[i]*lt[j]j).index
            criterion = false        
            for k = 1:t
                if k==i || k==j
                    continue
                end
                if !((min(i,k),max(i,k)) in B || (min(j,k),max(j,k)) in B) && divides(lt[k],lcm)
                    criterion = true
                    break
                end
            end
            if !criterion
                S = rem(S_polynomial(G[i],G[j]), G...)
                if S!=zero
                    t = t+1
                    push!(G, S/LC(S))
                    push!(lt, LT(G[t]))
                    union!(B, [(k,t) for k=1:t-1])
                    #println("\t",S)
                end
            end
        end
        delete!(B, (i,j))
    end
    
    # generate minimal Groebner basis
    to_be_deleted = Int[]
    for i=1:length(G)
        if i in to_be_deleted
            continue
        end
        for j=1:length(G)
            if j==i || j in to_be_deleted || i in to_be_deleted 
                continue
            end        
            if divides(LT(G[j]), LT(G[i]))
                push!(to_be_deleted, i)
            end
         end    
    end
    deleteat!(G, to_be_deleted)
    for i=1:length(G)
        G[i] = G[i]/LC(G[i])
    end
    
    # generate reduced Groebner basis
    t = length(G)
    for i=1:t
        G[i] = rem(G[i], vcat(G[1:i-1],G[i+1:t])...)
    end
        
    G
end

