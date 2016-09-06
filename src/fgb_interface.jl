function _write_fgb_input{K<:Field, N}(si, ff::Polynomial{K,N} ...)
    println(si, N, " ", length(ff))
    for f in ff
        println(si, length(f))
        h = lcm([den(c) for (i,c) in f])
        for (i,c) in f
            for k=1:N
                print(si, i[k], " ")
            end
            println(si, num(c*h))
        end
    end    
end

function _read_fgb_output(so, K::DataType)
    N = parse(Int, readuntil(so,' '))
    n_poly = parse(Int, readuntil(so,'\n')[1:end-1])
    F = Polynomial{K,N}[]
    for p=1:n_poly
        n_terms = parse(Int, readuntil(so,'\n')[1:end-1])
        f = Polynomial{K,N}([((UInt16[parse(UInt16, readuntil(so,' ')) 
                                      for i=1:N]...),
                    convert(K,parse(BigInt, readuntil(so,'\n')[1:end-1])))
                    for t=1:n_terms])
        push!(F,f)              
    end
    F
end

function groebner_fgb{K<:Field, N}(ff::Polynomial{K,N} ...)
    call_fgb = joinpath(dirname(@__FILE__),
               "..", "deps", "bin", "call_fgb")
    (so,si,pr) = readandwrite(`$call_fgb`)

    _write_fgb_input(si, ff...)

    G = _read_fgb_output(so, K)

    close(so)
    close(si)
    close(pr) 

    G
end

