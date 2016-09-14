cd(dirname(@__FILE__))

if (!ispath("bin"))
    run(`mkdir bin`)
end


if (!ispath("call_FGb"))
    download("http://www-polsys.lip6.fr/~jcf/FGb/C/@downloads/call_FGb6.maclinux.x64.tar.gz",
    "./call_FGb6.maclinux.x64.tar.gz")
    run(`tar xzvf call_FGb6.maclinux.x64.tar.gz`)
end

cd(joinpath(dirname(@__FILE__), "src"))
run(`make call_fgb`)
run(`mv call_fgb ../bin`)

run(`make call_giac`)
run(`mv call_giac ../bin`)
