# MultiPolynomials

My own multivariable polynomials package for Julia.

A reference implementation (not very efficient) of [Buchberger's algorithm](https://en.wikipedia.org/wiki/Buchberger%27s_algorithm) for computing [Gröbner bases](https://en.wikipedia.org/wiki/Gr%C3%B6bner_basis) is included. But there is also an interface to the
highly efficient [FGb library](http://www-polsys.lip6.fr/~jcf/FGb/index.html), see the example below.

##Installation
In a Julia notebook type
```julia
Pkg.clone("https://github.com/HaraldHofstaetter/MultiPolynomials.jl")
Pkg.build("MultiPolynomials")
```
##Examples
To get easy access to the examples, copy them into the home directory:
```julia
cp(joinpath(homedir(), ".julia/v0.4/MultiPolynomials/examples/"), joinpath(homedir(), "MultiPolynomials_examples"), remove_destination=true)
```
Then 'MultiPolynomials_examples' will be listed in the JuliaBox home screen. The examples contain among others
+ [FGb.ipynb](https://github.com/HaraldHofstaetter/MultiPolynomials.jl/blob/master/examples/FGb.ipynb)

