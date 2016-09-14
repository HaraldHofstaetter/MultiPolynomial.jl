# MultiPolynomials

My own multivariate polynomials package for Julia.

A reference implementation (not very efficient) of [Buchberger's algorithm](https://en.wikipedia.org/wiki/Buchberger%27s_algorithm) for computing [Gröbner bases](https://en.wikipedia.org/wiki/Gr%C3%B6bner_basis) is included. But there are also  interfaces to the
highly efficient [FGb](http://www-polsys.lip6.fr/~jcf/FGb/index.html)
and [Giac](https://www-fourier.ujf-grenoble.fr/~parisse/giac.html) libraries, see the example below.

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

