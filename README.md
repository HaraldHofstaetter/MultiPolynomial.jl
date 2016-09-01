# MultiPolynomials

My own multivariable polynomials package for Julia, including an implementation of [Buchberger's algorithm](https://en.wikipedia.org/wiki/Buchberger%27s_algorithm) for computing [Gröbner bases](https://en.wikipedia.org/wiki/Gr%C3%B6bner_basis).

##Installation
In a Julia notebook type
```julia
Pkg.clone("https://github.com/HaraldHofstaetter/MultiPolynomials.jl")
```
##Examples
To get easy access to the examples, copy them into the home directory:
```julia
cp(joinpath(homedir(), ".julia/v0.4/MultiPolynomials/examples/"), joinpath(homedir(), "MultiPolynomials_examples"), remove_destination=true)
```
