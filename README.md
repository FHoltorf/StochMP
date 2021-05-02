# StochMP

This repository contains the first draft of a Julia package that implements the computation of *hard* bounds on the moments (and related statistics) associated with stochastic jump and diffusion processes via convex optimization. The algorithmic approach is described in our preprint ["Tighter Bounds on transient moments of stochastic chemical systems"](https://arxiv.org/abs/2104.01309) (Holtorf and Barton, 2021) where we put an emphasis on the application to stochastic reaction networks. Currently, a more versatile, extended version of this package which is better integrated within the Julia ecosystem is under development (see [MarkovBounds.jl](https://github.com/FHoltorf/MarkovBounds.jl)).

The core of the code resides in the module StochMP which implements routines to assemble moment problems and the associated semidefinite programming relaxations specifically for bounding transient (or stationary) moments of jump and diffusion processes. In the case of jump processes which model stochastic chemical systems, these routines for example require only high level input data such as

* Stoichiometry matrix
* Reaction propensities 
* Reaction invariants as defined by a list of dependent species
* Polynomial inequalities that are feasible on the reachable set of the reaction network

characterizing the reaction network under investigation. Note, however, that *only* polynomial problem data is supported, i.e., arrival rates/reaction propensities, drift coefficient and diffusion matrix must be polynomial functions of the system state. 

### Examples

This repository contains several examples for the use of StochMP for the analysis of stochastic chemical systems. All examples are discussed in detail in our preprint ["Tighter Bounds on transient moments of stochastic chemical systems"](https://arxiv.org/abs/2104.01309) (Holtorf and Barton, 2021). Please refer to the preprint for more background information.

### How to use this code?

If you would like to use this code or run any of the example code, make sure the following dependencies of StochMP are installed

* LinearAlgebra.jl

* JuMP.jl
* SumOfSquares.jl
* DynamicPolynomials.jl
* MultivariatePolynomials.jl

Moreover, be sure to add the path of StochMP to Julia's LOAD_PATH:

```julia
	push!(LOAD_PATH, $PATH_TO_STOCHMP)
```



### References

Holtorf, Flemming, and Paul I. Barton. "Tighter bounds on transient moments of stochastic chemical systems." *arXiv preprint arXiv:2104.01309* (2021).







