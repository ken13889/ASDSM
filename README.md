# ASDSM

Paper on Arxiv:  A domain splitting strategy for solving PDEs

In this work we develop a novel domain splitting strategy for the solution of partial differential equations. 
Focusing on a uniform discretization of the $d$-dimensional advection-diffusion equation, 
our proposal is a two-level algorithm that merges the solutions obtained from the discretization of the 
equation over highly anisotropic submeshes to compute an initial approximation of the fine solution. 
The algorithm then iteratively refines the initial guess by leveraging the structure of the residual.
erforming costly calculations on anisotropic submeshes enable us to reduce the dimensionality of the problem by one, 
and the merging process, which involves the computation of solutions over disjoint domains, allows for parallel implementation.
