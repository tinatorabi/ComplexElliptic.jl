# ComplexElliptic

A package for evaluating [Jacobi elliptic functions](https://dlmf.nist.gov/22.2) for complex arguments.


The algorithm is the descending Landen transformation, which is explained in detail in L. Howell's MIT PhD thesis. Additional transformations where used from Gradshteyn & Ryzhik "Table of Integrals, Series, and Products", Fifth Edition, and M. Abramowitz and I.A. Stegun, "Handbook of Mathematical Functions," Dover, 1965.

This algorithm was originally implemented in MATLAB by Toby Driscoll. 

```jlcon
julia> import ComplexElliptic.ellipjc

julia> ellipjc(3+4im, 2)
(3.8611134365506516 - 27.032884005488256im, -19.4248392253911 + 2.7781705335555387im, 0.689812164575172 - 0.20129607887149814im)
```

The package can be installed as follows:

```jlcon
julia> ] add ComplexElliptic

```
