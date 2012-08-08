# rotations
=======================================================

This package implements tools for working with rotational data: it allows simulation from the most commonly used distributions on the sphere, it estimates different mean and median type estimators for the main direction of a number of rotations, it provides (bootstrap) confidence regions for the estimates, and it allows to visualize rotational data.

## Installation

To install in your R session use:
```
library(devtools)
install_github("rotations","heike")
library(rotations)
```
## Data generation

Simulating data in SO(3) is a two step process.  One first simulates a set of angels from an symmetric distribution about zero.  Then one can simulate a matrix from the uniform-axis random spin class of distributions based on those angles.

### Misorientation angle simulation
There are three angular distributions to simulate from.  They are all symmetric about zero and bounded between negative pi and pi.

* `rvmises(n,kappa)` simulates data from the von Mises circular distributions
* `rcayley(n,kappa)` simulates data from the Cayley distributions
* `rfisher(n,kappa)` simulates data from the von Mises-Fisher matrix distribution

### Matrix simulation
Given a sample of n misorientation angles, use the `genR` function to simulate a matrix from the corresponding matrix distribution.  The default *central orientation*, or mean, is the 3-by-3 identity matrix but any matrix in SO(3) can be the central orientation and specified by using the `S` option in `genR`.

## Central orientation estimation