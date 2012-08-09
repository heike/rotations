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
Four estimators of the central orientation and available:

* `rmedian(Rs)` estimates the matrix in SO(3) that minimizes the sum of first order Euclidean distances from each sample point.  It's referred to as the projected median because it finds the median in the traditional sense then projects it into SO(3).
* `HartleyL1(Rs)` finds the matrix in SO(3) minimizing the sum of first order Riemannian distances and is called the geometric median.  The algorithm used was proposed in Hartley et al. (2011)
* `arith.mean(Rs)` is called the projected mean and is defined as the matrix minimizing the sum of squared Euclidean distances.  Similar to `rmedian` it finds the mean matrix then projects it to SO(3).  See Moakher (2002) for details.
* `MantonL2(Rs)` estimates the central orientation by finding the matrix that minimizes the sum of squared Riemannian distances.  It is called the geometric mean and the algorithm used was proposed by Manton (2004).

## Visualization
Based on `ggplot2` by Hadley (2009) we have developed a method for visualizing random rotations in SO(3) in three-dimensions.  The function `eyeBall` plots the data on the sphere along with all four estimators of the central orientation if the option `show.estimates` is set to `TRUE`.  By default the data are centered around  the three-dimensional identity matrix, but any matrix in SO(3) can be used as the center with the `center` option.  Finally, since only one column of the three is displayed at a time the `column` option allows the user to choose which column to display.  Additional options can be passed to the call to `qplot` too.