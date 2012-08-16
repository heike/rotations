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
Given a sample of n misorientation angles, use the `genR` function to simulate a matrix from the corresponding matrix distribution.  The default *central orientation*, or mean, is the 3-by-3 identity matrix but any matrix in SO(3) can be the central orientation and specified by using the `S` option in `genR`.  The class of data resulting from this method is SO3.

## Central orientation estimation
Four estimators of the central orientation are available from these two functions:

* `median(Rs)` finds the matrix in SO(3) that minimizes the sum of first order distances from each sample point.  The `type` option allows the user to choose the Riemannian or Euclidean geometry to find the corresponding median-type estimator called the *intrinsic* or *projected* median respectively.  The algorithm used to compute the intrinsic median is described in Hartley et al. (2011).   Because the algorithms are iterative the `epsilon` and `maxIter` values control when convergence is reached or when to stop. 

* `mean(Rs)` finds the matrix in SO(3) that minimizes the sum of second order distances from each sample point.  The `type`, `epsilon` and `maxIter` options work as described in `median` above.  For a discussion on the projected mean estimator see Moakher (2002) and the intrinsic mean is discussed in Manton (2004).

## Visualization
Based on `ggplot2` by Wickham (2009) we have developed a method for visualizing random rotations in SO(3) in three-dimensions.  The function `eyeBall` plots the data on the sphere along with all four estimators of the central orientation if the option `show.estimates` is set to `TRUE`.  By default the data are centered around  the three-dimensional identity matrix, but any matrix in SO(3) can be used as the center with the `center` option.  Finally, since only one column of the three is displayed at a time the `column` option allows the user to choose which column to display.  Additional options can be passed to the call to `qplot` too.