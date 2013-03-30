# rotations
=======================================================

This package implements tools for working with rotational data: it allows simulation from the most commonly used distributions on the sphere, it estimates different mean and median type estimators for the main direction of a number of rotations, it provides (bootstrap) confidence regions for the estimates, and it allows to visualize rotational data.

## Installation

To install in your R session use:
```
library(devtools)
install_github("rotations","heike") #Do this once
library(rotations)
```
## Data generation

Three-dimensional rotations are determined uniquely by three numbers: two to define an axis and one to specify the rotation of all three dimensions about that axis.  This package follows the uniform-axis random-spin (*UARS*) framework of Bingham et al (2009a) in simulating the axis and angle of rotation independently.  As the UARS name suggests, the axis, sometimes called the *misorientation axis* is chosen uniformly on the unit sphere and the angle of rotation, sometimes called the *misorientation angle*, is choose according to some circular distribution symmetric about zero and in the interval $[-\pi,\pi)$.  The distribution chosen is usually referred to as the *angular distribution* as it describes the distributions of the misorientation angles.

To date there are four angular distributions from which to choose and two parameterizations of rotations.  The former are the uniform, von Mises circular, Cayley and the von Mises-Fisher matrix distributions.  One can describe a rotation in a 3-by-3 *rotation matrix* or a four dimensional vector with unit length called a *quaternion*.  Each will be described below. 

### Misorientation angle simulation
There are four angular distributions to simulate from.  They are all symmetric about zero and in the interval $[-\pi,\pi)$.

* `rhaar(n)`  simulates data from the uniform distribution on the sphere, also called the Haar measure
* `rvmises(n,kappa)` simulates data from the von Mises circular distributions
* `rcayley(n,kappa)` simulates data from the Cayley distributions
* `rfisher(n,kappa)` simulates data from the von Mises-Fisher matrix distribution

### Rotation simulation
Given a sample of n angles, use the `genR` function to simulate a rotation from the corresponding rotational distribution.  The default *central orientation*, or mean, is the identity but any rotation can be the central orientation and is specified by using the `S` option.  The three parameterizations are considered next.

#### Matrix parameterization
When the `space` argument in `genR` is set to `SO3` then a 3-by-3 matrix is formed as described by Rodrigues' formula.  Given axis $u$ and angle $r$ the *rotation matrix* $R(u,r)$, or simply $R$, is formed as follows:

$$R=uu^\top+(I_{3\times 3}-uu^\top)\cos(r)+\Phi(u)\sin(r)$$

where 
$$\Phi(u)=\begin{pmatrix} 0 & -u_3 & u_2 \\ u_3 & 0 & -u_1\\ -u_2 & u_1 & 0\\ \end{pmatrix}.$$ 

If the `S` option is left to the default then $R$ is returned, otherwise $SR$ is returned

#### Quaternion parameterization
An alternative to the rotation matrix is the quaternion, which is chosen by setting the `space` option to `Q4`.  This formulation is typically popular with mathematicians and engineers.  Given the axis $u$ and angle $r$ the quaternion is formed as follows

$$q=(\cos(r/2),\sin(r/2)u)^\top.$$

It's easy to see that a vector formed in such a fashion has unit length as $q^\top q=\cos(r/2)^2+\sin(r/2)^2(u_1^2+u_2^2+u_3^2)=\cos(r/2)^2+\sin(r/2)^2=1$.

## Central orientation estimation
For any of the three parameterizations above four estimators of the central orientation are available from the following two functions.

* `median(Rs)` finds the rotation that minimizes the sum of first order distances from each sample point.  The `type` option allows the user to choose the Riemannian or Euclidean geometry to find the corresponding median-type estimator called the *intrinsic* or *projected* median respectively.  The algorithm used to compute the intrinsic median is described in Hartley et al. (2011).   Because the algorithms are iterative the `epsilon` and `maxIter` values control when convergence is reached or when to stop. 

* `mean(Rs)` finds the rotation that minimizes the sum of second order distances from each sample point.  The `type`, `epsilon` and `maxIter` options work as described in `median` above.  For a discussion on the projected mean estimator see Moakher (2002) and the intrinsic mean is discussed in Manton (2004).

## Confidence Regions
For the projected mean estimator three methods to estimate a confidence region are available using the `region` function.  The three versions available are differentiated by the `type` option named for the first author in the paper the method first appeared: `Rancourt`, `Zhang` and `Fisher`.  The `Zhang` method is a bootstrap method and additional arguments are required to specify the number of bootstrap samples used (`m`), which defaults at 300.  The `Zhang` method can also be used to compute a confidence region for the intrinsic median by specifying `estimator=median` instead of the default `estimator=mean`.

## Visualization
Based on `ggplot2` by Wickham (2009) we have developed a method for visualizing random rotations in SO(3) in three-dimensions.  The function `plot` plots the data on the sphere along with all four estimators of the central orientation if the option `show_estimates` is set to `TRUE`.  By default the data are centered around  the three-dimensional identity matrix, but any matrix in SO(3) can be used as the center with the `center` option.  Finally, since only one column of the three is displayed at a time the `column` option allows the user to choose which column to display.  Additional options can be passed to the call to `qplot` too.