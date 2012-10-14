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

To date there are four angular distributions from which to choose and three parameterizations of rotations.  The former are the uniform, von Mises circular, Cayley and the von Mises-Fisher matrix distributions.  One can describe a rotation in a 3-by-3 *rotation matrix*, a four dimensional vector with unit length called a *quaternion* or a three-dimensional vector called an *Euler angle*.  Each will be described below. 

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

$$q=(\cos(r/2),\sin(r/2)u).$$

It's easy to see that a vector formed in such a fashion has unit length as $q^\top q=\cos(r/2)^2+\sin(r/2)^2(u_1^2+u_2^2+u_3^2)=\cos(r/2)^2+\sin(r/2)^2=1$.

#### Euler angles
We use the $z-x-z$ forumaltion of Euler angles, i.e. the triplet $r=(\alpha,\beta,\gamma)$ is a $z-x-z$ Euler angle if $\alpha,\gamma\in[0,2\pi]$ and $\beta\in[0,\pi]$ where the rotation $r$ is formed by first rotating the reference coordinate orientation ($S$ from before) counterclockwise about the $z$-axis $\alpha$ radians, then rotate the $x$-axis $\beta$ radians and finally about the $z$-axis again $\gamma$ radians.

## Central orientation estimation
For any of the three parameterizations above four estimators of the central orientation are available from the following two functions.  In truth, the algorithms are written for the matrix parameterization so the quaternion and Euler angle versions reparameterize to matrices, find the mean or median then reperamaterize to the original form.

* `median(Rs)` finds the rotation that minimizes the sum of first order distances from each sample point.  The `type` option allows the user to choose the Riemannian or Euclidean geometry to find the corresponding median-type estimator called the *intrinsic* or *projected* median respectively.  The algorithm used to compute the intrinsic median is described in Hartley et al. (2011).   Because the algorithms are iterative the `epsilon` and `maxIter` values control when convergence is reached or when to stop. 

* `mean(Rs)` finds the rotation that minimizes the sum of second order distances from each sample point.  The `type`, `epsilon` and `maxIter` options work as described in `median` above.  For a discussion on the projected mean estimator see Moakher (2002) and the intrinsic mean is discussed in Manton (2004).

## Confidence Regions
For each of the four estimators of the central orientation, a bootstrap confidence region is available through the `CIradius` function.  There are two versions available: one for the matrix representation and one for the quaternion representation.  The `fun` option allows for the choice of `mean` or `median` and all the options that go along with them, e.g., `type=intrinsic` or `projected`.  One can also set the number of bootstrap sample to draw (`B`), the size of each bootstrap draw (`m`) and the percentile wanted (`q`).

The procedure used is as follows:

0. Estimate $\widehat{S}$ from $(R_1,\dots,R_n)$ in Rs
1. Sample $R_1$,...,$R_m$ from Rs with replacement
2. Estimate $\widehat{S}^*$ from bootstrap sample
3. Compute $\hat{T}=d_r(\widehat{S},\widehat{S}^*)$
4. Repeat steps 1-3 B times
5. Report q% percentile of $\hat{T}$ to be the radius of the confidence 'cone'

In place of $d_r$, Bingham et. al. (2009b) used the maximum absolute angle between each axis of $\widehat{S}$ and $\widehat{S}^*$, call this $\alpha$.  Our method is slightly more conservative because $d_r>\alpha$.

## Visualization
Based on `ggplot2` by Wickham (2009) we have developed a method for visualizing random rotations in SO(3) in three-dimensions.  The function `plot` plots the data on the sphere along with all four estimators of the central orientation if the option `show_estimates` is set to `TRUE`.  By default the data are centered around  the three-dimensional identity matrix, but any matrix in SO(3) can be used as the center with the `center` option.  Finally, since only one column of the three is displayed at a time the `column` option allows the user to choose which column to display.  Additional options can be passed to the call to `qplot` too.