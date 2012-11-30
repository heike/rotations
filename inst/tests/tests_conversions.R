#library(testthat)
library(rotations)

u<-c(1/sqrt(2),0,1/sqrt(2))
r<-pi/2	

context("Conversions")
expect_equal(Q4(SO3(u,r)),Q4(u,r))
expect_equal(EA(SO3(u,r)),EA(u,r))
expect_equal(Q4(EA(u,r)),Q4(u,r))
expect_equal(EA(Q4(u,r)),EA(u,r))
expect_equal(SO3(EA(u,r)),SO3(u,r))
expect_equal(SO3(Q4(u,r)),SO3(u,r))
