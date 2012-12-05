theta <- acos(runif(1, -1, 1))
phi <- runif(1, -pi, pi)
u<- c(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta))
r<-rvmises(1)

R<-SO3(u,r)
Q<-Q4(u,r)
EAi<-EA(u,r)

context("Basics")
#Does the angle function extract the correct angle from the rotation
expect_equal(angle(R),abs(r))
expect_equal(angle(Q),abs(r))
expect_equal(angle(EAi),abs(r))

#Does the axis function extract the correct axis from the rotation
expect_equal(abs(axis2(R)),abs(u))
expect_equal(abs(axis2(Q)),abs(u))
expect_equal(abs(axis2(EAi)),abs(u))

#Can we recreate the rotation matrix from the axis and angle functions
expect_equal(SO3(axis2(R),angle(R)),as.SO3(matrix(R,3,3)))
expect_equal(Q4(axis2(Q),angle(Q)),Q)
expect_equal(EA(axis2(EAi),angle(EAi)),EAi)