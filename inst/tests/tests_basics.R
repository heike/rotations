#library(rotations)

r<-rvmises(1)
R<-genR(r)
Q<-genR(r,space="Q4")
EA<-genR(r,space="EA")

expect_equal(angle(R),abs(r))
expect_equal(angle(Q),abs(r))
expect_equal(angle(EA),abs(r))
