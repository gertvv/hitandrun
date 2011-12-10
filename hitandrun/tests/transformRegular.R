## Check that transformations are applied correctly in non-homogeneous case.

library(hitandrun)

# constraints: x_1 >= 0, x_2 >= 0, x_1 + x_2 <= 1
A <- rbind(c(-1, 0), c(0, -1), c(1, 1))
b <- c(0, 0, 1)
d <- c("<=", "<=", "<=")
constr <- list(constr=A, rhs=b, dir=d)

# take a point x0 within the polytope
x0 <- c(0.25, 0.25)

# transformation that flips x_1
transform <- rbind(c(-1, 0), c(0, 1)) 

# sample 10,000 points
samples <- har(x0, constr, 1E4, transform=transform)$samples

stopifnot(dim(samples) == c(1E4, 2))

# Check that x_1 <= 0
stopifnot(samples[,1] <= 0)
# Check that x_2 >= 0
stopifnot(samples[,2] >= 0)

# Check that -x_1 + x_2 <= 1
stopifnot(-samples[,1] + samples[,2] <= 1)
