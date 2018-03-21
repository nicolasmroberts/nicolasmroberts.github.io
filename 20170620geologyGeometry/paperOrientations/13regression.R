


### DATA ###

source("library/all.R")

# Load the Engebretson et al. (1984) data as cited by Prentice (1987). Times are in millions of years before present.
from2PsiThetaPhiDeg <- function(twoPsiThetaPhiDeg) {
  rho <- sin(twoPsiThetaPhiDeg[1] * degree / 2)
  phi <- twoPsiThetaPhiDeg[2] * degree
  theta <- twoPsiThetaPhiDeg[3] * degree
  v <- cartesianFromSpherical(c(rho, phi, theta))
  q0 <- cos(twoPsiThetaPhiDeg[1] * degree / 2)
  rotMatrixFromQuaternion(c(q0, v[3], v[1], v[2]))
}
farpacTs <- c(0, 5, 9, 17, 28, 37, 43, 48, 56, 61, 66, 71, 85, 119, 127, 135, 145, 165)
farpac2PsiThetaPhiDegs <- list(c(0, 0, 0), c(3.1, -19, 207), c(6.3, -69, 186), c(12.6, -79, 188), c(20, -81, 205), c(32.6, -88, 155), c(43.5, -87, 85), c(51.4, -86, 70), c(61.7, -85, 49), c(65, -86, 32), c(68.1, -86, 352), c(71.3, -86, 13), c(80.3, -85, 315), c(107.4, -78, 285), c(108.8, -75, 275), c(110.5, -71, 269), c(114.2, -68, 269), c(121.7, -63, 270))
farpacRs <- lapply(farpac2PsiThetaPhiDegs, from2PsiThetaPhiDeg)

# This is Fig. 13 without any regression on it.
rotEqualAnglePlot(points=farpacRs, colors=hues(farpacTs))



### GEODESIC REGRESSION ###

# Do geodesic regression. Check that err is 0 (if not, increase steps) and the minimum eigenvalue is positive.
farpacGeod <- rotGeodesicRegression(farpacTs, farpacRs, numSteps=1000)
farpacGeod$error
farpacGeod$minEigenvalue

# Check the quality of the fit, both in the R^2 statistic and graphically. When I did this, I got R^2 = 0.555.
farpacGeod$rSquared
farpacCurve <- lapply(seq(from=min(farpacTs), to=max(farpacTs), length.out=100), farpacGeod$prediction)
rotEqualAnglePlot(points=farpacRs, curves=list(farpacCurve), colors=hues(farpacTs))

# Do a permutation test. This might take several minutes. Check that the errors are 0 and the minimum eigenvalues are positive. Construct a p-value for the original result, based on the fraction of permuted data sets that produce a larger R^2. When I did this, I got p = 0 based on 1000 permutations.
perms <- rotGeodesicRegressionPermutations(farpacTs, farpacRs, numPerms=1000, numSteps=10000)
length(perms)
sum(perms > farpacGeod$rSquared)
sum(perms > farpacGeod$rSquared) / length(perms)

# So in summary we have a highly significant result that explains R^2 = 0.555 of the variance in the data. But the geodesic is obviously not capturing the tendency in the data very well. We need to consider more complicated curves.



### FARALLON-PACIFIC KERNEL REGRESSION ###

# Mathematica gave me 0.07508508. R gave me 0.03976778. That's weird...
h <- 0.07508507972982359
#h <- rotBandwidthForKernelRegression(farpacTs, farpacRs, dnorm)
h <- 0.03976778

# Do the kernel regression. This might take one minute. Check that errors are 0 and the minimum eigenvalues are positive.
farpacKerns <- lapply(seq(from=min(farpacTs), to=max(farpacTs), length.out=100),
                      rotKernelRegression, farpacTs, farpacRs, h, numSteps=1000)
sapply(farpacKerns, function(regr) regr$error)
sapply(farpacKerns, function(regr) {regr$minEigenvalue > 0})

# Plot the regression curve.
farpacKernCurve <- lapply(farpacKerns, function(regr) regr$r)
rotEqualAnglePlot(points=farpacRs, curves=list(farpacKernCurve), colors=hues(farpacTs))

# With R's version of h, R^2 is 0.98, which is great.
farpacKernRSq <- rotRsquaredForKernelRegression(farpacTs, farpacRs, h, numSteps=1000)
farpacKernRSq

# Do a permutation test. (For me, 100 permutations takes about half an hour and 1,000 permutations takes about 5 hours.) Construct a p-value for the original result, based on the fraction of permuted data sets that produce a larger R^2. When I did this, I got p = 0 based on 1,000 permutations.
perms <- rotKernelRegressionPermutations(farpacTs, farpacRs, h, numPerms=100, numSteps=1000)
length(perms)
sum(perms > farpacKernRSq$rSquared)
sum(perms > farpacKernRSq$rSquared) / length(perms)

geoFileFromData(fileName="farPacRots18.csv", dataFrame=data.frame(t(sapply(farpacRs, as.numeric))))


