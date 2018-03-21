


# Copyright 2016-2017 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# The complete 3D orientation of an object can be described using three perpendicular unit vectors in right-handed order. And actually the third vector is redundant, as it is simply the cross product of the other two. For example, the orientation of a duck in flight could be described by one vector pointing from its heart toward its head, and a second vector pointing from its heart out of its back. The particular convention that you use doesn't matter, as long as you use it consistently. You put these three vectors into the rows of a matrix R, which is then special orthogonal. This matrix is the rotation that rotates the vectors back to the x-, y-, and z-coordinate axes. In this file, 'rotation matrix' means 'real 3x3 matrix that is special orthogonal'.



### SYMMETRY GROUPS ###

# Frequently orientations are subject to some symmetry group G, which is a finite set of rotations (satisfying certain properties). For any rotation Q in G, the rotation matrix Q %*% R represents the same orientation as R does. Here are some common symmetry groups.

#' When trivial symmetry is used, the orientations are simply rotations. This case is so important that we have separate code, in rotations.R, for doing it. But let's include the trivial group, for completeness.
oriTrivialGroup <- list(diag(c(1, 1, 1)))

#' Ray-in-plane symmetry is applicable to faults-with-slip-directions, such as slickensides (and certain minerals).
oriRayInPlaneGroup <- list(
  diag(c(1, 1, 1)),
  diag(c(-1, 1, -1)))

#' Line-in-plane symmetry is applicable to foliation-lineations, cylindrical fold orientations, triaxial ellipsoid orientations, and earthquake focal mechanisms (and olivine).
oriLineInPlaneGroup <- list(
  diag(c(1, 1, 1)),
  diag(c(1, -1, -1)),
  diag(c(-1, 1, -1)),
  diag(c(-1, -1, 1)))

#' Trigonal trapezohedral is the point group of alpha-quartz.
oriTrigonalTrapezohedralGroup <- list(
  diag(c(1, 1, 1)),
  rotMatrixAboutZ(pi * 2 / 3),
  rotMatrixAboutZ(pi * 4 / 3), 
  diag(c(1, -1, -1)),
  rotMatrixAboutZ(pi * 2 / 3) %*% diag(c(1, -1, -1)) %*% t(rotMatrixAboutZ(pi * 2 / 3)),
  rotMatrixAboutZ(pi * 4 / 3) %*% diag(c(1, -1, -1)) %*% t(rotMatrixAboutZ(pi * 4 / 3)))

#' Hexagonal trapezohedral is the point group of beta-quartz.
oriHexagonalTrapezohedralGroup <- list(
  diag(c(1, 1, 1)),
  rotMatrixAboutZ(pi * 1 / 3),
  rotMatrixAboutZ(pi * 2 / 3),
  diag(c(-1, -1, 1)), 
  rotMatrixAboutZ(pi * 4 / 3), 
  rotMatrixAboutZ(pi * 5 / 3), 
  diag(c(1, -1, -1)),
  rotMatrixAboutZ(pi * 1 / 3) %*% diag(c(1, -1, -1)) %*% t(rotMatrixAboutZ(pi * 1 / 3)),
  rotMatrixAboutZ(pi * 2 / 3) %*% diag(c(1, -1, -1)) %*% t(rotMatrixAboutZ(pi * 2 / 3)),
  diag(c(-1, -1, 1)) %*% diag(c(1, -1, -1)) %*% diag(c(-1, -1, 1)),
  rotMatrixAboutZ(pi * 4 / 3) %*% diag(c(1, -1, -1)) %*% t(rotMatrixAboutZ(pi * 4 / 3)),
  rotMatrixAboutZ(pi * 5 / 3) %*% diag(c(1, -1, -1)) %*% t(rotMatrixAboutZ(pi * 5 / 3)))

#' Replicating a list of orientations into their multiple rotation representatives.
#' 
#' @param rs A list of rotation matrices. The representative rotations Rs.
#' @param group A list of rotation matrices. The symmetry group G.
#' @return A list of rotation matrices. The set G Rs. For each orientation G R, its |G| representatives are spread out in this list.
oriSymmetrizedRotations <- function(rs, group) {
  unlist(lapply(group, function(g) lapply(rs, function(r) g %*% r)), recursive=FALSE, use.names=FALSE)
}



### MISCELLANEOUS METHODS ###

#' The distance between two orientations as points in SO(3) / G.
#'
#' @param r A rotation matrix.
#' @param q A rotation matrix.
#' @param group A list of rotation matrices. The symmetry group G.
#' @return A real number (in the interval [0, pi]). The distance from G R to G R.
oriDistance <- function(r, q, group) {
  qRT <- tcrossprod(q, r)
  trGQRT <- max(sapply(group, function(g) tr(g %*% qRT)))
	arcCos((trGQRT - 1) / 2)
}

#' Diameter of a set of orientations.
#' 
#' @param rs A list of rotation matrices.
#' @param group A list of rotation matrices. The symmetry group G.
#' @return A real number (non-negative).
oriDiameter <- function(rs, group) {
  f <- function(i) {
    max(sapply(1:(i - 1), function(j) oriDistance(rs[[i]], rs[[j]], group)))
  }
  max(sapply(2:(length(rs)), f))
}

#' Selecting an orientation representative near a given rotation.
#' 
#' @param r A rotation matrix.
#' @param center A rotation matrix.
#' @param group A list of rotation matrices. The symmetry group G.
#' @return A rotation matrix. The element of G R that is closest to center.
oriNearestRepresentative <- function(r, center, group) {
  gRs <- lapply(group, function(g) g %*% r)
  trCTGRs <- sapply(gRs, function(gr) tr(crossprod(center, gr)))
  i <- which.max(trCTGRs)
  gRs[[i]]
}

#' Selecting orientation representatives near a given rotation.
#' 
#' @param rs A list of rotation matrices.
#' @param center A rotation matrix.
#' @param group A list of rotation matrices. The symmetry group G.
#' @return A list of rotation matrices. For each R, the element of G R that is closest to center.
oriNearestRepresentatives <- function(rs, center=rs[[1]], group) {
  lapply(rs, oriNearestRepresentative, center, group)
}

#' The Frechet (geodesic L^2) variance of a set of orientations.
#'
#' @param rs A list of rotation matrices.
#' @param center A rotation matrix. Often the mean of the rs.
#' @param group A list of rotation matrices. The symmetry group G.
#' @return A real number (non-negative). The variance of the points G R about the point G center.
oriVariance <- function(rs, center, group) {
  dists <- sapply(rs, oriDistance, center, group)
  sum(dists^2) / (2 * length(rs))
}

#' The Frechet (geodesic L^2) mean of a set of orientations as points in SO(3) / G.
#'
#' An iterative algorithm for computing the Frechet mean --- the orientation that minimizes the Frechet variance. The iterations continue until change-squared of epsilon is achieved or numSteps iterations have been used.
#' @param rs A list of rotation matrices.
#' @param group A list of rotation matrices. The symmetry group G.
#' @param numSeeds A real number (positive integer). How many seeds to try, in trying to find a global optimum.
#' @param numSteps A real number (positive integer). Bound on how many iterations to use.
#' @return A list consisting of $mean (a special orthogonal real 3x3 matrix), $variance (a real number), $changeSquared (a real number), and $numSteps (a non-negative integer). changeSquared is the square of the size of the final step. numSteps is the number of iterations used.
oriMeanVariance <- function(rs, group, numSeeds=1, numSteps=1000) {
  seeds <- sample(rs, numSeeds)
  # No variance is ever as large as 5.
  best <- c(5)
  for (seed in seeds) {
    rBar <- seed
    changeSquared <- epsilon + 1.0
    k <- 0
    while (changeSquared >= epsilon && k < numSteps) {
      w <- diag(c(0, 0, 0))
      for (r in rs) {
        rBarTGRs <- lapply(group, function(g) crossprod(rBar, g %*% r))
        i <- which.max(sapply(rBarTGRs, tr))
        w <- w + rotLog(rBarTGRs[[i]])
      }
      w <- w / length(rs)
      rBar <- rBar %*% rotExp(w)
      changeSquared <- tr(crossprod(w, w))
      k <- k + 1
    }
    var <- oriVariance(rs, rBar, group)
    if (var < best[[1]])
      best <- list(var, rBar, changeSquared, k)
  }
  list(variance=best[[1]], mean=best[[2]], changeSquared=best[[3]], numSteps=best[[4]])
}

#' The Frechet (geodesic L^2) mean. Convenience shortcut for oriMeanVariance.
#' 
#' @param rs A list of rotation matrices.
#' @param group A list of rotation matrices. The symmetry group G.
#' @param numSeeds A real number (positive integer). How many seeds to try, in trying to find a global optimum.
#' @param numSteps A real number (positive integer). Bound on how many iterations to use.
#' @return A rotation matrix. The Frechet mean.
oriFrechetMean <- function(rs, group, numSeeds=1, numSteps=1000) {
  oriMeanVariance(rs, group, numSeeds=numSeeds, numSteps=numSteps)$mean
}

#' Projected arithmetic mean of a set of orientations as points in SO(3) / G.
#'
#' An iterative algorithm for computing the projected arithmetic mean of orientations rather than rotations. See Bachmann et al. (2010). The iterations continue until the change is small or the allowed number of iterations has been exhausted.
#' @param rs A list of rotation matrices.
#' @param group A list of rotation matrices. The symmetry group G.
#' @param numSeeds A real number (positive integer). The number of seeds to try, in trying to find a global optimum.
#' @param numSteps A real number (positive integer). Bound on how many iterations to use.
#' @return A list consisting of $values (as in rotMeanScatter), $rotations (as in rotMeanScatter), $changeSquared (a real number), and $numSteps (a non-negative integer). changeSquared is the square of the size of the final step. numSteps is the number of iterations used.
oriMeanScatter <- function(rs, group, numSeeds=1, numSteps=1000) {
  seeds <- sample(rs, numSeeds)
  best <- c(-2)
  for (seed in seeds) {
    rBar <- seed
    # No concentration is ever bigger than 1.
    concenOld <- -6
    concenNew <- -4
    k <- 0
    while ((concenNew - concenOld)^2 >= epsilon && k < numSteps) {
      rsNear <- oriNearestRepresentatives(rs, rBar, group)
      meanScatter <- rotMeanScatter(rsNear)
      concenOld <- concenNew
      concenNew <- meanScatter$values[[1]]
      rBar <- meanScatter$rotations[[1]]
      k <- k + 1
    }
    if (concenNew > best[[1]])
      best <- list(concenNew, meanScatter, (concenNew - concenOld)^2, k)
  }
  list(values=best[[2]]$values, rotations=best[[2]]$rotations, changeSquared=best[[3]], numSteps=best[[4]])
}

#' Projected mean orientation. Convenience shortcut for oriMeanScatter.
#' 
#' @param rs A list of rotation matrices.
#' @param group A list of rotation matrices. The symmetry group G.
#' @param numSeeds A real number (positive integer). The number of seeds to try, in trying to find a global optimum.
#' @param numSteps A real number (positive integer). Bound on how many iterations to use.
#' @return A rotation matrix.
oriProjectedMean <- function(rs, group, numSeeds=1, numSteps=1000) {
  oriMeanScatter(rs, group, numSeeds, numSteps)$rotations[[1]]
}



### TANGENT SPACE METHODS ###

# The tangent space methods all operate in the tangent space at a point. To apply them to orientations, simply pick the representative for each orientation that is closest to the point, using oriNearestRepresentative(s). As always, tangent space methods are best for data that are tightly concentrated.

# When working modulo a symmetry group, there is an additional upper bound on how spread out the data can be. Let r be the minimum distance between any two group elements. In other words, r is the distance from the identity to the nearest other group element. See rotSeparation. Let d be the diameter of the data set (oriDiameter). Then we can safely work with representatives as long as d < r / 2. For then oriNearestRepresentative never leaves the chosen set of representatives.



### BOOTSTRAPPING ###

#' Inference for mean orientations, based on non-parametric bootstrapping.
#' 
#' This function bootstraps the Frechet orientation mean, returning various pieces of information. The first is a list of the bootstrapped means. The user should oriEqualAnglePlot and oriEqualVolumePlot this list, to make sure that the means form a fairly tight ellipsoidal cluster. If so, then the second piece of information may be used: An R function that takes an orientation R0 as input, as produces as output a p-value for the hypothesis that the mean of the population is R0. This p-value is the fraction of means that are farther from the mean of means than R0 is, based on the Mahalanobis distance of the means.
#' @param rs A list of rotation matrices. The data.
#' @param numBoots A real number (positive integer). The number of bootstrap samples.
#' @param group A list of rotation matrices. The symmetry group G.
#' @return A list consisting of $bootstraps (a list of rotation matrices), $pvalueOri (an R function from {orientation matrices} to {real numbers}), and everything returned by rotMahalanobisInference.
oriBootstrapInference <- function(rs, numBoots, group) {
  boots <- replicate(numBoots, oriFrechetMean(sample(rs, length(rs), replace=TRUE), group), simplify=FALSE)
  bootMean <- oriFrechetMean(boots, group, numSeeds=10, numSteps=10000)
  boots <- oriNearestRepresentatives(boots, bootMean, group)
  infer <- rotMahalanobisInference(boots, bootMean)
  pfunc <- function(r) {
    max(sapply(group, function(g) infer$pvalue(g %*% r)))
  }
  infer$pvalueOri <- pfunc
  infer$bootstraps <- boots
  infer
}



### REGRESSION ###

# Like oriGeodesicRegression, but does not rescale its x-data ahead of time.
oriNativeGeodesicRegression <- function(xs, rs, group, numSteps=100) {
  # Let Q be the R whose x is closest to zero.
  q <- rs[[which.min(as.numeric(xs)^2)]]
  # Define the function to be minimized.
  n <- length(rs)
  e <- function(mw) {
    b <- rotExp(rotAntisymmetricFromVector(mw[4:6])) %*% q
    m <- rotAntisymmetricFromVector(mw[1:3])
    f <- function(i) {
      oriDistance(rs[[i]], rotExp(xs[[i]] * m) %*% b, group)^2
    }
    sum(sapply(1:n, f)) / (2 * n)
  }
  # Find the minimum, using the constant geodesic Q as the seed.
  seed <- c(0, 0, 0, 0, 0, 0)
  solution <- optim(seed, e, hessian=TRUE, control=list(maxit=numSteps))
  # Report diagnostic information.
  eigvals <- eigen(solution$hessian, symmetric=TRUE, only.values=TRUE)$values
  m <- rotAntisymmetricFromVector(solution$par[1:3])
  b <- rotExp(rotAntisymmetricFromVector(solution$par[4:6])) %*% q
  rBar <- oriFrechetMean(rs, group)
  rSq <- 1 - solution$value / oriVariance(rs, rBar, group)
  prediction <- function(x) {
    rotExp(x * m) %*% b
  }
  list(m=m, b=b, error=solution$convergence, minEigenvalue=min(eigvals), rSquared=rSq, prediction=prediction)
}

#' Best-fit geodesic curve in SO(3) / G.
#'
#' Returns the best-fit geodesic R(x) = (exp (x m)) b for the given (x, R) pairs. Minimizes the sum of squared distances in the R-direction only, as in ordinary linear regression --- not in the x-direction.
#' @param xs A list of real numbers.
#' @param rs A list of rotation matrices of same length as xs.
#' @param group A list of rotation matrices. The symmetry group G.
#' @param numSteps A real number (positive integer). Limits the number of iterations to use.
#' @return A list consisting of $m (anti-symmetric 3x3 real matrix), $b (rotation matrix), $error (real number), $minEigenvalue (real number), $rSquared (real number), $prediction (R function). error is 0 if and only if the minimization succeeds. If error == 1, try increasing numSteps, to 1000 say. minEigenvalue is the least eigenvalue of the Hessian at the solution; worry about the result if it is non-positive 0. rSquared is the R^2 statistic measuring the amount of variance captured, from none (0) to all (1). prediction takes a real number x as input, and returns the predicted rotation matrix R(x) as output.
oriGeodesicRegression <- function(xs, rs, group, numSteps=100) {
  x0 <- min(xs)
  x1 <- max(xs)
  regr <- oriNativeGeodesicRegression(scales(xs), rs, group, numSteps)
  mNew <- regr$m / (x1 - x0)
  bNew <- rotExp(regr$m * -x0 / (x1 - x0)) %*% regr$b
  prediction <- function(x) {
    regr$prediction((x - x0) / (x1 - x0))
  }
  list(m=mNew, b=bNew, error=regr$error, minEigenvalue=regr$minEigenvalue, rSquared=regr$rSquared, prediction=prediction)
}

#' Permutation test for significance of a geodesic regression.
#' 
#' # Returns up to numPerms R^2 values. May be fewer than numPerms, because of some regressions failing. Let n be the dimension of this vector and g the number of R^2 values greater than the R^2 for the original regression of the data. Let p = g / n. Small values of p indicate that the dependency detected by the regression is meaningful. Uses an internal rescaling, just like oriGeodesicRegression.
#' @param xs A vector of real numbers.
#' @param rs A list of rotation matrices, of same length as xs.
#' @param numPerms A real number (positive integer). The number of permutations to try.
#' @param group A list of rotation matrices. The symmetry group G.
#' @param numSteps A real number (positive integer). Limits the number of iterations to use.
#' @return A vector of real numbers, of dimension at most numPerms.
oriGeodesicRegressionPermutations <- function(xs, rs, numPerms, group, numSteps=10000) {
  ys <- scales(xs)
  f <- function(i) {
    print(i / numPerms)
    regr <- oriNativeGeodesicRegression(sample(ys, size=length(ys)), rs, group, numSteps)
    c(regr$error, regr$minEigenvalue, regr$rSquared)
  }
  perms <- sapply(1:numPerms, f)
  perms[3,][perms[1,] == 0 & perms[2,] > 0]
}

#' Kernel regression to fit an orientation R(x) to (x, R) data.
#' 
#' This function interpolates/extrapolates an orientation R for a given x-value, based on a given set of (x, R) data. The chosen bandwidth h may have a substantial effect on the results. See oriBandwidthForKernelRegression.
#' @param x A real number. The x-value at which to predict the orientation R.
#' @param xs A vector of real numbers. The x-values at which R(x) is known.
#' @param rs A list of rotation matrices. The values of R corresponding to xs.
#' @param h A real number (positive). The bandwidth. The kernel k is scaled by h, meaning that k is changed to function(x) {k(x / h) / h}.
#' @param group A list of rotation matrices. The symmetry group G.
#' @param k An R function from {real numbers} to {real numbers}. The kernel, not scaled by bandwidth.
#' @param numSteps A real number (positive integer). Limits the number of iterations to use.
#' @return A list consisting of $r (the rotation matrix R(x)), $error (which should be 0; if it is 1, then increase numSteps), and $minEigenvalue (which should be positive; if not, then the minimization has failed).
oriKernelRegression <- function(x, xs, rs, h, group, k=dnorm, numSteps=100) {
  # Let Q be the Ri whose xi is closest to x.
  q <- rs[[which.min((as.numeric(xs) - x)^2)]]
  # Define the function to be minimized.
  kh <- function(x) {k(x / h) / h}
  e <- function(w) {
    r <- rotExp(rotAntisymmetricFromVector(w)) %*% q
    f <- function(i) {
      kh(x - xs[[i]]) * oriDistance(rs[[i]], r, group)^2
    }
    sum(sapply(1:length(xs), f)) / sum(sapply(xs, function(y) kh(x - y)))
  }
  # Minimize the function with respect to R, using Q as the seed.
  seed <- c(0.0, 0.0, 0.0)
  solution <- optim(seed, e, hessian=TRUE, control=list(maxit=numSteps))
  # Report diagnostic information.
  eigvals <- eigen(solution$hessian, symmetric=TRUE, only.values=TRUE)$values
  r <- rotExp(rotAntisymmetricFromVector(solution$par)) %*% q
  list(r=r, error=solution$convergence, minEigenvalue=min(eigvals))
}

#' Cross-validation algorithm for choosing the bandwidth for kernel regression.
#'
#' @param xs A vector of real numbers. The x-values at which R(x) is known. Assumed to have length >= 3.
#' @param rs A list of rotation matrices. The values of R corresponding to xs.
#' @param group A list of rotation matrices. The symmetry group G.
#' @param k An R function from {real numbers} to {real numbers}. The kernel, not scaled by bandwidth.
#' @param numSteps A real number (positive integer). Limits the number of iterations to use.
#' @return A real number (positive). The bandwidth to use in kernel regression.
oriBandwidthForKernelRegression <- function(xs, rs, group, k=dnorm, numSteps=100) {
    # Build jackknifed lists ahead of time.
    n <- length(xs)
    xjs <- listOmitting(xs)
    rjs <- listOmitting(rs)
    # Define the function to minimize.
    g <- function(h) {
        rjhs <- lapply(1:n, function(j) oriKernelRegression(xs[[j]], xjs[[j]], rjs[[j]], h=h, group=group, k=k, numSteps=numSteps))
        sum(sapply(1:n, function(j) oriDistance(rs[[j]], rjhs[[j]]$r, group)^2))
    }
    # Minimize the function with respect to h.
    solution <- optimize(g, interval=c(0, pi))
    solution$minimum
}

#' R^2 statistic for kernel regression.
#'
#' @param xs A vector of real numbers. The x-values at which R(x) is known.
#' @param rs A list of rotation matrices. The values of R corresponding to xs.
#' @param h A real number (positive). The bandwidth. The kernel k is scaled by h, meaning that k is changed to function(x) {k(x / h) / h}.
#' @param group A list of rotation matrices. The symmetry group G.
#' @param k An R function from {real numbers} to {real numbers}. The kernel, not scaled by bandwidth.
#' @param numSteps A real number (positive integer). Limits the number of iterations to use.
#' @return A real number (usually between 0 and 1). The R^2 statistic.
oriRsquaredForKernelRegression <- function(xs, rs, h, group, k=dnorm, numSteps=100) {
    n <- length(xs)
    qs <- lapply(xs, function(x) oriKernelRegression(x, xs, rs, h=h, group=group, k=k, numSteps=numSteps)$r)
    e <- sum(sapply(1:n, function(i) oriDistance(rs[[i]], qs[[i]], group)^2)) / (2 * n)
    1 - e / oriVariance(rs, oriFrechetMean(rs, group), group)
}


### PLOTTING ###

#' X-Z-X Euler angle plot of orientations as symmetric sets of rotations.
#'
#' This function is much like rotEulerAnglePlot, but it plots |G| symmetric copies of the points, curves, and triangles. This plot does not have an equal-angle or equal-volume property. Indeed, it is extremely distorted near certain gimbal lock lines. Curves and triangles are not supported.
#' @param points A list of rotation matrices.
#' @param group A list of rotation matrices. The symmetry group G.
#' @param showBoundary Logical. Whether to show the bounding box.
#' @param ... Plotting options to be passed to plot3D. If colors is used, its length should equal that of points.
#' @return NULL.
oriEulerAnglePlot <- function(points, group, showBoundary=FALSE, ...) {
  rotEulerAnglePlot(points=oriSymmetrizedRotations(points, group), showBoundary=showBoundary, ...)
}

#' Axis-angle plot of orientations as symmetric sets of rotations.
#' 
#' This function is much like rotAxisAnglePlot, but it plots |G| symmetric copies of the points, curves, and triangles. Curves that cross the boundary are automatically clipped. On the other hand, triangles are not clipped. Thus you probably don't want to use any triangles.
#' @param group A list of rotation matrices. The symmetry group G.
#' @param points A list of rotation matrices.
#' @param curves A list of lists of rotation matrices.
#' @param triangles A list of length-3 lists of rotation matrices.
#' @param ... Plotting options to be passed to rotBallPlot. If colors is used, its length should equal that of points.
#' @return NULL.
oriAxisAnglePlot <- function(group, points=list(), curves=list(), triangles=list(), ...) {
  pointss <- oriSymmetrizedRotations(points, group)
  f <- function(g, ps) {lapply(ps, function(p) g %*% p)}
  curvess <- Reduce(c, lapply(group, function(g) lapply(curves, function(ps) f(g, ps))))
  triangless <- Reduce(c, lapply(group, function(g) lapply(triangles, function(ps) f(g, ps))))
  rotAxisAnglePlot(pointss, curvess, triangless, ...)
}

#' Equal-angle plot of orientations as symmetric sets of rotations.
#'
#' This function is much like rotEqualAnglePlot, but it plots |G| symmetric copies of the points, curves, and triangles.
#' @param group A list of rotation matrices. The symmetry group G.
#' @param points A list of rotation matrices.
#' @param curves A list of lists of rotation matrices.
#' @param triangles A list of length-3 lists of rotation matrices.
#' @param ... Plotting options to be passed to rotBallPlot. If colors is used, its length should equal that of points.
#' @return NULL.
oriEqualAnglePlot <- function(group, points=list(), curves=list(), triangles=list(), ...) {
  pointss <- oriSymmetrizedRotations(points, group)
  f <- function(g, ps) {lapply(ps, function(p) g %*% p)}
  curvess <- Reduce(c, lapply(group, function(g) lapply(curves, function(ps) f(g, ps))))
  triangless <- Reduce(c, lapply(group, function(g) lapply(triangles, function(ps) f(g, ps))))
  rotEqualAnglePlot(pointss, curvess, triangless, ...)
}

#' Equal-volume plot of orientations as symmetric sets of rotations.
#' 
#' This function is much like rotEqualVolumePlot, but it plots |G| symmetric copies of the points, curves, and triangles. Curves that cross the boundary are automatically clipped. On the other hand, triangles are not clipped. Thus you probably don't want to use any triangles.
#' @param group A list of rotation matrices. The symmetry group G.
#' @param points A list of rotation matrices.
#' @param curves A list of lists of rotation matrices.
#' @param triangles A list of length-3 lists of rotation matrices.
#' @param ... Plotting options to be passed to rotBallPlot. If colors is used, its length should equal that of points.
#' @return NULL.
oriEqualVolumePlot <- function(group, points=list(), curves=list(), triangles=list(), ...) {
  pointss <- oriSymmetrizedRotations(points, group)
  f <- function(g, ps) {lapply(ps, function(p) g %*% p)}
  curvess <- Reduce(c, lapply(group, function(g) lapply(curves, function(ps) f(g, ps))))
  triangless <- Reduce(c, lapply(group, function(g) lapply(triangles, function(ps) f(g, ps))))
  rotEqualVolumePlot(pointss, curvess, triangless, ...)
}

#' Rodrigues plot of orientations as symmetric sets of rotations.
#'
#' This function is much like rotRodriguesPlot, but it plots |G| symmetric copies of the points, curves, and triangles.
#' @param group A list of rotation matrices. The symmetry group G.
#' @param points A list of rotation matrices.
#' @param curves A list of lists of rotation matrices.
#' @param triangles A list of length-3 lists of rotation matrices.
#' @param ... Plotting options to be passed to rotBallPlot. If colors is used, its length should equal that of points.
#' @return NULL.
oriRodriguesPlot <- function(group, points=list(), curves=list(), triangles=list(), ...) {
  pointss <- oriSymmetrizedRotations(points, group)
  f <- function(g, ps) {lapply(ps, function(p) g %*% p)}
  curvess <- Reduce(c, lapply(group, function(g) lapply(curves, function(ps) f(g, ps))))
  triangless <- Reduce(c, lapply(group, function(g) lapply(triangles, function(ps) f(g, ps))))
  rotRodriguesPlot(pointss, curvess, triangless, ...)
}



### HIGHER-LEVER PLOTS ###

#' Equal-angle plot of geodesic regression results.
#' 
#' The user may choose to draw extra spurs, joining the Rs to their corresponding predictions on the curve.
#' @param xs A vector of real numbers.
#' @param rs A list of rotation matrices.
#' @param regr The output of oriGeodesicRegression with those xs and rs.
#' @param group A list of rotation matrices. The symmetry group G.
#' @param colors A list of strings (colors). Used to color the points only.
#' @param spurs A list of rotation matrices. When not NULL, it's usually equal to Rs.
#' @param ... Other plotting options to be passed to rotEqualAnglePlot.
#' @return NULL.
oriGeodesicRegressionPlot <- function(xs, rs, regr, group, colors=hues(xs), spurs=NULL, ...) {
  rotA <- rotExp(min(xs) * regr$m) %*% regr$b
  rotB <- rotExp(mean(range(xs)) * regr$m) %*% regr$b
  rotC <- rotExp(max(xs) * regr$m) %*% regr$b
  curves <- list(rotGeodesicPoints(rotA, rotB, 20), rotGeodesicPoints(rotB, rotC, 20))
  if (!is.null(spurs)) {
    spurCurves <- lapply(1:length(xs), function(i) {rotGeodesicPoints(spurs[[i]], regr$prediction(xs[[i]]), 10)})
    curves <- c(curves, spurCurves)
  }
  oriEqualAnglePlot(group=group, points=rs, curves=curves, colors=colors, ...)
}


