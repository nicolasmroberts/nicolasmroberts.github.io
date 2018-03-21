


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# A ray is a unit 3D vector, expressed in Cartesian coordinates as an R vector of the form u = c(x, y, z), where x^2 + y^2 + z^2 == 1. A ray is tantamount to a line with a preferred direction along that line. This file offers R functions for dealing with rays. Some are preceded by detailed documentation. The others are not intended for use by end-users.



### MISCELLANY ###

#' Vector normalization (to have length 1).
#' 
#' @param v A d-dimensional vector. Cannot be of length 0.
#' @return A d-dimensional vector of length 1.
rayNormalized <- function(v) {
  v / sqrt(sum(v^2))
}

#' Projection of a vector onto the plane perpendicular to another vector.
#' 
#' For example, when v = c(0, 0, 1), returns the upward-most ray perpendicular to pole.
#' @param pole A d-dimensional vector perpendicular to the plane of interest. Must have non-zero length. Need not be unit.
#' @param v A d-dimensional vector. Should not be parallel to pole. Need not be unit.
#' @return A unit d-dimensional vector.
rayOrthogonalProjection <- function(pole, v) {
  rayNormalized(v - pole * dot(v, pole) / dot(pole, pole))
}

#' A ray perpendicular to a given vector, deterministically.
#' 
#' The result is deterministic but arbitrary. Due to the hairy ball theorem, the result cannot be chosen to depend smoothly on the input.
#' @param v A 3-dimensional vector. Need not be unit.
#' @return A ray, perpendicular to v.
rayOrthogonal <- function(v) {
  if (abs(v[1]) < 0.5)
    w <- c(1, 0, 0)
  else if (abs(v[2]) < 0.5)
    w <- c(0, 1, 0)
  else
    w <- c(0, 0, 1)
  rayOrthogonalProjection(v, w)
}

#' A ray perpendicular to a given vector, probabilistically.
#' 
#' In theory, the returned ray is uniformly chosen on the circle's worth of rays perpendicular to the given vector.
#' @param v A 3-dimensional vector. Need not be unit.
#' @return A ray, perpendicular to v.
rayOrthogonalUniform <- function(v) {
  rayOrthogonalProjection(v, rayUniform())
}

#' The great circle perpendicular to the given ray.
#' 
#' @param pole A ray.
#' @param numSteps The number of line segments to be used in the approximation of the great circle.
#' @return A list of numSteps+1 rays. The first and last one are identical. Otherwise, they are evenly spaced on the circle perpendicular to the given pole, and in order (either clockwise or counter-clockwise, depending on your viewpoint).
rayGreatCircle <- function(pole, numSteps=50) {
  v <- rayOrthogonal(pole)
  w <- cross(pole, v)
  angles <- (0:numSteps) * (2 * pi / numSteps)
  lapply(angles, function(a) {cos(a) * v + sin(a) * w})
}

#' A small circle with the given ray as its pole.
#' 
#' @param v A ray. The pole of the small circle.
#' @param r A real number. The radius of the circle around the pole, in radians, measured along the surface of the sphere. For example, the Arctic circle would be at radius r == 23.5 * degree from the north pole.
#' @return A list of numPoints+1 rays. The first and last one are identical. Otherwise, they are evenly spaced on the small circle, and in order (either clockwise or counter-clockwise, depending on your viewpoint).
raySmallCircle <- function(v, r, numPoints=50) {
  # Make a rotation matrix with v as its third column.
  v2 <- rayOrthogonal(v)
  rotation <- cbind(v2, cross(v, v2), v)
  # I'm not sure what the best way to handle NA values is.
  if (is.na(r))
    r <- 0
  # Make the small circle near the north pole and rotate it to near v.
  lapply(0:numPoints, function(s) as.numeric(rotation %*% cartesianFromSpherical(c(1, r, 2 * pi * s / numPoints))))
}

#' Geodesic between two points on the unit sphere.
#' 
#' Like rayGreatCircle, but specifying the great circle in terms of two points on it, rather than its pole.
#' @param u A ray.
#' @param v A ray.
#' @param numSteps The number of line segments to be used in the approximation of the great circle arc from u to v.
#' @return A list of numSteps+1 rays. The first is u and the last is v. They are evenly spaced and in order.
rayGeodesicPoints <- function(u, v, numSteps=10) {
  dotProduct <- dot(u, v)
  if (dotProduct >= 1)
    replicate(numSteps + 1, u, simplify=FALSE)
  else {
    angle <- arcCos(dotProduct)
    perp <- rayNormalized(v - dotProduct * u)
    lapply(0:numSteps, function(i) {cos(angle * i / numSteps) * u + sin(angle * i / numSteps) * perp})
  }
}

# Returns the ray that is fraction s of the way from u to v.
# Performs badly when u == -v.
rayInterpolation <- function(u, v, s) {
  dotProduct <- dot(u, v)
  if (dotProduct >= 1)
    u
  else {
    angle <- arcCos(dotProduct)
    perp <- rayNormalized(v - dotProduct * u)
    cos(s * angle) * u + sin(s * angle) * perp
  }
}

# Given a nonempty curve in the unit sphere (a list of unit 3D Cartesian vectors).
# Returns a list of curves, where each curve is entirely z >= 0 or entirely z <= 0.
# Whenever it crosses the z == 0 equator, ends and begins adjacent curves there.
rayCurvesUpperLower <- function(us) {
  curves <- list()
  signs <- c()
  curve <- list(us[[1]])
  currentSign <- sign(us[[1]][[3]])
  for (u in us[2:length(us)]) {
    if (sign(u[[3]]) == 0 || sign(u[[3]]) == currentSign)
      # Continue the current curve in the current hemisphere.
      curve[[length(curve) + 1]] <- u
    else if (currentSign == 0) {
      # Continue the current curve and commit to a hemisphere.
      curve[[length(curve) + 1]] <- u
      currentSign <- sign(u[[3]])
    } else {
      v <- curve[[length(curve)]]
      if (sign(v[[3]]) == 0) {
        # The current curve has already ended.
        curves[[length(curves) + 1]] <- curve
        signs[[length(signs) + 1]] <- currentSign
        # Start the new curve at v or -v followed by u.
        if (dot(v, u) > 0)
          curve <- list(v, u)
        else
          curve <- list(-v, u)
        currentSign <- sign(u[[3]])
      } else {
        # Find equatorial point between v and u.
        p <- rayNormalized(v + (v[[3]] / (v[[3]] - u[[3]])) * (u - v))
        # End current curve with whichever of p, -p is closer to v.
        if (dot(v, p) > 0)
          curve[[length(curve) + 1]] <- p
        else
          curve[[length(curve) + 1]] <- -p
        curves[[length(curves) + 1]] <- curve
        signs[[length(signs) + 1]] <- currentSign
        # Begin next curve with whichever is closer to u, followed by u.
        if (dot(u, p) > 0)
          curve <- list(p, u)
        else
          curve <- list(-p, u)
        currentSign <- sign(u[[3]])
      }
    }
  }
  # Finish the current curve and return.
  curves[[length(curves) + 1]] <- curve
  signs[[length(signs) + 1]] <- currentSign
  list(curves=curves, signs=signs)
}

# A list of four triangles that partition the unit sphere.
rayTetrahedron <- list(
  list(c(1, 1, 1) / sqrt(3), c(1, -1, -1) / sqrt(3), c(-1, 1, -1) / sqrt(3)),
  list(c(-1, -1, 1) / sqrt(3), c(1, -1, -1) / sqrt(3), c(-1, 1, -1) / sqrt(3)),
  list(c(-1, 1, -1) / sqrt(3), c(-1, -1, 1) / sqrt(3), c(1, 1, 1) / sqrt(3)),
  list(c(1, -1, -1) / sqrt(3), c(-1, -1, 1) / sqrt(3), c(1, 1, 1) / sqrt(3)))

# Given a triangle of unit vectors in counterclockwise order when viewed from outside.
# Returns a list of 3^numNonAdapt triangles of unit vectors, each in counterclockwise order again.
rayRefinedTriangle <- function(tri, numNonAdapt) {
  if (numNonAdapt == 0)
    list(tri)
  else {
    w1 <- rayNormalized(tri[[1]] + tri[[2]])
    w2 <- rayNormalized(tri[[2]] + tri[[3]])
    w3 <- rayNormalized(tri[[3]] + tri[[1]])
    unlist(list(
      rayRefinedTriangle(list(w1, w2, w3), numNonAdapt - 1),
      rayRefinedTriangle(list(tri[[1]], w1, w3), numNonAdapt - 1),
      rayRefinedTriangle(list(tri[[2]], w2, w1), numNonAdapt - 1),
      rayRefinedTriangle(list(tri[[3]], w3, w2), numNonAdapt - 1)),
      recursive=FALSE, use.names=FALSE)
  }
}

#' Triangular approximation to the unit sphere.
#' 
#' @param numNonAdapt A real number (non-negative integer). The number of times to refine the base triangulation of the sphere.
#' @return A list of triangles, where each triangle is a list of three rays. The number of triangles is 4^(1 + numNonAdapt), so each triangle has spherical area pi / 4^numNonAdapt.
rayTetrahedralSphere <- function(numNonAdapt) {
  unlist(
    lapply(rayTetrahedron, rayRefinedTriangle, numNonAdapt),
    recursive=FALSE, use.names=FALSE) 
}

#' Uniformly random points on the unit sphere.
#' 
#' @param n A real number (positive integer) or NULL.
#' @return If n is NULL, then a single ray. If n is a positive integer, then a list of n rays.
rayUniform <- function(n=NULL) {
  if (is.null(n))
    cartesianFromHorizontal(c(runif(1, 0, 2 * pi), runif(1, -1, 1)))
  else
    replicate(n, rayUniform(), simplify=FALSE)
}



### EXTRINSIC METHODS ###

#' Extrinsic mean and scalar scatter of a set of rays.
#' 
#' Scatter varies between 0 (tight concentration) and 1 (wide dispersion). This scatter is denoted 1 - R-bar in Mardia and Jupp (2000), p. 163. Arguably the preferred measure of scatter should be 2 * (1 - R-bar) or 1 - R-bar^2, but this function implements neither of those.
#' @param us A list of rays.
#' @return A list with elements $mean (ray) and $scatter (a real number between 0 and 1). 
rayMeanScatter <- function(us) {
  u <- arithmeticMean(us)
  r <- sqrt(dot(u, u))
  list(mean=(u / r), scatter=(1 - r))
}

#' Extrinsic mean of a set of rays.
#' 
#' Convenience shortcut for rayMeanScatter(us)$mean.
#' @param us A list of rays.
#' @return A ray.
rayProjectedMean <- function(us) {
  rayNormalized(arithmeticMean(us))
}

#' Bootstrapped extrinsic mean with percentile confidence region and hypothesis tests.
#' 
#' The inference is based on percentiles of Mahalanobis distance in the tangent space at the mean of the bootstrapped means. The user should check that the bootstrapped means form a tight ellipsoidal cluster, before taking such a region seriously.
#' @param ls A list of rays.
#' @param numBoots A real number (positive integer). The number of bootstrapped means to compute. 10,000 might be a good value.
#' @param ... Other arguments to be passed to the underlying rayMahalanobisInference function.
#' @return A list. See rayMahalanobisInference.
rayBootstrapInference <- function(ls, numBoots, ...) {
  boots <- replicate(numBoots, rayProjectedMean(sample(ls, length(ls), replace=TRUE)), simplify=FALSE)
  bootMean <- rayProjectedMean(boots)
  rayMahalanobisInference(boots, bootMean, ...)
}



### DIFFERENTIAL GEOMETRY ###

#' Distance between two points on the unit sphere.
#' 
#' @param u A ray.
#' @param v A ray.
#' @return A real number. The angular difference between the two rays, in radians. Always between 0 and pi, inclusive.
rayDistance <- function(u, v) {
  arcCos(dot(u, v))
}

#' L^2 variance of a set of rays about a point.
#' 
#' I'm not sure about the weighting on this. If you change it, change regression too.
#' @param us A list of rays.
#' @param center A ray. Usually some kind of mean of the us.
#' @return A real number (in the interval [0, pi^2 / 2]).
rayVariance <- function(us, center) {
  sum(sapply(us, function(u) rayDistance(u, center)^2)) / (2 * length(us))
}

#' The exponential map on the unit sphere.
#' 
#' Regards v as a vector in the tangent space of the unit sphere at the point p. Returns a point q on the unit sphere, a distance of |v| from p, in the direction of v. Partial inverse to rayLog.
#' @param p A ray.
#' @param v A 3-dimensional vector, perpendicular to p, not necessarily unit.
#' @return A ray.
rayExp <- function(p, v) {
  normV <- sqrt(dot(v, v))
  if (normV == 0)
    p
  else {
    r <- rbind(p, v / normV, cross(p, v / normV))
    as.numeric(t(r) %*% c(cos(normV), sin(normV), 0))
  }
}

#' Inverse exponential map on the unit sphere.
#' 
#' Returns a vector v in the tangent space to the unit sphere at p, such that rayExp(p, v) = q.
#' @param p A ray.
#' @param q A ray.
#' @return A 3-dimensional vector, perpendicular to p, not necessarily unit.
rayLog <- function(p, q) {
  normV <- arcCos(dot(p, q))
  w <- rayNormalized(cross(p, q))
  v <- normV * cross(w, p)
  if (dot(v, q) >= 0)
    v
  else
    -v
}

# Tests that those are inverses.
#p <- rayUniform()
#v <- runif(1, min=0.1, max=3.1) * rayOrthogonalUniform(p)
#sqrt(dot(v, v))
#rayDistance(p, rayExp(p, v))
#v
#rayLog(p, rayExp(p, v))

#' Wrapping a plane around the unit sphere.
#' 
#' This function composes the exponential map with a non-canonical isomorphism to the plane R^2. The isomorphism is specified by the user through a rotation matrix R. The first row of R is regarded as a point p on the unit sphere. The other two rows define an isomorphism between R^2 and the tangent plane to the unit sphere at p. w is mapped through this isomorphism into the tangent plane, and then into the sphere via the exponential map. Inverse to rayTangentVectorFromPoint.
#' @param w A 2-dimensional real vector.
#' @param rotation A 3x3 real matrix (special orthogonal).
#' @return A ray.
rayPointFromTangentVector <- function(w, rotation) {
  v <- as.numeric(t(rotation) %*% c(0, w))
  q <- rayExp(rotation[1,], v)
  q
  # Does this other code do the same thing?
  #v <- c(sqrt(1 - vec[[1]]^2 - vec[[2]]^2), vec)
  #as.numeric(t(rotation) %*% v)
}

#' Unwrapping unit sphere into a tangent plane.
#' 
#' Inverse to rayPointFromTangentVector.
#' @param q A ray.
#' @param rotation A 3x3 real matrix (special orthogonal).
#' @return A 2-dimensional real vector.
rayTangentVectorFromPoint <- function(q, rotation) {
  v <- rayLog(rotation[1,], q)
  w <- as.numeric(rotation %*% v)[2:3]
  w
  # Does this faster code do the same thing?
  #as.numeric(rotation %*% q)[2:3]
}

# These tests demonstrate that the preceding two functions are inverses.
#p <- rayUniform()
#q <- rayUniform()
#perp <- rayOrthogonalUniform(p)
#rot <- rbind(p, perp, cross(p, perp))
#q
#w <- rayTangentVectorFromPoint(q, rot)
#rayPointFromTangentVector(w, rot)
#w
#rayTangentVectorFromPoint(rayPointFromTangentVector(w, rot), rot)

#' Elliptical percentile confidence region based on Mahalanobis distance in a tangent space.
#' 
#' @param us A list of rays.
#' @param uBar A ray. Usually something like rayProjectedMean(us).
#' @param alpha A real number, either 0.05 or 0.01.
#' @param numPoints A real number (non-negative integer). The resolution with which to sample the boundary curve of the (1 - alpha)* 100% confidence region.
#' @param doIsotropic Logical. If TRUE, forces the inverse covariance to the identity matrix and hence the region to be circular.
#' @return A list. $bootstraps is us. $center is uBar. $covarInv is the inverse covariance matrix in the tangent space, which is just the identity if doIsotropic is TRUE. rotation is the rotation used in rayTangentVectorFromPoint, etc. $q000, $q025, $q050, $q075, $q095, $q099 are quantiles of Mahalanobis norm. $pvalue is an R function that assigns to any ray its p-value. If numPoints > 0, then the list also has elements $alpha and $points. $alpha is alpha. $points is a list of numPoints + 1 rays describing the boundary of the confidence region, in order, with the first and last points identical.
rayMahalanobisInference <- function(us, uBar, alpha=0.05, numPoints=0, doIsotropic=FALSE) {
  # Map the rays into the tangent space at the mean.
  perp <- rayOrthogonalUniform(uBar)
  rot <- rbind(uBar, perp, cross(uBar, perp))
  vs <- lapply(us, rayTangentVectorFromPoint, rot)
  # Compute the usual covariance stuff.
  if (doIsotropic)
    covarInv <- diag(c(1, 1))
  else {
    covar <- arithmeticMean(lapply(vs, function(v) {outer(v, v)}))
    covarInv <- solve(covar)
  }
  norms <- sapply(vs, function(v) {sqrt(v %*% covarInv %*% v)})
  empiricalCDF <- ecdf(norms)
  # Build the p-value function.
  f <- function(u) {
    v <- rayTangentVectorFromPoint(u, rot)
    1 - empiricalCDF(sqrt(v %*% covarInv %*% v))
  }
  # Compute a few popular percentiles.
  qs <- quantile(norms, probs=c(0.00, 0.25, 0.50, 0.75, 0.95, 0.99, 1.00), names=FALSE)
  result <- list(bootstraps=us, pvalue=f, center=uBar, covarInv=covarInv, rotation=rot,
                 q000=qs[[1]], q025=qs[[2]], q050=qs[[3]], q075=qs[[4]], q095=qs[[5]], q099=qs[[6]], q100=qs[[7]])
  if (numPoints > 0) {
    # Prepare to generate points bounding the confidence region.
    if (alpha != 0.05 && alpha != 0.01)
      alpha <- 0.05
    if (alpha == 0.01)
      q <- qs[[6]]
    else
      q <- qs[[5]]
    # Make the ellipse in the 2D tangent space.
    eig <- eigen(covarInv, symmetric=TRUE)
    semiaxes <- eig$values^(-1 / 2)
    circle <- lapply(0:numPoints, function(s) c(cos(s * 2 * pi / numPoints), sin(s * 2 * pi / numPoints)))
    vs <- lapply(circle, function(v) as.numeric(eig$vectors %*% (q * semiaxes * v)))
    # Embed the tangent space in 3D and wrap it.
    result$points <- lapply(vs, rayPointFromTangentVector, rot)
    result$alpha <- alpha
  }
  result
}



### FISHER DISTRIBUTION ###

# should do more with Fisher distribution:
# unbiased estimation of isotropic Fisher parameters (Mardia, 1972, p. 250)
# isotropic Fisher goodness of fit (Mardia, 1972, p. 252)
# uniformity of rays 256
# Watson-Williams test p. 214
# Watson-Williams two-sample test 263
# assumes Fisher concentrations equal, so also test that (266)

# Generates one ray from the Fisher distribution centered at the third column of the rotation matrix.
rayFisherHelper <- function(rot, kappa) {
  # The azimuthal coordinate theta is uniform on the unit circle.
  theta <- runif(1, min=-pi, max=pi)
  # Meanwhile, phi is distributed as f = [(k / (2 sinh k)) e^(k cos phi) sin phi] on [0, pi].
  f <- function(phi) {
    (kappa / (2 * sinh(kappa))) * exp(kappa * cos(phi)) * sin(phi)
  }
  # Bound f.
  aa <- exp((sqrt(1 + 4 * kappa^2) - 1) / 2)
  bb <- sqrt(sqrt(1 + 4 * kappa^2) - 1)
  cc <- 2 * sqrt(2) * sinh(kappa)
  bound <- aa * bb / cc
  # Perform acceptance-rejection sampling.
  phi <- runif(1, min=0, max=pi)
  y <- runif(1, min=0, max=bound)
  while (y > f(phi)) {
    phi <- runif(1, min=0, max=pi)
    y <- runif(1, min=0, max=bound)
  }
  as.numeric(rot %*% cartesianFromSpherical(c(1, phi, theta)))
}

#' Sampling points from the Fisher distribution.
#' 
#' Also called von Mises-Fisher distribution or Langevin distribution. Uses expressions on p. 172 of Mardia and Jupp (2000), with a naive acceptance-rejection sampling algorithm. As kappa increases, this algorithm gets less and less efficient --- for example, about 19 tries per acceptance when kappa == 100.
#' @param mu A ray. The center of the distribution.
#' @param kappa A real number (positive). The concentration of the distribution.
#' @param n A real number (positive integer) or NULL.
#' @return If n is NULL, then a single ray. If n is a positive integer, then a list of n rays.
rayFisher <- function(mu, kappa, n=NULL) {
  # Make a rotation matrix with mu as its third column.
  nu <- rayOrthogonalUniform(mu)
  rot <- cbind(cross(nu, mu), nu, mu)
  if (is.null(n)) {
    rayFisherHelper(rot, kappa)
  } else
    replicate(n, rayFisherHelper(rot, kappa), simplify=FALSE)
}

# From Mardia and Jupp (2000), Appendix 3.2.
rayFisherMLEKappaHats <- c(
  0.000, 0.030, 0.060, 0.090, 0.120, 0.150, 0.180, 0.211, 0.241, 0.271,
  0.302, 0.332, 0.363, 0.394, 0.425, 0.456, 0.488, 0.519, 0.551, 0.583,
  0.615, 0.647, 0.680, 0.713, 0.746, 0.780, 0.814, 0.848, 0.883, 0.918,
  0.953, 0.989, 1.025, 1.062, 1.100, 1.137, 1.176, 1.215, 1.255, 1.295,
  1.336, 1.378, 1.421, 1.464, 1.508, 1.554, 1.600, 1.647, 1.696, 1.746,
  1.797, 1.849, 1.903, 1.958, 2.015, 2.074, 2.135, 2.198, 2.263, 2.330,
  2.401, 2.473, 2.549, 2.628, 2.711, 2.798, 2.888, 2.984, 3.085, 3.191,
  3.304, 3.423, 3.551, 3.687, 3.832, 3.989, 4.158, 4.341, 4.541, 4.759,
  4.998, 5.262, 5.555, 5.882, 6.250, 6.667, 7.143, 7.692, 8.333, 9.091,
  10.000, 11.111, 12.500, 14.286, 16.667, 20.000, 25.000, 33.333, 50.000, 100.000)
rayFisherMLEInterpolation <- approxfun(x=seq(from=0.00, to=0.99, by=0.01), y=rayFisherMLEKappaHats)

#' Maximum likelihood estimation of the Fisher parameters.
#' 
#' Based on Mardia and Jupp (2000, p. 198).
#' @param xs A list of rays.
#' @return A list with members $muHat (the mean ray, identical to rayProjectedMean), $rBar (non-negative real number), and $kappaHat (a positive real number).
rayFisherMLE <- function(xs) {
  xBar <- arithmeticMean(xs)
  rBar <- sqrt(dot(xBar, xBar))
  x0 <- xBar / rBar
  if (rBar > 0.9)
    kappaHat <- 1 / (1 - rBar)
  else
    kappaHat <- rayFisherMLEInterpolation(rBar)
  list(muHat=x0, rBar=rBar, kappaHat=kappaHat)
}

# Test.
#mu <- rayUniform()
#us <- rayFisher(mu=mu, kappa=10, n=100)
#rayFisherMLE(us)
#mu

# From Mardia and Jupp (2000, Appendix 3.1).
rayFisherConfidenceKappas <- c(
  seq(from=0.0, to=4.9, by=0.1),
  seq(from=5.0, to=9.8, by=0.2),
  seq(from=10.0, to=12.5, by=0.5),
  c(13.0, 14.0, 15.0, 20.0, 30.0, 40.0, 50.0, 100.0))
rayFisherConfidenceDeltas <- c(
  154.2, 152.9, 151.5, 150.0, 148.4, 146.6, 144.8, 142.8, 140.8, 138.6,
  136.3, 133.9, 131.4, 128.9, 126.2, 123.6, 120.9, 118.3, 115.6, 113.0,
  110.4, 107.9, 105.4, 103.1, 100.8, 98.6, 96.5, 94.5, 92.6, 90.8,
  89.0, 87.4, 85.8, 84.3, 82.8, 81.4, 80.1, 78.8, 77.6, 76.5,
  75.4, 74.3, 73.3, 72.3, 71.3, 70.4, 69.6, 68.7, 67.9, 67.1,
  66.4, 64.9, 63.6, 62.3, 61.1, 60.0, 58.9, 57.9, 56.9, 56.0,
  55.1, 54.3, 53.5, 52.7, 52.0, 51.3, 50.6, 50.0, 49.3, 48.7,
  48.2, 47.6, 47.1, 46.5, 46.0, 45.5, 44.4, 43.3, 42.3, 41.4,
  40.5, 39.7, 38.2, 36.8, 31.8, 25.8, 22.3, 19.9, 14.0) * degree
rayFisherConfidenceInterpolation <- approxfun(x=rayFisherConfidenceKappas, y=rayFisherConfidenceDeltas)

#' Confidence region for the Fisher distribution mean.
#' 
#' Theoretically better than rayFisherLargeSampleConfidence? But in my tests it seems no better. Coverage rates tend to be too small for small sample sizes n. n = 30 delivers around 93% coverage. n = 100 delivers adequate coverage. Based on Eq. (10.4.26) from Mardia and Jupp (2000).
#' @param xs A list of rays.
#' @return A list consisting of $angle95 (a real number in [0, pi]) and everything from rayFisherMLE.
rayFisherConfidence095 <- function(xs) {
  n <- length(xs)
  mle <- rayFisherMLE(xs)
  kappa <- n * mle$kappaHat * mle$rBar
  if (kappa > 100)
    mle$angle95 <- 140.2 * degree * kappa^(-0.5)
  else
    mle$angle95 <- rayFisherConfidenceInterpolation(kappa)
  mle
}

rayFisherConfidenceExperiment <- function(N, kappa, n, alpha=0.05) {
  f <- function(kappa, n, alpha) {
    mu <- rayUniform()
    us <- rayFisher(mu, kappa, n)
    uBar <- rayProjectedMean(us)
    delta <- rayFisherConfidence095(us)$angle95
    mu %*% uBar > cos(delta)
  }
  pHat <- sum(replicate(N, f(kappa, n, alpha))) / N
  se <- standardErrorProportion(N, pHat)
  c(pHat, pHat - 2 * se, pHat + 2 * se)
}
# These experiments show that kappa doesn't affect the accuracy much. But n = 30 is a bit small.
#rayFisherConfidenceExperiment(N=2000, kappa=1, n=10) # 0.8765000 0.8617862 0.8912138
#rayFisherConfidenceExperiment(N=2000, kappa=10, n=10) # 0.9005000 0.8871135 0.9138865
#rayFisherConfidenceExperiment(N=2000, kappa=100, n=10) # 0.9015000 0.8881735 0.9148265
#rayFisherConfidenceExperiment(N=2000, kappa=1, n=30) # 0.9270000 0.9153663 0.9386337
#rayFisherConfidenceExperiment(N=2000, kappa=10, n=30) # 0.9290000 0.9175144 0.9404856
#rayFisherConfidenceExperiment(N=2000, kappa=100, n=30) # 0.9240000 0.9121489 0.9358511
#rayFisherConfidenceExperiment(N=2000, kappa=1, n=100) # 0.9435000 0.9331745 0.9538255
#rayFisherConfidenceExperiment(N=2000, kappa=10, n=100) # 0.9405000 0.9299208 0.9510792
#rayFisherConfidenceExperiment(N=2000, kappa=100, n=100) # 0.9525000 0.9429875 0.9620125
#rayFisherConfidenceExperiment(N=2000, kappa=1, n=300) # 0.9490000 0.9391614 0.9588386
#rayFisherConfidenceExperiment(N=2000, kappa=10, n=300) # 0.9485000 0.9386159 0.9583841
#rayFisherConfidenceExperiment(N=2000, kappa=100, n=300) # 0.9525000 0.9429875 0.9620125

# This is for comparing our results to Allmendinger and Cardozo.
#mu <- rayUniform()
#u10s <- rayFisher(mu, 6, 10)
#geoTrendPlungeDegFileFromRays(u10s, "u10s.txt")
#u100s <- rayFisher(mu, 6, 100)
#geoTrendPlungeDegFileFromRays(u100s, "u100s.txt")
#rayFisherConfidence095(u10s)$angle95 / degree
#rayFisherLargeSampleConfidence(u10s)$angle / degree
#rayFisherConfidence095(u100s)$angle95 / degree
#rayFisherLargeSampleConfidence(u100s)$angle / degree

#' Asymptotic confidence region for the Fisher distribution mean.
#' 
#' Assumes large sample size n. When the sample size is small, this function's reported angle is typically too small. n = 10 is too small (only 90% coverage). n = 30 is a bit small (about 94% coverage). Even n = 55 is borderline. But n = 70 is quite adequate. Uses Eq. (10.4.31) from Mardia and Jupp (2000).
#' @param xs A list of rays.
#' @param alpha A real number (in [0, 1]). The significance level --- for example, 0.05 for 95% confidence.
#' @return A list with members $muHat (a ray, identical to rayProjectedMean), $rBar (a non-negative real number), $angle (a real number in [0, pi]). $angle is the radius of the confidence region, in radians, measured along the surface of the sphere.
rayFisherLargeSampleConfidence <- function(xs, alpha=0.05) {
  xBar <- arithmeticMean(xs)
  rBar <- sqrt(dot(xBar, xBar))
  x0 <- xBar / rBar
  tMatrix <- arithmeticMean(lapply(xs, function(x) outer(x, x)))
  n <- length(xs)
  sinDelta <- sqrt(-log(alpha) * (1 - x0 %*% tMatrix %*% x0) / (n * rBar^2))
  list(muHat=x0, rBar=rBar, angle=arcSin(sinDelta))
}

rayFisherLargeSampleConfidenceExperiment <- function(N, kappa, n, alpha=0.05) {
  f <- function(kappa, n, alpha) {
    mu <- rayUniform()
    us <- rayFisher(mu, kappa, n)
    uBar <- rayProjectedMean(us)
    delta <- rayFisherLargeSampleConfidence(us, alpha)$angle
    mu %*% uBar > cos(delta)
  }
  pHat <- sum(replicate(N, f(kappa, n, alpha))) / N
  se <- standardErrorProportion(N, pHat)
  c(pHat, pHat - 2 * se, pHat + 2 * se)
}
# These experiments show that kappa doesn't matter much to the accuracy. But n = 30 is a bit small.
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=1, n=10) # 0.8945000 0.8807618 0.9082382
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=1, n=30) # 0.9385000 0.9277559 0.9492441
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=1, n=100) # 0.9530000 0.9435352 0.9624648
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=1, n=300) # 0.9510000 0.9413461 0.9606539
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=10, n=10) # 0.902500 0.889234 0.915766
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=10, n=30) # 0.9375000 0.9266747 0.9483253
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=10, n=100) # 0.9550000 0.9457291 0.9642709
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=10, n=300) # 0.9535000 0.9440832 0.9629168
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=100, n=10) # 0.9045000 0.8913562 0.9176438
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=100, n=30) # 0.9370000 0.9261344 0.9478656
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=100, n=100) # 0.9455000 0.9353482 0.9556518
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=100, n=300) # 0.9435000 0.9331745 0.9538255
# Even n = 55 may be a bit small.
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=1, n=55) # 0.9500000 0.9402532 0.9597468
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=10, n=55) # 0.9400000 0.9293793 0.9506207
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=100, n=55) # 0.9335000 0.9223575 0.9446425
# But n = 70 is quite adequate.
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=1, n=70) # 0.9465000 0.9364364 0.9565636
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=10, n=70) # 0.9485000 0.9386159 0.9583841
#rayFisherLargeSampleConfidenceExperiment(N=2000, kappa=100, n=70) # 0.9475000 0.9375257 0.9574743

# See Tauxe (2010, p. 214).
rayFisherTauxe <- function(xs, alpha=0.05) {
  resultant <- Reduce("+", xs)
  r <- sqrt(dot(resultant, resultant))
  xBar <- resultant / r
  n <- length(xs)
  kappa <- (n - 1) / (n - r)
  angle <- arcCos(1 - (alpha^(1 / (1 - n)) - 1) * (n - r) / r)
  list(muHat=xBar, kappaHat=kappa, angle=angle)
}

rayFisherTauxeExperiment <- function(N, kappa, n, alpha=0.05) {
  f <- function(kappa, n, alpha) {
    mu <- rayUniform()
    us <- rayFisher(mu, kappa, n)
    tauxe <- rayFisherTauxe(us, alpha)
    mu %*% tauxe$muHat > cos(tauxe$angle)
  }
  pHat <- sum(replicate(N, f(kappa, n, alpha))) / N
  se <- standardErrorProportion(N, pHat)
  c(pHat, pHat - 2 * se, pHat + 2 * se)
}
# These experiments show that n doesn't affect the accuracy much. kappa == 1 is too dispersed, but kappa == 3 is fine.
#rayFisherTauxeExperiment(N=2000, kappa=1, n=5) # 0.9065000 0.8934802 0.9195198
#rayFisherTauxeExperiment(N=2000, kappa=3, n=5) # 0.9475000 0.9375257 0.9574743
#rayFisherTauxeExperiment(N=2000, kappa=10, n=5) # 0.9440000 0.9337176 0.9542824
#rayFisherTauxeExperiment(N=2000, kappa=30, n=5) # 0.9510000 0.9413461 0.9606539
#rayFisherTauxeExperiment(N=2000, kappa=100, n=5) # 0.9495000 0.9397072 0.9592928
#rayFisherTauxeExperiment(N=2000, kappa=1, n=10) # 0.9005000 0.8871135 0.9138865
#rayFisherTauxeExperiment(N=2000, kappa=3, n=10) # 0.9470000 0.9369809 0.9570191
#rayFisherTauxeExperiment(N=2000, kappa=10, n=10) # 0.9580000 0.9490294 0.9669706
#rayFisherTauxeExperiment(N=2000, kappa=30, n=10) # 0.9495000 0.9397072 0.9592928
#rayFisherTauxeExperiment(N=2000, kappa=100, n=10) # 0.9390000 0.9282968 0.9497032
#rayFisherTauxeExperiment(N=2000, kappa=1, n=30) # 0.8605000 0.8450055 0.8759945
#rayFisherTauxeExperiment(N=2000, kappa=3, n=30) # 0.9520000 0.9424401 0.9615599
#rayFisherTauxeExperiment(N=2000, kappa=10, n=30) # 0.9550000 0.9457291 0.9642709
#rayFisherTauxeExperiment(N=2000, kappa=30, n=30) # 0.9420000 0.9315467 0.9524533
#rayFisherTauxeExperiment(N=2000, kappa=100, n=30) # 0.9510000 0.9413461 0.9606539
#rayFisherTauxeExperiment(N=2000, kappa=1, n=100) # 0.8735000 0.8586341 0.8883659
#rayFisherTauxeExperiment(N=2000, kappa=3, n=100) # 0.9535000 0.9440832 0.9629168
#rayFisherTauxeExperiment(N=2000, kappa=10, n=100) # 0.9465000 0.9364364 0.9565636
#rayFisherTauxeExperiment(N=2000, kappa=30, n=100) # 0.9545000 0.9451802 0.9638198
#rayFisherTauxeExperiment(N=2000, kappa=100, n=100) # 0.9525000 0.9429875 0.9620125
#rayFisherTauxeExperiment(N=2000, kappa=1, n=300) # 0.8795000 0.8649412 0.8940588
#rayFisherTauxeExperiment(N=2000, kappa=3, n=300) # 0.9500000 0.9402532 0.9597468
#rayFisherTauxeExperiment(N=2000, kappa=10, n=300) # 0.9525000 0.9429875 0.9620125
#rayFisherTauxeExperiment(N=2000, kappa=30, n=300) # 0.9460000 0.9358922 0.9561078
#rayFisherTauxeExperiment(N=2000, kappa=100, n=300) # 0.9475000 0.9375257 0.9574743

# Tauxe's kappaHats are usually smaller than the other MLE's kappaHats (except when kappa == 1; then Tauxe's is often larger).
#us <- rayFisher(mu=rayUniform(), kappa=1, n=10)
#rayFisherMLE(us)$kappaHat
#rayFisherTauxe(us)$kappaHat



### KENT DISTRIBUTION ###

# There is some code in kent.R, but some of it works quite poorly.



### REGRESSION AND CURVE FITTING ###

# Altered from lineGeodesicRegression.
rayGeodesicRegression <- function(xs, ls, numSteps=1000, numPoints=0) {
  # Let l0 be the l whose x is closest to zero.
  l0 <- ls[[which.min(as.numeric(xs)^2)]]
  # Define the function to be minimized.
  n <- length(ls)
  e <- function(wb) {
    a <- rotExp(rotAntisymmetricFromVector(wb[1:3]))
    f <- function(i) {
      rayDistance(ls[[i]], as.numeric(a %*% c(cos(wb[[4]] * xs[[i]]), sin(wb[[4]] * xs[[i]]), 0)))^2
    }
    sum(sapply(1:n, f)) / (2 * n)
  }
  # Find the minimum, using the constant geodesic l0 as the seed.
  seed <- c(0, 0, 0, 0)
  solution <- optim(seed, e, hessian=TRUE, control=list(maxit=numSteps))
  # Report diagnostic information.
  eigvals <- eigen(solution$hessian, symmetric=TRUE, only.values=TRUE)$values
  rot <- rotExp(rotAntisymmetricFromVector(solution$par[1:3]))
  a <- solution$par[[4]]
  rSq <- 1 - solution$value / rayVariance(ls, rayProjectedMean(ls))
  result <- list(a=a, rotation=rot, error=solution$convergence, minEigenvalue=min(eigvals), rSquared=rSq)
  result$prediction <- function(x) {
    as.numeric(rot %*% c(cos(a * x), sin(a * x), 0))
  }
  if (numPoints >= 1)
    result$points <- lapply(seq(from=min(xs), to=max(xs), length.out=(numPoints + 1)), result$prediction)
  result
}

rayRescaledGeodesicRegression <- function(xs, ls, numSteps=1000, numPoints=0) {
  x0 <- min(xs)
  x1 <- max(xs)
  # Perform regression in scaled coordinates.
  regr <- rayGeodesicRegression(scales(xs), ls, numSteps, numPoints=0)
  # Scale the results back into the original coordinates.
  a <- regr$a / (x1 - x0)
  col1 <- regr$rotation[,1] * cos(a * x0) - regr$rotation[,2] * sin(a * x0)
  col2 <- regr$rotation[,1] * sin(a * x0) + regr$rotation[,2] * cos(a * x0)
  rot <- rotProjectedMatrix(cbind(col1, col2, cross(col1, col2)))
  result <- list(a=a, rotation=rot, error=regr$error, minEigenvalue=regr$minEigenvalue, rSquared=regr$rSquared)
  result$prediction <- function(x) {
    as.numeric(rot %*% c(cos(a * x), sin(a * x), 0))
  }
  if (numPoints >= 1)
    result$points <- lapply(seq(from=x0, to=x1, length.out=numPoints), result$prediction)
  result
}

rayGeodesicRegressionPermutations <- function(xs, ls, numPerms, numSteps=10000) {
  ys <- scales(xs)
  f <- function(i) {
    print(i / numPerms)
    regr <- rayGeodesicRegression(sample(ys, size=length(ys)), ls, numSteps)
    c(regr$error, regr$minEigenvalue, regr$rSquared)
  }
  perms <- sapply(1:numPerms, f)
  perms[3,][perms[1,] == 0 & perms[2,] > 0]
}

# Test.
#r <- rayUniform()
#rs <- lapply(rayGeodesicPoints(r, rayOrthogonalUniform(r), numSteps=10), rayFisher, 20)
#xs <- seq(from=0, to=1, length.out=length(rs))
#regr <- rayRescaledGeodesicRegression(xs, rs, numSteps=1000, numPoints=10)
#regr$error
#regr$minEigenvalue
#regr$rSquared
#regr$a / degree
#geoTrendPlungeDegFromCartesian(regr$rotation[,3])
#rayEqualAreaPlot(rs, colors=hues(xs), curves=list(regr$points))
#perms <- rayGeodesicRegressionPermutations(xs, rs, 1000)
#length(perms)
#sum(perms > regr$rSquared) / length(perms)

#' Fitting a small circle to some points on the unit sphere.
#' 
#' @param us A list of rays.
#' @return A list consisting of $pole (a ray, the pole to the small circle), $angle (the distance from the pole to the small circle, in radians), $error (0 if and only if minimization succeeds; if 1, then increase numSteps), and $minEigenvalue (worry about the result if this isn't positive).
rayBestFitSmallCircle <- function(us, numSeeds=5, numSteps=1000) {
  f <- function(phiTheta) {
    p <- cartesianFromSpherical(c(1, phiTheta))
    angles <- sapply(us, function(u) arcCos(dot(p, u)))
    var(angles)
  }
  # Find the minimum, starting from a few seeds.
  best <- list(value=(2 * pi^2))
  for (i in 1:numSeeds) {
    seed <- sphericalFromCartesian(rayUniform())[2:3]
    solution <- optim(seed, f, lower=c(0, -pi), upper=c(pi, pi), hessian=TRUE, control=list(maxit=numSteps), method="L-BFGS-B")
    if (solution$value <= best$value)
      best <- solution
  }
  # Report diagnostic information.
  eigvals <- eigen(best$hessian, symmetric=TRUE, only.values=TRUE)$values
  p <- cartesianFromSpherical(c(1, best$par))
  angles <- sapply(us, function(u) arcCos(dot(p, u)))
  a <- mean(angles)
  if (a > pi / 2) {
    a <- pi - a
    p <- -p
  }
  list(pole=p, angle=a, error=best$convergence, minEigenvalue=min(eigvals))
}

# Test with exact or inexact data. Works well.
#pole <- rayUniform()
#angle <- runif(1, 0, pi)
#us <- raySmallCircle(pole, angle, numPoints=10)
##us <- lapply(us, function(u) rayNormalized(u + rnorm(3, 0, 0.1)))
#rayFit <- rayBestFitSmallCircle(us, numSeeds=10)
#rayFit
#rayEqualAreaPlot(us, curves=list(raySmallCircle(rayFit$pole, rayFit$angle)))

raySmallCircleRegression <- function(xs, us, numSeeds=5, numSteps=1000, numPoints=0) {
  f <- function(phiThetaAlphaTauSigma) {
    pole <- cartesianFromSpherical(c(1, phiThetaAlphaTauSigma[1:2]))
    uOf0 <- cartesianFromSpherical(c(1, phiThetaAlphaTauSigma[4:5]))
    pred <- function(x) {
      as.numeric(rotMatrixFromAxisAngle(c(pole, x * phiThetaAlphaTauSigma[[3]])) %*% uOf0)
    }
    preds <- lapply(xs, pred)
    dists <- mapply(rayDistance, preds, us)
    dot(dists, dists) / (2 * length(us))
  }
  # Find the minimum, starting from a few seeds.
  best <- list(value=(2 * pi^2))
  for (i in 1:numSeeds) {
    seed <- c(sphericalFromCartesian(rayUniform())[2:3], runif(1, -pi, pi), sphericalFromCartesian(rayUniform())[2:3])
    solution <- optim(seed, f, lower=c(0, -pi, 0, 0, -pi), upper=c(pi, pi, pi, pi, pi),
                      hessian=TRUE, control=list(maxit=numSteps), method="L-BFGS-B")
    if (solution$value <= best$value)
      best <- solution
  }
  # Report results.
  eigvals <- eigen(best$hessian, symmetric=TRUE, only.values=TRUE)$values
  pole <- cartesianFromSpherical(c(1, best$par[1:2]))
  angle <- best$par[[3]]
  uOf0 <- cartesianFromSpherical(c(1, best$par[4:5]))
  var <- rayVariance(us, rayProjectedMean(us))
  rSquared <- 1 - best$value / var
  pred <- function(x) {as.numeric(rotMatrixFromAxisAngle(c(pole, x * angle)) %*% uOf0)}
  results <- list(error=best$convergence, minEigenvalue=min(eigvals),
                  pole=pole, angle=angle, rSquared=rSquared, prediction=pred)
  if (numPoints >= 1) {
    ys <- seq(from=min(xs), to=max(xs), length.out=numPoints)
    results$points <- lapply(ys, pred)
  }
  results
}

rayRescaledSmallCircleRegression <- function(xs, ls, ...) {
  # Perform regression in scaled coordinates.
  x0 <- min(xs)
  x1 <- max(xs)
  results <- raySmallCircleRegression(scales(xs), ls, ...)
  # Scale the results back into the original coordinates.
  results$rescaledPrediction <- results$prediction
  results$prediction <- function(x) {results$rescaledPrediction((x - x0) / (x1 - x0))}
  results$angle <- results$angle / (x1 - x0)
  results
}

raySmallCircleRegressionPermutations <- function(xs, ls, numPerms, ...) {
  ys <- scales(xs)
  f <- function(i) {
    print(i / numPerms)
    regr <- raySmallCircleRegression(sample(ys, size=length(ys)), ls, ...)
    c(regr$error, regr$minEigenvalue, regr$rSquared)
  }
  perms <- sapply(1:numPerms, f)
  perms[3,][perms[1,] == 0 & perms[2,] > 0]
}

# Demonstration and testing. Permutations don't work very well yet, for some reason.
raySmallCircleRegressionTest <- function(n=10, kappa=100, numPerms=100, numSteps=10000) {
  pole <- rayUniform()
  angle <- runif(1, 0, pi)
  uOf0 <- rayUniform()
  xs <- runif(n, -1, 1)
  xs <- sort(xs)
  pred <- function(x) as.numeric(rotMatrixFromAxisAngle(c(pole, x * angle)) %*% uOf0)
  preds <- lapply(xs, pred)
  #rayEqualAreaPlot(preds, colors=hues(xs))
  us <- lapply(preds, rayFisher, kappa)
  #rayEqualAreaPlot(us, colors=hues(xs))
  #print(xs)
  regr <- raySmallCircleRegression(xs, us, numPoints=20)
  print(c(regr$error, regr$minEigenvalue, regr$rSquared, regr$angle))
  print(regr$pole)
  rayEqualAreaPlot(us, colors=hues(xs), curves=list(regr$points))
  regr <- rayRescaledSmallCircleRegression(xs, us, numPoints=20)
  print(c(regr$error, regr$minEigenvalue, regr$rSquared, regr$angle))
  print(regr$pole)
  rayEqualAreaPlot(us, colors=hues(xs), curves=list(regr$points))
  rSquareds <- raySmallCircleRegressionPermutations(xs, us, numPerms=numPerms, numSteps=numSteps)
  print(c(length(rSquareds), sum(rSquareds > regr$rSquared)))
}
#raySmallCircleRegressionTest()



### PLOTTING ###

#' Equal-area plot of rays.
#'
#' @param points A list of rays, to be plotted as points.
#' @param curves A list of lists of rays. Each list of rays is regarded as a curve to be plotted. Each curve is automatically subidivided at the boundary. Lower-hemisphere parts appear solid, while upper-hemisphere appear dotted.
#' @param colors A character vector. Colors to be passed to the underlying R point-plotting function. They are truncated or recycled as needed.
#' @param shapes A character vector. Shapes to be assigned to the plotted points: "c", "s", or "t" for circle, square, or triangle. Lower-hemisphere (upper-) points are plotted using filled (unfilled) versions of the shapes. They are truncated or recycled as needed.
#' @return NULL.
rayEqualAreaPlot <- function(points=list(), curves=list(), colors=c("black"), shapes=c("c")) {
  # Convert the chosen shapes into their underlying shape codes, considering hemisphere.
  f <- function(i) {
    j <- (i - 1) %% length(shapes) + 1
    if (points[[i]][[3]] > 0) {
      if (shapes[[j]] == "s")
        0
      else if (shapes[[j]] == "t")
        2
      else
        1
    } else {
      if (shapes[[j]] == "s")
        15
      else if (shapes[[j]] == "t")
        17
      else
        19
    }
  }
  underShapes <- sapply(1:length(points), f)
  # All points will be positioned as if on the lower hemisphere.
  underPoints <- lapply(points, lower)
  # Subdivide curves into lower- and upper-hemispherical parts.
  underCurves <- list()
  styles <- c()
  for (curve in curves) {
    curvesSigns <- rayCurvesUpperLower(curve)
    for (i in 1:length(curvesSigns$curves))
      if (curvesSigns$signs[[i]] == 1) {
        underCurves[[length(underCurves) + 1]] <- lapply(curvesSigns$curves[[i]], function(v) -v)
        styles[[length(styles) + 1]] <- "dotted"
      } else {
        underCurves[[length(underCurves) + 1]] <- curvesSigns$curves[[i]]
        styles[[length(styles) + 1]] <- "solid"
      }
  }
  plotEqualArea(points=underPoints, curves=underCurves, colors=colors, shapes=underShapes, styles=styles)
}

#' Equal-area plot of rays, with a circle about each point.
#' 
#' @param points A list of rays, to be plotted as points.
#' @param radii A vector of real numbers. The radii of circles to be drawn about the points. Radii are measured in radians, along the surface of the unit sphere.
#' @param colors A character vector. See rayEqualAreaPlot.
#' @param shapes A character vector. See rayEqualAreaPlot.
#' @return NULL.
rayEqualAreaRadiusPlot <- function(points=list(), radii=c(), colors=c("black"), shapes=c("c")) {
  rayEqualAreaPlot(points, colors=colors, curves=thread(raySmallCircle, points, radii))
}

#' Convenience shortcut for plotting two sets of rays.
#' @param pointsA A list of rays.
#' @param pointsB A list of rays.
#' @param colorA Character. A color, of the sort used in all R graphics routines.
#' @param colorB Character. A color, of the sort used in all R graphics routines.
#' @param curves A list of lists of rays. Curves to be plotted.
#' @return NULL.
rayEqualAreaPlotTwo <- function(pointsA, pointsB, colorA="red", colorB="cyan", curves=list()) {
  rayEqualAreaPlot(points=c(pointsA, pointsB), curves=curves,
                   colors=c(replicate(length(pointsA), colorA), replicate(length(pointsB), colorB)))
}

#' Convenience shortcut for plotting three sets of rays.
#' @param pointsA A list of rays.
#' @param pointsB A list of rays.
#' @param pointsC A list of rays.
#' @param colorA Character. A color, of the sort used in all R graphics routines.
#' @param colorB Character. A color, of the sort used in all R graphics routines.
#' @param colorC Character. A color, of the sort used in all R graphics routines.
#' @param curves A list of lists of rays. Curves to be plotted.
#' @return NULL.
rayEqualAreaPlotThree <- function(pointsA, pointsB, pointsC, colorA="red", colorB="green", colorC="blue", curves=list()) {
  rayEqualAreaPlot(points=c(pointsA, pointsB, pointsC), curves=curves,
                   colors=c(replicate(length(pointsA), colorA),
                            replicate(length(pointsB), colorB),
                            replicate(length(pointsC), colorC)))
}



### TWO DIMENSIONS ###

# Some of the above machinery still works: rayNormalized, rayMeanScatter, rayProjectedMean, rayDistance.

#' Uniformly random points on the unit circle.
#' 
#' @param n A real number (positive integer) or NULL.
#' @return If n is NULL, then a single 2D ray. If n is a positive integer, then a list of n 2D rays.
ray2DUniform <- function(n=NULL) {
  if (is.null(n)) {
    angle <- runif(1, min=-pi, max=pi)
    c(cos(angle), sin(angle))
  }
  else
    replicate(n, rayUniform(), simplify=FALSE)
}

#' Random points on the unit circle, drawn from the wrapped normal distribution.
#' 
#' @param mean A 2D ray.
#' @param sd A real number (positive). The standard deviation sigma of the underlying normal distribution, in radians.
#' @param n A real number (positive integer) or NULL.
#' @return If n is NULL, then a single 2D ray. If n is a positive integer, then a list of n 2D rays.
ray2DWrappedNormal <- function(mean, sd, n=NULL) {
  if (is.null(n))
    ray2DWrappedNormal(mean, sd, 1)[[1]]
  else {
    angle <- atan2(mean[[2]], mean[[1]])
    angles <- rnorm(n, mean=angle, sd=sd)
    lapply(angles, function(a) c(cos(a), sin(a)))
  }
}

#' Bootstrapped extrinsic mean for 2D rays.
#' 
#' The inference is based on percentiles of distance from the mean of the bootstrapped means. The user should check that the bootstrapped means form a unimodal symmetric distribution, before taking such a region seriously.
#' @param ls A list of 2D rays.
#' @param numBoots A real number (positive integer). The number of bootstrapped means to compute. 10,000 might be a good value.
#' @return A list ($center, $bootstraps, $pvalue, $q000, $q025, $q050, $q075, $q095, $q099, $q100). $bootstraps are the bootstraps. $center is their mean. The other fields are quantiles of distance from the mean, in radians, among the bootstraps. For example, a 95% confidence region consists of all 2D rays within $q095 of $center. $pvalue is an R function from 2D rays to real numbers, assigning a p-value to any given null hypothesis for the mean.
ray2DBootstrapInference <- function(ls, numBoots) {
  boots <- replicate(numBoots, rayProjectedMean(sample(ls, length(ls), replace=TRUE)), simplify=FALSE)
  bootMean <- rayProjectedMean(boots)
  dists <- sapply(boots, rayDistance, bootMean)
  empiricalCDF <- ecdf(dists)
  # Build the p-value function.
  f <- function(u) {
    1 - empiricalCDF(rayDistance(u, bootMean))
  }
  # Compute a few popular percentiles.
  qs <- quantile(dists, probs=c(0.00, 0.25, 0.50, 0.75, 0.95, 0.99, 1.00), names=FALSE)
  list(center=bootMean, bootstraps=boots, pvalue=f,
       q000=qs[[1]], q025=qs[[2]], q050=qs[[3]], q075=qs[[4]], q095=qs[[5]], q099=qs[[6]], q100=qs[[7]])
}

#' Rose plot for 2D rays.
#' 
#' The data are binned. (Binning starts at the positive x-axis and proceeds counterclockwise. In a future release maybe we could add an offset, to start the binning elsewhere.) Each datum can be given a weight, which is tantamount to repeating the datum in the data set. Each petal in the Rose plot represents the total weight in a bin, either by area or by length. Concentric circles indictate the 10%, 20%, 30%, etc. weight levels.
#' @param rays A list of 2D or 3D real vectors. In each one, only the 2D projection (the first two components) is used. It must be non-zero but need not be unit.
#' @param weights A vector of real numbers. The weights to attach to the rays.
#' @param numBins The number of bins to use. For example, numBins == 36 means 10-degree-wide bins.
#' @param inner A real number (non-negative). This parameter reserves a circle of empty space at the center of the plot. More precisely, the plot takes place between radius inner and radius 1. inner == 0 seems traditional in geology, but I feel that inner == 0.25 (say) is prettier and easier to read.
#' @param areal Logical. If FALSE, then each petal's length is proportional to its bin's weight. If TRUE, then each petal's area is proportional to its bin's weight. areal == FALSE seems traditional in geology, but all data visualization advice I've ever seen suggests that areal == TRUE is the right choice.
#' @return NULL.
ray2DRosePlot <- function(rays, weights=replicate(length(rays), 1), numBins=36, inner=0, areal=TRUE) {
  # Bin the data, filling each bin with its fraction of the total weight.
  bins <- replicate(numBins, 0)
  for (i in 1:length(rays)) {
    heading <- atan2(rays[[i]][[2]], rays[[i]][[1]]) %% (2 * pi)
    bin <- heading %/% (2 * pi / numBins) + 1
    bins[[bin]] <- bins[[bin]] + weights[[i]]
  }
  bins <- bins / sum(weights)
  # Prepare to scale bins and benchmarks.
  if (areal)
    radius <- function(w) {sqrt(w * (1 - inner^2) + inner^2)}
  else
    radius <- function(w) {w * (1 - inner) + inner}
  # Draw benchmarks.
  plot.new()
  plot.window(xlim=c(-1, 1), ylim=c(-1, 1))
  xs <- cos((0:360) * (2 * pi / 360))
  ys <- sin((0:360) * (2 * pi / 360))
  for (w in seq(from=0, to=1, by=0.1))
    lines(radius(w) * xs, radius(w) * ys)
  # Draw the petals.
  for (bin in 1:numBins) {
    headings <- seq(from=((bin - 1) * 2 * pi / numBins), to=(bin * 2 * pi / numBins), by=(2 * pi / 360))
    xs <- radius(bins[[bin]]) * cos(headings)
    ys <- radius(bins[[bin]]) * sin(headings)
    heading <- bin * 2 * pi / numBins
    r <- radius(max(bins[[bin]], bins[[(bin %% numBins) + 1]]))
    xs <- c(xs, inner * cos(heading), r * cos(heading))
    ys <- c(ys, inner * sin(heading), r * sin(heading))
    lines(xs, ys)
  }
}


