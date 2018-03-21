


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# A line is expressed as a unit 3D vector in Cartesian coordinates, as an R vector u = c(x, y, z) where x^2 + y^2 + z^2 == 1. It is crucial to remember that u and -u represent the same line. Many of the line functions here are implemented in terms of the ray functions of rays.R. In particular, when lines are tightly concentrated, you can often ignore their 'negative copies' and treat them as rays, with no appreciable effect on the statistics.



### MISCELLANY ###

#' Distance between two lines, as the angle (in radians) between them.
#' 
#' @param u A line.
#' @param v A line.
#' @return A real number. The distance between the two lines, in radians. Between 0 and pi / 2, inclusive.
lineDistance <- function(u, v) {
  arcCos(abs(dot(u, v)))
}

#' L^2 variance of a set of lines about a point.
#' 
#' I'm not sure about the weighting on this. If you change it, change regression too.
#' @param us A list of lines.
#' @param center A line. Usually some kind of mean of the us.
#' @return A real number (in the interval [0, pi^2 / 8]).
lineVariance <- function(us, center) {
  sum(sapply(us, function(u) lineDistance(u, center)^2)) / (2 * length(us))
}

#' Geodesic curve between two points on the unit sphere.
#' 
#' @param u A line.
#' @param v A line.
#' @param numSteps The number of line segments to be used in the approximation of the great circle arc from u to v.
#' @return A list of numSteps+1 lines. The first is u and the last is v. They are evenly spaced and in order.
lineGeodesicPoints <- function(u, v, numSteps=10) {
  if (dot(u, v) < 0)
    rayGeodesicPoints(u, -v, numSteps)
  else 
    rayGeodesicPoints(u, v, numSteps)
}

#' Uniformly random lines.
#' 
#' @param n A real number (positive integer) or NULL.
#' @return If n is NULL, then a single line. If n is a positive integer, then a list of n lines.
lineUniform <- function(n=NULL) {
  if (is.null(n))
    lower(rayUniform())
  else
    lapply(rayUniform(n), lower)
}

# Bingham test of uniformity, Mardia and Jupp (2000, p. 232)!!!
# Gine test of uniformity, p. 233



### EXTRINSIC METHODS ###

#' Projected mean and scatter of a set of lines.
#' 
#' @param us A list of lines.
#' @return A list with elements $values, $vectors. $vectors is a 3x3 real matrix whose columns are the principal directions of dispersion: mean, other direction along main girdle, pole to girdle. $values is a vector of three numbers, nonnegative, summing to 1, and decreasing. They capture the strength of the concentration about the principal directions.
lineMeanScatter <- function(us) {
  tMatrices <- lapply(us, function(u) outer(u, u))
  tMatrix <- arithmeticMean(tMatrices)
  eig <- eigen(tMatrix, symmetric=TRUE)
  eig
}

#' Shortcut convenience function for the projected mean.
#' 
#' @param us A list of lines.
#' @return A line.
lineProjectedMean <- function(us) {
  eig <- lineMeanScatter(us)
  eig$vectors[,1]
}

#' Bootstrapped projected mean with percentile confidence region and hypothesis tests.
#' 
#' The inference is based on percentiles of Mahalanobis distance in the tangent space at the mean of the bootstrapped means. The user should check that the bootstrapped means form a tight ellipsoidal cluster, before taking such a region seriously.
#' @param ls A list of lines.
#' @param numBoots A real number (positive integer). The number of bootstrapped means to compute. 10,000 might be a good value.
#' @param ... Other arguments to be passed to the underlying rayMahalanobisInference function.
#' @return A list. See rayMahalanobisInference for most of it. But the $pvalue function treats its input as a ray. An added $pvalueLine function treats its input as a line, so use that.
lineBootstrapInference <- function(ls, numBoots, ...) {
  boots <- replicate(numBoots, lineProjectedMean(sample(ls, length(ls), replace=TRUE)), simplify=FALSE)
  bootMean <- lineProjectedMean(boots)
  boots <- lapply(boots, function(u) {if (dot(u, bootMean) >= 0) u else -u})
  inf <- rayMahalanobisInference(boots, bootMean, ...)
  inf$pvalueLine <- function(l) {
    if (dot(l, inf$center) < 0)
      inf$pvalue(-l)
    else
      inf$pvalue(l)
  }
  inf
}



### WATSON DISTRIBUTION ###

# Helper function for lineWatson. Returns one line sampled from the Watson distribution.
lineWatsonHelper <- function(f, bound, rot) {
  phi <- runif(1, min=0, max=pi)
  y <- runif(1, min=0, max=bound)
  while (y > f(phi)) {
    phi <- runif(1, min=0, max=pi)
    y <- runif(1, min=0, max=bound)
  }
  theta <- runif(1, min=-pi, max=pi)
  as.numeric(rot %*% cartesianFromSpherical(c(1, phi, theta)))
}

#' Normalizing constant for the Watson distribution.
#' 
#' That is, the number C such that C exp(kappa (mu^T u)^2) integrates to 1. By the way, C == 1 / 1F1(0.5, 1.5, kappa).
#' @param kappa A real number. The concentration parameter.
#' @return A real number.
lineWatsonNormalizer <- function(kappa) {
  if (kappa >= 0)
    (2 * sqrt(kappa)) / (sqrt(pi) * as.numeric(erfi(sqrt(kappa))))
  else
    (2 * sqrt(-kappa)) / (sqrt(pi) * as.numeric(erf(sqrt(-kappa))))
}

#' Sampling lines from the Watson distribution.
#' 
#' A naive acceptance-rejection sampling algorithm, based on bounding the density (with respect to the distance from mu) with a constant. For large kappa, this method grows inefficient. For kappa == 100, about 13 tries are needed per success. For kappa == -100, about 18 tries are needed.
#' @param mu A line. The mean of the distribution (if kappa > 0) or the pole to the girdle of the distribution (if kappa < 0).
#' @param kappa A real number. The concentration parameter.
#' @param n A real number (positive integer) or NULL.
#' @return If n is NULL, then a single line. If n is a positive integer, then a list of n lines.
lineWatson <- function(mu, kappa, n=NULL) {
  nu <- rayOrthogonalUniform(mu)
  rot <- cbind(cross(nu, mu), nu, mu)
  cc <- lineWatsonNormalizer(kappa)
  f <- function(phi) {exp(kappa * cos(phi)^2) * sin(phi) * cc / 2}
  if (kappa <= 0.5)
    bound <- cc / 2
  else
    bound <- exp(kappa - 1 / 2) * (2 * kappa)^(-1 / 2) * cc / 2
  if (is.null(n))
    lineWatsonHelper(f, bound, rot)
  else
    replicate(n, lineWatsonHelper(f, bound, rot), simplify=FALSE)
}

lineWatsonMLED3s <- c(0.001, seq(from=0.005, to=0.330, by=0.005), 0.333,
                      seq(from=0.34, to=0.99, by=0.01), 0.995, 0.999)
lineWatsonMLEKappaHats <- c(
  -500, -100, -50, -33.33, -25, -20, -16.67, -14.29, -12.25, -11.11, -9.992, -9.087, -8.327, -7.681, -7.126,
  -6.641, -6.215, -5.836, -5.495, -5.188, -4.908, -4.651, -4.415, -4.196, -3.993, -3.802, -3.624, -3.457,
  -3.298, -3.148, -3.006, -2.870, -2.741, -2.617, -2.499, -2.385, -2.275, -2.170, -2.068, -1.970, -1.874, -1.782, -1.692,
  -1.605, -1.520, -1.438, -1.357, -1.279, -1.202, -1.127, -1.053, -0.982, -0.911, -0.842, -0.774, -0.708,
  -0.642, -0.578, -0.514, -0.452, -0.390, -0.330, -0.270, -0.211, -0.152, -0.095, -0.038, -0.004,
  0.075, 0.184, 0.292, 0.398, 0.503, 0.606, 0.708, 0.809, 0.909, 1.008, 1.106, 1.204, 1.302, 1.399, 1.497,
  1.594, 1.692, 1.790, 1.888, 1.987, 2.087, 2.188, 2.289, 2.392, 2.496, 2.602, 2.709, 2.819,
  2.930, 3.044, 3.160, 3.280, 3.402, 3.529, 3.659, 3.764, 3.934, 4.079, 4.231, 4.389, 4.556, 4.731, 4.917,
  5.115, 5.326, 5.552, 5.797, 6.063, 6.354, 6.676, 7.035, 7.438, 7.897, 8.426, 9.043, 9.776,
  10.654, 11.746, 13.112, 14.878, 17.242, 20.560, 25.546, 33.866, 50.521, 100.510, 200.5, 1000.5)
lineWatsonMLEInterpolation <- approxfun(x=lineWatsonMLED3s, y=lineWatsonMLEKappaHats)

#' Maximum likelihood estimation of the Watson distribution parameters.
#' 
#' From Mardia and Jupp (2000, Section 10.3.2).
#' @param xs A list of lines.
#' @param shape NULL or character, either 'bipolar' or 'girdle'. If NULL, then this function chooses automatically.
#' @return A list with members $muHat (a line, the MLE of the mean), $kappaHat (a real number, the MLE of the concentration, shape (character, either 'bipolar' or 'girdle'), d3 (a positive real number, the D3 from which kappaHat was computed), and $eigenvalues (the eigenvalues of the T-bar matrix, in descending order).
lineWatsonMLE <- function(xs, shape=NULL) {
  # In eigen, the eigenvectors are descending and the eigenvectors are unit.
  tBar <- arithmeticMean(lapply(xs, function(x) outer(x, x)))
  eig <- eigen(tBar, symmetric=TRUE)
  # Pick a shape if necessary.
  if (is.null(shape)) {
    if (eig$values[[1]] - eig$values[[2]] >= eig$values[[2]] - eig$values[[3]])
      shape <- "bipolar"
    else
      shape <- "girdle"
  }
  # Shape determines which eigenvector is picked for muHat.
  if (shape == "bipolar")
    muHat <- eig$vectors[,1]
  else
    muHat <- eig$vectors[,3]
  # Get kappaHat from the lookup table.
  d3 <- as.numeric(muHat %*% tBar %*% muHat)
  if (d3 > 0.9)
    kappaHat <- 1 / (1 - d3)
  else if (d3 < 0.05)
    kappaHat <- -1 / (2 * d3)
  else
    kappaHat <- lineWatsonMLEInterpolation(d3)
  list(muHat=muHat, kappaHat=kappaHat, shape=shape, d3=d3, eigenvalues=eig$values)
}

# Test that the MLE and the sampling are both working. Try both positive and negative kappa.
#mu <- lineUniform()
#kappa <- -100
#n <- 1000
#mu
#lineWatsonMLE(lineWatson(mu, kappa, n))

#' One-sample inference about the mean of the Watson distribution.
#' 
#' Assumes large concentration --- either kappa >> 0 or kappa << 0. From Mardia and Jupp (2000, Section 10.7.3).
#' @param xs A list of lines.
#' @param alpha A real number, between 0 and 1. The significance level for the confidence region.
#' @param shape NULL or character, either 'bipolar' or 'girdle'. If NULL, then this function chooses automatically.
#' @return A list with members $shape, $tBar, $rhs, $pvalue. $shape is either 'bipolar' or 'girdle'. If 'bipolar', then the confidence region consists of all lines u such that u^T %*% $tBar %*% u > $rhs. If 'girdle', then the confidence region consists of all lines u such that u^T %*% $tBar %*% u < $rhs. $pvalue is an R function that takes as input a line u0 and produces as output a real number in [0, 1] --- the p-value for the null hypothesis that the Watson mean is u0.
lineWatsonInference <- function(xs, alpha=0.05, shape=NULL) {
  # In eigen, the eigenvectors are descending and the eigenvectors are unit.
  tBar <- arithmeticMean(lapply(xs, function(x) outer(x, x)))
  eig <- eigen(tBar, symmetric=TRUE)
  # Pick a shape if necessary.
  if (is.null(shape)) {
    if (eig$values[[1]] - eig$values[[2]] >= eig$values[[2]] - eig$values[[3]])
      shape <- "bipolar"
    else
      shape <- "girdle"
  }
  n <-  length(xs)
  if (shape == "bipolar") {
    t1 <- eig$values[[1]]
    rhs <- t1 + (t1 - 1) * qf(alpha, 2, 2 * n - 2, lower.tail=FALSE) / (n - 1)
    func <- function(u0) {
      f <- as.numeric((t1 - u0 %*% tBar %*% u0) * (n - 1) / (1 - t1))
      1 - pf(f, 2, 2 * n - 2)
    }
  } else {
    t3 <- eig$values[[3]]
    rhs <- t3 * (1 + qf(alpha, 2, n - 2, lower.tail=FALSE) * 2 / (n - 2))
    func <- function(u0) {
      f <- as.numeric((u0 %*% tBar %*% u0 - t3) * (n - 2) / (2 * t3))
      1 - pf(f, 2, n - 2)
    }
  }
  list(shape=shape, tBar=tBar, rhs=rhs, pvalue=func)
}

# Tests. Works well for kappa >= 10 and n >= 10. Works well for kappa <= -10 and n >= 100. Works poorly in other cases, including all |kappa| == 1 cases.
lineWatsonInferenceExperiment <-  function(kappa, n, N) {
  f <- function(kappa, n) {
    mu <- lineUniform()
    us <- lineWatson(mu, kappa, n)
    inf <- lineWatsonInference(us)
    lhs <- as.numeric(mu %*% inf$tBar %*% mu)
    if (inf$shape == "bipolar")
      c(lhs - inf$rhs, inf$pvalue(mu))
    else
      c(inf$rhs - lhs, inf$pvalue(mu))
  }
  results <- replicate(N, f(kappa, n))
  coverageRateRegions <- sum(results[1,] > 0) / N
  coverageRatePvalues <- sum(results[2,] > 0.05) / N
  mismatches1 <- sum(results[1,] > 0 & results[2,] < 0.05)
  mismatches2 <- sum(results[1,] < 0 & results[2,] > 0.05)
  c(coverageRateRegions, coverageRatePvalues, mismatches1, mismatches2)
}
#lineWatsonInferenceExperiment(kappa=-100, n=1000, N=2000) # 0.943 0.943 0.000 0.000
#lineWatsonInferenceExperiment(kappa=-10, n=1000, N=2000) # 0.947 0.947 0.000 0.000
#lineWatsonInferenceExperiment(kappa=-1, n=1000, N=2000) # 0.7735 0.7735 0.0000 0.0000
#lineWatsonInferenceExperiment(kappa=1, n=1000, N=2000) # 0.8295 0.8295 0.0000 0.0000
#lineWatsonInferenceExperiment(kappa=10, n=1000, N=2000) # 0.969 0.969 0.000 0.000
#lineWatsonInferenceExperiment(kappa=100, n=1000, N=2000) # 0.9485 0.9485 0.0000 0.0000
#lineWatsonInferenceExperiment(kappa=-100, n=100, N=2000) # 0.9515 0.9515 0.0000 0.0000
#lineWatsonInferenceExperiment(kappa=-10, n=100, N=2000) # 0.9485 0.9485 0.0000 0.0000
#lineWatsonInferenceExperiment(kappa=-1, n=100, N=2000) # 0.5155 0.5155 0.0000 0.0000
#lineWatsonInferenceExperiment(kappa=1, n=100, N=2000) # 0.7095 0.7095 0.0000 0.0000
#lineWatsonInferenceExperiment(kappa=10, n=100, N=2000) # 0.96 0.96 0.00 0.00
#lineWatsonInferenceExperiment(kappa=100, n=100, N=2000) # 0.948 0.948 0.000 0.000
#lineWatsonInferenceExperiment(kappa=-100, n=30, N=2000) # 0.9085 0.9085 0.0000 0.0000
#lineWatsonInferenceExperiment(kappa=-10, n=30, N=2000) # 0.8795 0.8795 0.0000 0.0000
#lineWatsonInferenceExperiment(kappa=-1, n=30, N=2000) # 0.3645 0.3645 0.0000 0.0000
#lineWatsonInferenceExperiment(kappa=1, n=30, N=2000) # 0.5425 0.5425 0.0000 0.0000
#lineWatsonInferenceExperiment(kappa=10, n=30, N=2000) # 0.9565 0.9565 0.0000 0.0000
#lineWatsonInferenceExperiment(kappa=100, n=30, N=2000) # 0.962 0.962 0.000 0.000
#lineWatsonInferenceExperiment(kappa=-100, n=10, N=2000) # 0.636 0.636 0.000 0.000
#lineWatsonInferenceExperiment(kappa=-10, n=10, N=2000) # 0.5625 0.5625 0.0000 0.0000
#lineWatsonInferenceExperiment(kappa=-1, n=10, N=2000) # 0.4575 0.4575 0.0000 0.0000
#lineWatsonInferenceExperiment(kappa=1, n=10, N=2000) # 0.6255 0.6255 0.0000 0.0000
#lineWatsonInferenceExperiment(kappa=10, n=10, N=2000) # 0.9655 0.9655 0.0000 0.0000
#lineWatsonInferenceExperiment(kappa=100, n=10, N=2000) # 0.9465 0.9465 0.0000 0.0000

#' Multi-sample inference about the mean of the Watson distribution.
#' 
#' Assumes large sample sizes. Assumes that all of the data sets are Watson-distributed with the same unknown concentration kappa. From Mardia and Jupp (2000, Section 10.7.4).
#' @param xss A list of lists of lines.
#' @param shape NULL or character, either 'bipolar' or 'girdle'. If NULL, then this function chooses automatically.
#' @return A real number. The p-value for the null hypothesis that the means of the q distributions are all equal.
lineLargeMultiSampleWatsonInference <- function(xss, shape=NULL) {
  q <- length(xss)
  ns <- sapply(xss, length)
  tBars <- lapply(xss, function(xs) arithmeticMean(lapply(xs, function(x) outer(x, x))))
  n <- sum(ns)
  tBar <- arithmeticMean(lapply(unlist(xss, recursive=FALSE), function(x) outer(x, x)))
  eig <- eigen(tBar, symmetric=TRUE)
  # Pick a shape if necessary.
  if (is.null(shape)) {
    if (eig$values[[1]] - eig$values[[2]] >= eig$values[[2]] - eig$values[[3]])
      shape <- "bipolar"
    else
      shape <- "girdle"
  }
  if (shape == "bipolar") {
    t1 <- eig$vectors[,1]
    eHatT4 <- sum(sapply(xss, function(xs) sum(sapply(xs, function(x) dot(x, t1)^4)))) / n
    t1s <- sapply(tBars, function(mat) eigen(mat, symmetric=TRUE, only.values=TRUE)$values[[1]])
    t1 <- eig$values[[1]]
    chiSq <- (3 * t1 - 1) * (dot(ns, t1s) - n * t1) / (2 * (t1 - eHatT4))
    1 - pchisq(chiSq, 2 * q - 2)
  } else {
    t3 <- eig$vectors[,3]
    eHatT4 <- sum(sapply(xss, function(xs) sum(sapply(xs, function(x) dot(x, t3)^4)))) / n
    t3s <- sapply(tBars, function(mat) eigen(mat, symmetric=TRUE, only.values=TRUE)$values[[3]])
    t3 <- eig$values[[3]]
    chiSq <- (1 - 3 * t3) * (n * t3 - dot(ns, t3s)) / (2 * (t3 - eHatT4))
    1 - pchisq(chiSq, 2 * q - 2)
  }
}

#' Multi-sample inference about the mean of the Watson distribution.
#' 
#' Assumes tight concentration. Tests suggest that it works for |kappa| >= 10, but not for |kappa| == 1. Assumes that all of the data sets are Watson-distributed with the same unknown concentration kappa. From Mardia and Jupp (2000, Section 10.7.4).
#' @param xss A list of lists of lines.
#' @param shape NULL or character, either 'bipolar' or 'girdle'. If NULL, then this function chooses automatically.
#' @return A real number. The p-value for the null hypothesis that the means of the q distributions are all equal.
lineConcentratedMultiSampleWatsonInference <- function(xss, shape=NULL) {
  q <- length(xss)
  ns <- sapply(xss, length)
  tBars <- lapply(xss, function(xs) arithmeticMean(lapply(xs, function(x) outer(x, x))))
  n <- sum(ns)
  tBar <- arithmeticMean(lapply(unlist(xss, recursive=FALSE), function(x) outer(x, x)))
  eig <- eigen(tBar, symmetric=TRUE)
  # Pick a shape if necessary.
  if (is.null(shape)) {
    if (eig$values[[1]] - eig$values[[2]] >= eig$values[[2]] - eig$values[[3]])
      shape <- "bipolar"
    else
      shape <- "girdle"
  }
  if (shape == "bipolar") {
    t1s <- sapply(tBars, function(mat) eigen(mat, symmetric=TRUE, only.values=TRUE)$values[[1]])
    t1 <- eig$values[[1]]
    f <- (dot(ns, t1s) - n * t1) * (n - q) / (dot(ns, 1 - t1s) * (q - 1))
    1 - pf(f, 2 * q - 2, 2 * n - 2 * q)
  } else {
    t3s <- sapply(tBars, function(mat) eigen(mat, symmetric=TRUE, only.values=TRUE)$values[[3]])
    t3 <- eig$values[[3]]
    f <- (n * t3 - dot(ns, t3s)) * (n - 2 * q) / (dot(ns, t3s) * (2 * q - 2))
    1 - pf(f, 2 * q - 2, n - 2 * q)
  }
}

# Tests. The large-sample-size method is always too conservative, even at n = 3000? The high-concentration method works for |kappa| >= 10, for all sample sizes, but not for |kappa| == 1.
lineMultiSampleWatsonInferenceExperiment <-  function(kappa, n, N, q) {
  f <- function(kappa, n, q) {
    mu <- lineUniform()
    uss <- replicate(q, lineWatson(mu, kappa, n), simplify=FALSE)
    large <- lineLargeMultiSampleWatsonInference(uss)
    concen <- lineConcentratedMultiSampleWatsonInference(uss)
    c(large, concen)
  }
  results <- replicate(N, f(kappa, n, q))
  rateLarge <- sum(results[1,] > 0.05) / N
  rateConcen <- sum(results[2,] > 0.05) / N
  c(rateLarge, rateConcen)
}
#lineMultiSampleWatsonInferenceExperiment(kappa=-100, n=10, N=2000, q=2) # 0.9935 0.9490
#lineMultiSampleWatsonInferenceExperiment(kappa=-10, n=10, N=2000, q=2) # 0.9975 0.9460
#lineMultiSampleWatsonInferenceExperiment(kappa=-3, n=10, N=2000, q=2) # 0.9990 0.9335
#lineMultiSampleWatsonInferenceExperiment(kappa=-1, n=10, N=2000, q=2) # 0.999 0.866
#lineMultiSampleWatsonInferenceExperiment(kappa=1, n=10, N=2000, q=2) # 0.9995 0.8855
#lineMultiSampleWatsonInferenceExperiment(kappa=3, n=10, N=2000, q=2) # 0.9985 0.9615
#lineMultiSampleWatsonInferenceExperiment(kappa=10, n=10, N=2000, q=2) # 0.9995 0.9615
#lineMultiSampleWatsonInferenceExperiment(kappa=100, n=10, N=2000, q=2) # 0.997 0.953
#lineMultiSampleWatsonInferenceExperiment(kappa=-100, n=30, N=2000, q=2) # 0.9975 0.9515
#lineMultiSampleWatsonInferenceExperiment(kappa=-10, n=30, N=2000, q=2) # 0.9980 0.9505
#lineMultiSampleWatsonInferenceExperiment(kappa=-3, n=30, N=2000, q=2) # 0.9970 0.9285
#lineMultiSampleWatsonInferenceExperiment(kappa=-1, n=30, N=2000, q=2) # 0.9985 0.7645
#lineMultiSampleWatsonInferenceExperiment(kappa=1, n=30, N=2000, q=2) # 0.9985 0.8115
#lineMultiSampleWatsonInferenceExperiment(kappa=3, n=30, N=2000, q=2) # 0.9990 0.9605
#lineMultiSampleWatsonInferenceExperiment(kappa=10, n=30, N=2000, q=2) # 0.9965 0.9545
#lineMultiSampleWatsonInferenceExperiment(kappa=100, n=30, N=2000, q=2) # 0.9975 0.9460
#lineMultiSampleWatsonInferenceExperiment(kappa=-100, n=100, N=2000, q=2) # 0.999 0.961
#lineMultiSampleWatsonInferenceExperiment(kappa=-10, n=100, N=2000, q=2) # 0.9975 0.9480
#lineMultiSampleWatsonInferenceExperiment(kappa=-3, n=100, N=2000, q=2) # 0.9965 0.9220
#lineMultiSampleWatsonInferenceExperiment(kappa=-1, n=100, N=2000, q=2) # 0.9990 0.7505
#lineMultiSampleWatsonInferenceExperiment(kappa=1, n=100, N=2000, q=2) # 0.996 0.795
#lineMultiSampleWatsonInferenceExperiment(kappa=3, n=100, N=2000, q=2) # 0.9975 0.9565
#lineMultiSampleWatsonInferenceExperiment(kappa=10, n=100, N=2000, q=2) # 0.9985 0.9625
#lineMultiSampleWatsonInferenceExperiment(kappa=100, n=100, N=2000, q=2) # 0.9970 0.9495
#lineMultiSampleWatsonInferenceExperiment(kappa=-100, n=300, N=2000, q=2) # 0.9980 0.9485
#lineMultiSampleWatsonInferenceExperiment(kappa=-10, n=300, N=2000, q=2) # 0.9980 0.9555
#lineMultiSampleWatsonInferenceExperiment(kappa=-3, n=300, N=2000, q=2) # 0.9985 0.9410
#lineMultiSampleWatsonInferenceExperiment(kappa=-1, n=300, N=2000, q=2) # 0.9980 0.7735
#lineMultiSampleWatsonInferenceExperiment(kappa=1, n=300, N=2000, q=2) # 0.9975 0.8135
#lineMultiSampleWatsonInferenceExperiment(kappa=3, n=300, N=2000, q=2) # 0.9945 0.9655
#lineMultiSampleWatsonInferenceExperiment(kappa=10, n=300, N=2000, q=2) # 0.9975 0.9555
#lineMultiSampleWatsonInferenceExperiment(kappa=100, n=300, N=2000, q=2) # 0.9965 0.9550
#lineMultiSampleWatsonInferenceExperiment(kappa=-100, n=1000, N=2000, q=2) # 0.9975 0.9440
#lineMultiSampleWatsonInferenceExperiment(kappa=-10, n=1000, N=2000, q=2) # 0.9990 0.9575
#lineMultiSampleWatsonInferenceExperiment(kappa=-3, n=1000, N=2000, q=2) # 0.998 0.936
#lineMultiSampleWatsonInferenceExperiment(kappa=-1, n=1000, N=2000, q=2) # 0.9955 0.7835
#lineMultiSampleWatsonInferenceExperiment(kappa=1, n=1000, N=2000, q=2) # 0.9960 0.8115
#lineMultiSampleWatsonInferenceExperiment(kappa=3, n=1000, N=2000, q=2) # 0.9970 0.9625
#lineMultiSampleWatsonInferenceExperiment(kappa=10, n=1000, N=2000, q=2) # 0.9995 0.9625
#lineMultiSampleWatsonInferenceExperiment(kappa=100, n=1000, N=2000, q=2) # 0.9955 0.9530
#lineMultiSampleWatsonInferenceExperiment(kappa=-100, n=3000, N=2000, q=2) # 0.9990 0.9495 36minutes
#lineMultiSampleWatsonInferenceExperiment(kappa=-10, n=3000, N=2000, q=2) # 0.9995 0.9540
#lineMultiSampleWatsonInferenceExperiment(kappa=-3, n=3000, N=2000, q=2) # 0.9985 0.9260
#lineMultiSampleWatsonInferenceExperiment(kappa=-1, n=3000, N=2000, q=2) # 0.9975 0.7820
#lineMultiSampleWatsonInferenceExperiment(kappa=1, n=3000, N=2000, q=2) # 0.9965 0.8040
#lineMultiSampleWatsonInferenceExperiment(kappa=3, n=3000, N=2000, q=2) # 0.9975 0.9615
#lineMultiSampleWatsonInferenceExperiment(kappa=10, n=3000, N=2000, q=2) # 0.9965 0.9555
#lineMultiSampleWatsonInferenceExperiment(kappa=100, n=3000, N=2000, q=2) # 0.9990 0.9535



### BINGHAM DISTRIBUTION ###

#' Simulation from the Bingham distribution.
#' 
#' Uses the function rbingham in package Directional. In that convention, the Bingham probability density is proportional to exp(-x^T A x), not exp(x^T A x). The mean is the eigenvector of A with least eigenvalue. The main direction of the dispersion is toward the eigenvector of A with the intermediate eigenvalue. The pole to the dispersion is the eigenvector with greatest eigenvalue.
#' @param a A symmetric real 3x3 matrix.
#' @param n A real number (positive integer) or NULL.
#' @return If n is NULL, then a single line. If n is a positive integer, then a list of n lines.
lineBingham <- function(a, n=NULL) {
  if (is.null(n))
    lower(as.numeric(rbingham(1, a)))
  else {
    us <- rbingham(n, a)
    if (n == 1)
      list(lower(as.numeric(us)))
    else
      lapply(1:nrow(us), function(i) lower(us[i,]))
  }
}

# Bingham MLE lookup table from Appendix C of Tauxe's books.
lineBinghamK1Table <- cbind(
  c(-25.55, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-25.56, -13.11, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-25.58, -13.14, -9.043, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-25.6, -13.16, -9.065, -7.035, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-25.62, -13.18, -9.080, -7.042, -5.797, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-25.63, -13.19, -9.087, -7.041, -5.789, -4.917, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-25.64, -13.20, -9.087, -7.033, -5.773, -4.896, -4.231, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-25.65, -13.20, -9.081, -7.019, -5.752, -4.868, -4.198, -3.659, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-25.65, -13.19, -9.068, -6.999, -5.726, -4.836, -4.160, -3.616, -3.160, NA, NA, NA, NA, NA, NA, NA),
  c(-25.64, -13.18, -9.05, -6.974, -5.694, -4.799, -4.118, -3.570, -3.109, -2.709, NA, NA, NA, NA, NA, NA),
  c(-25.63, -13.17, -9.027, -6.944, -5.658, -4.757, -4.071, -3.518, -3.053, -2.649, -2.289, NA, NA, NA, NA, NA),
  c(-25.61, -23.14, -8.999, -6.910, -5.618, -4.711, -4.021, -3.463, -2.993, -2.584, -2.220, -1.888, NA, NA, NA, NA),
  c(-25.59, -13.12, -8.966, -6.870, -5.573, -4.661, -3.965, -3.403, -2.928, -2.515, -2.146, -1.809, -1.497, NA, NA, NA),
  c(-25.57, -13.09, -8.928, -6.827, -5.523, -4.606, -3.906, -3.338, -2.859, -2.441, -2.066, -1.724, -1.406, -1.106, NA, NA),
  c(-25.54, -13.05, -8.886, -6.778, -5.469, -4.547, -3.842, -3.269, -2.785, -2.361, -1.981, -1.634, -1.309, -1.002, -0.708, NA),
  c(-25.50, -13.01, -8.839, -6.725, -5.411, -4.484, -3.773, -3.195, -2.706, -2.277, -1.891, -1.537, -1.206, -0.891, -0.588, -0.292),
  c(-25.46, -12.96, -8.788, -6.668, -5.348, -4.415, -3.699, -3.116, -2.621, -2.186, -1.794, -1.433, -1.094, -0.771, -0.459, -0.152),
  c(-25.42, -12.91, -8.731, -6.606, -5.280, -4.342, -3.620, -3.032, -2.531, -2.089, -1.690, -1.322, -0.974, -0.642, NA, NA),
  c(-25.37, -12.86, -8.670, -6.539, -5.207, -4.263, -3.536, -2.941, -2.434, -1.986, -1.579, -1.202, NA, NA, NA, NA),
  c(-25.31, -12.80, -8.604, -6.466, -5.126, -4.179, -3.446, -2.845, -2.330, -1.874, NA, NA, NA, NA, NA, NA),
  c(-25.5, -12.73, -8.532, -6.388, -5.045, -4.089, -3.349, -2.741, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-25.19, -12.66, -8.454, -6.305, -4.955, -3.992, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-25.12, -12.58, -8.371, -6.215, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))

# Bingham MLE lookup table from Appendix C of Tauxe's books.
lineBinghamK2Table <- cbind(
  c(-25.55, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-13.09, -13.11, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-8.996, -9.019, -9.043, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-6.977, -6.999, -7.020, -7.035, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-5.760, -5.777, -5.791, -5.798, -5.797, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-4.923, -4.934, -4.941, -4.941, -4.933, -4.917, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-4.295, -4.301, -4.301, -4.294, -4.279, -4.258, -4.231, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-3.796, -3.796, -3.790, -3.777, -3.756, -3.729, -3.697, -3.659, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-3.381, -3.375, -3.363, -3.345, -3.319, -3.287, -3.249, -3.207, -3.160, NA, NA, NA, NA, NA, NA, NA),
  c(-3.025, -3.014, -2.997, -2.973, -2.942, -2.905, -2.863, -2.816, -2.765, -2.709, NA, NA, NA, NA, NA, NA),
  c(-2.712, -2.695, -2.673, -2.644, -2.609, -2.568, -2.521, -2.470, -2.414, -2.354, -2.289, NA, NA, NA, NA, NA),
  c(-2.431, -2.410, -2.382, -2.349, -2.309, -2.263, -2.212, -2.157, -2.097, -2.032, -1.963, -1.888, NA, NA, NA, NA),
  c(-2.175, -2.149, -2.117, -2.078, -2.034, -1.984, -1.929, -1.869, -1.805, -1.735, -1.661, -1.582, -1.497, NA, NA, NA),
  c(-1.939, -1.908, -1.871, -1.828, -1.779, -1.725, -1.665, -1.601, -1.532, -1.458, -1.378, -1.294, -1.203, -1.106, NA, NA),
  c(-1.718, -1.682, -1.641, -1.596, -1.540, -1.481, -1.417, -1.348, -1.274, -1.195, -1.110, -1.020, -0.923, -0.819, -0.708, NA),
  c(-1.510, -1.470, -1.423, -1.371, -1.313, -1.250, -1.181, -1.108, -1.028, -0.944, -0.853, -0.756, -0.653, -0.541, -0.421, -0.292), 
  c(-1.312, -1.267, -1.216, -1.159, -1.096, -1.028, -0.955, -0.876, -0.791, -0.701, -0.604, -0.500, -0.389, -0.269, -0.140, 0.000),
  c(-1.123, -1.073, -1.017, -9.555, -0.887, -0.814, -0.736, -0.651, -0.561, -0.464, -0.360, -0.249, -0.129, 0.000, NA, NA),
  c(-0.940, -0.885, -0.824, -0.757, -0.684, -0.606, -0.522, -0.432, -0.335, -0.231, -0.120, 0.000, NA, NA, NA, NA),
  c(-0.762, -0.702, -0.636, -0.564, -0.486, -0.402, -0.312, -0.215, -0.111, -0.000, NA, NA, NA, NA, NA, NA),
  c(-0.589, -0.523, -0.452, -0.374, -0.290, -0.200, -0.104, 0.000, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-0.418, -0.347, -0.270, -0.186, -0.097, 0.000, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  c(-0.250, -0.173, -0.090, 0.000, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))

# Computes K1, K2 from omega1, omega2, using the lookup tables.
lineBinghamK1K2MLE <- function(omega1, omega2) {
  omega1s <- seq(from=0.02, to=0.32, by=0.02)
  omega2s <- seq(from=0.02, to=0.46, by=0.02)
  point <- matrix(c(omega1, omega2), 1, 2)
  k1 <- interp.surface(list(x=omega1s, y=omega2s, z=lineBinghamK1Table), point)
  k2 <- interp.surface(list(x=omega1s, y=omega2s, z=lineBinghamK2Table), point)
  c(k1, k2)
}

#' Maximum likelihood estimation of the Bingham distribution parameters.
#' 
#' Based on Tauxe (2010). May not work if the data are too concentrated or too dispersed? Results are great for n == 1000, so-so for n == 100, and poor for n == 10. The Bingham probability density is proportional to exp(-x^T A x), not exp(x^T A x). I prefer to normalize A so that tr A == 0.
#' @param us A list of lines.
#' @param weights A vector of real numbers (positive), of length equal to us. They need not sum to 1; the function automatically normalizes them to do so.
#' @return A list with members $a (symmetric 3x3 real matrix A), $values (a real 3D vector; the eigenvalues of A), and $vectors (a rotation matrix; the eigenvectors of A are the columns). The values are in descending order and sum to zero.
lineTauxeBinghamMLE <- function(us, weights=replicate(length(us), 1)) {
  # In eigen, the eigenvectors are descending and the eigenvectors are unit.
  #tBar <- arithmeticMean(lapply(us, function(u) outer(u, u)))
  w <- weights / sum(weights)
  ts <- lapply(1:length(us), function(i) {w[[i]] * outer(us[[i]], us[[i]])})
  tBar <- Reduce("+", ts)
  eig <- eigen(tBar, symmetric=TRUE)
  # Find MLE concentration parameters k1, k2, based on eigvals omega1 <= omega2 <= omega3.
  omega1 <- eig$values[[3]]
  omega2 <- eig$values[[2]]
  omega3 <- eig$values[[1]]
  k1k2 <- lineBinghamK1K2MLE(omega1, omega2)
  # Compute the eigenvalues of A. Normalize so that tr A == 0.
  values <- c(-k1k2[[1]], -k1k2[[2]], 0)
  values <- values - sum(values) / 3
  # Repackage into matrix A.
  vectors <- cbind(eig$vectors[,3], eig$vectors[,2], eig$vectors[,1])
  if (det(vectors) < 0)
    vectors[,3] <- -vectors[,3]
  a <- vectors %*% diag(values) %*% t(vectors)
  list(a=a, values=values, vectors=vectors)
}

# Test.
#vectors <- rotUniform()
#values <- c(8, 4, 0)
#a <- vectors %*% diag(values) %*% t(vectors)
#us <- lineBingham(a, n=1000)
#mle <- lineTauxeBinghamMLE(us)
#a
#mle$a
#vectors
#mle$vectors
#oriDistance(t(vectors), t(mle$vectors), group=oriLineInPlaneGroup)
#values - min(values)
#mle$values - min(mle$values)

#' 95% confidence region for the mean of the Bingham distribution.
#' 
#' This function sometimes fails, if the data set is too concentrated, dispersed, or small? This is not a particularly high-quality implementation of the technique.
#' @param ls A list of lines.
#' @return A list with members $directions, $scatter, and $angles095. The first two members are identical to $vectors and $values in lineProjectedMean. $angles095 is a pair of real numbers. They describe the 95% confidence region, as two distances from the mean toward the two other principal dispersions, measured in radians along the unit sphere's surface. If the inference fails, then the angles are NA.
lineBinghamInference095 <- function(us) {
  # Compute eigensystem of (1 / n) SUM u_i u_i^T. Tauxe leaves out the n.
  n <- length(us)
  oriTsr <- arithmeticMean(lapply(us, function(u) outer(u, u)))
  eig <- eigen(oriTsr, symmetric=TRUE)
  # Find MLE concentration parameters k1, k2, based on eigvals omega1 <= omega2 <= omega3.
  omega1 <- eig$values[[3]]
  omega2 <- eig$values[[2]]
  omega3 <- eig$values[[1]]
  k1k2 <- lineBinghamK1K2MLE(omega1, omega2)
  # Use simple 95% confidence expressions from Tauxe Appendix C. She has another n here.
  k1 <- k1k2[[1]]
  k2 <- k1k2[[2]]
  epsilon32 <- 1.22 / (k2 * (omega2 - omega3))
  epsilon31 <- 1.22 / (k1 * (omega1 - omega3))
  list(directions=eig$vectors,
       scatter=eig$values,
       angles095=c(epsilon32, epsilon31))
}



### REGRESSION ##

#' Fitting a great circle arc to a set of lines.
#' 
#' @param xs A vector of real numbers. Values of the independent scalar variable x.
#' @param ls A list of lines, of the same length as xs.
#' @param numSteps A real number (positive integer). A bound on the number of numerical optimization iterations allowed. If you find that the results have $error != 0, then try increasing numSteps.
#' @param numPoints A real number (positive integer). The number of points along the regression curve requested.
#' @return A list (a, rotation, error, minEigenvalue, rSquared, prediction). The great circle curve is v(x) = R [cos(a x) sin(a x) 0]^T, where R == rotation. The third column of R is the pole of the great circle. a is the angle of rotation about that pole per unit x, in radians. So a / degree is the angle in degrees. error should be 0 and minEigenvalue should be positive. Otherwise there was some problem in the optimization. rSquared is the squared correlation coefficient, 0 <= R^2 <= 1. If numPoints > 0 then the return list also has a member $points consisting of numPoints+1 lines. prediction is an R function that takes an x value (real number) as input and produces the line v(x) as output.
lineGeodesicRegression <- function(xs, ls, numSteps=1000, numPoints=0) {
  # Let l0 be the l whose x is closest to zero.
  l0 <- ls[[which.min(as.numeric(xs)^2)]]
  # Define the function to be minimized.
  n <- length(ls)
  e <- function(wb) {
    a <- rotExp(rotAntisymmetricFromVector(wb[1:3]))
    f <- function(i) {
      lineDistance(ls[[i]], as.numeric(a %*% c(cos(wb[[4]] * xs[[i]]), sin(wb[[4]] * xs[[i]]), 0)))^2
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
  rSq <- 1 - solution$value / lineVariance(ls, lineProjectedMean(ls))
  result <- list(a=a, rotation=rot, error=solution$convergence, minEigenvalue=min(eigvals), rSquared=rSq)
  result$prediction <- function(x) {
    as.numeric(rot %*% c(cos(a * x), sin(a * x), 0))
  }
  if (numPoints >= 1)
    result$points <- lapply(seq(from=min(xs), to=max(xs), length.out=(numPoints + 1)), result$prediction)
  result
}

#' Fitting a great circle arc to a set of lines.
#' 
#' Similar to regressionLineGeodesicCurve, but internally rescales x to [0, 1] for better performance. Try this function first, and use the other function only if you find reason to.
#' @param xs See lineGeodesicRegression.
#' @param ls See lineGeodesicRegression.
#' @param numSteps See lineGeodesicRegression.
#' @param numPoints See lineGeodesicRegression.
#' @return See lineGeodesicRegression.
lineRescaledGeodesicRegression <- function(xs, ls, numSteps=1000, numPoints=0) {
  x0 <- min(xs)
  x1 <- max(xs)
  # Perform regression in scaled coordinates.
  regr <- lineGeodesicRegression(scales(xs), ls, numSteps, numPoints=0)
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

# Returns up to numPerms R^2 values, for comparison to the original.
#' Permutation test for geodesic regression of lines.
#' 
#' You request some number of permutations. This function attempts that many. Sometimes some of them fail, due to the optimization process. If you're getting many failures, then consider increasing numSteps. Anyway, after computing these R^2 values, you will want to know how many there are and how many exceed the R^2 for the original regression. The ratio of the latter to the former is your p-value for the null hypothesis that there is no geodesic dependence of the ls on the xs. In other words, small values of p indicate that the regression is meaningful.
#' @param xs See lineGeodesicRegression.
#' @param ls See lineGeodesicRegression.
#' @param numPerms A real number (positive integer). The number of permutations to try, say 1,000 or 10,000.
#' @param numSteps See lineGeodesicRegression.
#' @return A vector of real numbers, all between 0 and 1 inclusive. The length (dimension) of this vector cannot exceed numPerms. The R^2 value for each of the permutation tests in which the optimization succeeded.
lineGeodesicRegressionPermutations <- function(xs, ls, numPerms, numSteps=10000) {
  ys <- scales(xs)
  f <- function(i) {
    print(i / numPerms)
    regr <- lineGeodesicRegression(sample(ys, size=length(ys)), ls, numSteps)
    c(regr$error, regr$minEigenvalue, regr$rSquared)
  }
  perms <- sapply(1:numPerms, f)
  perms[3,][perms[1,] == 0 & perms[2,] > 0]
}

lineSmallCircleRegression <- function(xs, us, numSeeds=5, numSteps=1000, numPoints=0) {
  f <- function(phiThetaAlphaTauSigma) {
    pole <- cartesianFromSpherical(c(1, phiThetaAlphaTauSigma[1:2]))
    uOf0 <- cartesianFromSpherical(c(1, phiThetaAlphaTauSigma[4:5]))
    pred <- function(x) {
      angle <- x * phiThetaAlphaTauSigma[[3]]
      as.numeric(rotMatrixFromAxisAngle(c(pole, angle)) %*% uOf0)
    }
    preds <- lapply(xs, pred)
    dists <- mapply(lineDistance, preds, us)
    dot(dists, dists) / (2 * length(us))
  }
  # Find the minimum, starting from a few seeds.
  best <- list(value=(2 * pi^2))
  for (i in 1:numSeeds) {
    seed <- c(sphericalFromCartesian(lineUniform())[2:3], runif(1, -pi, pi), 
              sphericalFromCartesian(lineUniform())[2:3])
    solution <- optim(seed, f, lower=c(0, -pi, -pi, 0, -pi), upper=c(pi, pi, pi, pi, pi),
                      hessian=TRUE, control=list(maxit=numSteps), method="L-BFGS-B")
    if (solution$value <= best$value)
      best <- solution
  }
  # Report results.
  eigvals <- eigen(best$hessian, symmetric=TRUE, only.values=TRUE)$values
  pole <- cartesianFromSpherical(c(1, best$par[1:2]))
  angle <- best$par[[3]]
  uOf0 <- cartesianFromSpherical(c(1, best$par[4:5]))
  var <- lineVariance(us, lineProjectedMean(us))
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

lineRescaledSmallCircleRegression <- function(xs, ls, ...) {
  # Perform regression in scaled coordinates.
  x0 <- min(xs)
  x1 <- max(xs)
  results <- lineSmallCircleRegression(scales(xs), ls, ...)
  # Scale the results back into the original coordinates.
  results$rescaledPrediction <- results$prediction
  results$prediction <- function(x) {results$rescaledPrediction((x - x0) / (x1 - x0))}
  results$angle <- results$angle / (x1 - x0)
  results
}

lineSmallCircleRegressionPermutations <- function(xs, ls, numPerms, ...) {
  ys <- scales(xs)
  f <- function(i) {
    print(i / numPerms)
    regr <- lineSmallCircleRegression(sample(ys, size=length(ys)), ls, ...)
    c(regr$error, regr$minEigenvalue, regr$rSquared)
  }
  perms <- sapply(1:numPerms, f)
  perms[3,][perms[1,] == 0 & perms[2,] > 0]
}



### WELLNER (1979) TWO-SAMPLE TEST ###

#' Wellner's T-statistic (Wellner, 1979), which measures how different two sets of lines are.
#' 
#' @param xs A list of lines.
#' @param ys A list of lines.
#' @return A real number (nonnegative). Zero if the two data sets are identical.
lineWellner <- function(xs, ys) {
    # Prepare the number of xs, the number of ys, and the dimension of the hypersphere.
    m <- length(xs)
    n <- length(ys)
    p <- length(xs[[1]]) - 1
    # Compute m Tx and n Ty.
    f <- function(v) {outer(v, v)}
    mTx <- Reduce("+", lapply(xs, f))
    nTy <- Reduce("+", lapply(ys, f))
    # Compute Wellner's statistic.
    txMinusTy <- mTx / m - nTy / n
    tr(txMinusTy %*% txMinusTy) * m * n * (p + 1) * (p + 3) / (2 * (m + n))
}

# Also computes Wellner's T-statistic, with inputs presented in a manner suitable for fast permutation testing.
lineWellnerShortcut <- function(xys, choices, mTxPlusnTy) {
    # Prepare the number of xs, the number of ys, and the dimension of the hypersphere.
    m <- length(choices)
    n <- length(xys) - m
    p <- length(xys[[1]]) - 1
    # Compute m Tx and n Ty.
    mTx <- Reduce("+", lapply(choices, function(i) {outer(xys[[i]], xys[[i]])}))
    nTy <- mTxPlusnTy - mTx
    # Compute Wellner's statistic.
    txMinusTy <- mTx / m - nTy / n
    tr(txMinusTy %*% txMinusTy) * m * n * (p + 1) * (p + 3) / (2 * (m + n))
}

#' Two-sample test, based on permutations and Wellner's T-statistic (Wellner, 1979).
#' 
#' @param xs A list of lines.
#' @param ys A list of lines.
#' @param numPerms A real number (positive integer). The number of permutations, say 1,000 or 10,000.
#' @return A real number, between 0 and 1 inclusive. The fraction of tests in which T exceeds the original T for the data. You can interpret this as a p-value for the null hypothesis that the two populations are identical (not just that their means are identical). In other words, small values of p indicate that the distinction between the two populations is meaningful.
lineWellnerInference <- function(xs, ys, numPerms) {
    # Precompute.
    m <- length(xs)
    n <- length(ys)
    xys <- c(xs, ys)
    mTxPlusnTy <- Reduce("+", lapply(xys, function(v) {outer(v, v)}))
    # Compute Wellner's statistic for the actual data.
    t <- lineWellnerShortcut(xys, 1:m, mTxPlusnTy)
    # Compute Wellner's statistic for permuted data, numPerms times.
    ts <- replicate(numPerms, lineWellnerShortcut(xys, sample.int(m + n, m), mTxPlusnTy))
    # What proportion of permuted results are more extreme than the one observed?
    greaterThan <- Filter(function (u) {u > t}, ts)
    length(greaterThan) / numPerms
}

# Tests.
lineWellnerInferenceExperimentFuncAs <- function(eastValues, westValues, sameCenter, sameDirection) {
  function() {
    # Choose the centers. If sameDirection, then force sameCenter.
    east3 <- rayUniform()
    if (sameCenter || sameDirection)
      west3 <- east3
    else
      west3 <- rayUniform()
    # Choose the directions.
    east1 <- rayOrthogonalUniform(east3)
    if (sameDirection)
      west1 <- east1
    else
      west1 <- rayOrthogonalUniform(west3)
    eastRot <- cbind(east1, cross(east3, east1), east3)
    westRot <- cbind(west1, cross(west3, west1), west3)
    eastA <- eastRot %*% diag(eastValues) %*% t(eastRot)
    westA <- westRot %*% diag(westValues) %*% t(westRot)
    list(eastA, westA)
  }
}
lineWellnerInferenceExperiment <- function(funcAs, eastN, westN, numPerms, numTrials) {
  f <- function(i) {
    print(i / numTrials)
    as <- funcAs()
    easts <- lineBingham(as[[1]], eastN)
    wests <- lineBingham(as[[2]], westN)
    lineWellnerInference(easts, wests, numPerms=numPerms)
  }
  ps <- sapply(1:numTrials, f)
  pHat <- sum(ps < 0.05) / numTrials
  c(pHat, 2 * standardErrorProportion(numTrials, pHat))
}
#lineWellnerInferenceExperiment(lineWellnerInferenceExperimentFuncAs(c(1, 0, -1), c(1, 0, -1), TRUE, TRUE),
#                               10, 10, numPerms=1000, numTrials=2000) # 0.05450000 +- 0.01015182
#lineWellnerInferenceExperiment(lineWellnerInferenceExperimentFuncAs(c(2, 0, -2), c(2, 0, -2), TRUE, TRUE),
#                               10, 10, numPerms=1000, numTrials=2000) # 0.048000000 +- 0.009559916
#lineWellnerInferenceExperiment(lineWellnerInferenceExperimentFuncAs(c(5, 0, -5), c(5, 0, -5), TRUE, TRUE),
#                               10, 10, numPerms=1000, numTrials=2000) # 0.055 +- 0.01019559
#lineWellnerInferenceExperiment(lineWellnerInferenceExperimentFuncAs(c(5, 0, -1), c(5, 0, -1), TRUE, TRUE),
#                               10, 10, numPerms=1000, numTrials=2000) # 0.040500000 +- 0.008815866
#lineWellnerInferenceExperiment(lineWellnerInferenceExperimentFuncAs(c(1, 0, -5), c(1, 0, -5), TRUE, TRUE),
#                               10, 10, numPerms=1000, numTrials=2000) # 0.051500000 +- 0.009884103
#lineWellnerInferenceExperiment(lineWellnerInferenceExperimentFuncAs(c(1, 0, -1), c(1, 0, -1), TRUE, TRUE),
#                               100, 10, numPerms=1000, numTrials=2000) # 
#lineWellnerInferenceExperiment(lineWellnerInferenceExperimentFuncAs(c(2, 0, -2), c(2, 0, -2), TRUE, TRUE),
#                               100, 10, numPerms=1000, numTrials=2000) # 
#lineWellnerInferenceExperiment(lineWellnerInferenceExperimentFuncAs(c(5, 0, -5), c(5, 0, -5), TRUE, TRUE),
#                               100, 10, numPerms=1000, numTrials=2000) # 
#lineWellnerInferenceExperiment(lineWellnerInferenceExperimentFuncAs(c(5, 0, -1), c(5, 0, -1), TRUE, TRUE),
#                               100, 10, numPerms=1000, numTrials=2000) # 0.051500000 +- 0.009884103
#lineWellnerInferenceExperiment(lineWellnerInferenceExperimentFuncAs(c(1, 0, -5), c(1, 0, -5), TRUE, TRUE),
#                               100, 10, numPerms=1000, numTrials=2000) # 
#lineWellnerInferenceExperiment(lineWellnerInferenceExperimentFuncAs(c(1, 0, -1), c(1, 0, -1), TRUE, TRUE),
#                               100, 100, numPerms=1000, numTrials=2000) # 0.052500000 0.009974342
#lineWellnerInferenceExperiment(lineWellnerInferenceExperimentFuncAs(c(2, 0, -2), c(2, 0, -2), TRUE, TRUE),
#                               100, 100, numPerms=1000, numTrials=2000) # 0.052500000 0.009974342
#lineWellnerInferenceExperiment(lineWellnerInferenceExperimentFuncAs(c(5, 0, -5), c(5, 0, -5), TRUE, TRUE),
#                               100, 100, numPerms=1000, numTrials=2000) # 
#lineWellnerInferenceExperiment(lineWellnerInferenceExperimentFuncAs(c(5, 0, -1), c(5, 0, -1), TRUE, TRUE),
#                               100, 100, numPerms=1000, numTrials=2000) # 0.044000000 +- 0.009172132
#lineWellnerInferenceExperiment(lineWellnerInferenceExperimentFuncAs(c(1, 0, -5), c(1, 0, -5), TRUE, TRUE),
#                               100, 100, numPerms=1000, numTrials=2000) # 0.048000000 0.009559916



### PLOTTING FUNCTIONS ###

# Given 0 <= s <= 1, returns the line that is fraction s of the way from u to v. If u and v are perpendicular, then the ambiguity is resolved by treating them like rays.
lineInterpolation <- function(u, v, s=0.5) {
  if (dot(u, v) < 0)
    rayInterpolation(u, -v, s)
  else
    rayInterpolation(u, v, s)
}

# Returns a list of five boxes. Each box is a list of four vertices, in counter-clockwise order when viewed from above the sphere. Each vertex is a vector of four numbers: a unit 3D vector, with the value of f at that vector appended.
lineFunctionPlotBase <- function(f, ...) {
  a <- sqrt(2) / 2
  us <- list(c(a, 0, -a), c(0, a, -a), c(-a, 0, -a), c(0, -a, -a), c(1, 0, 0), c(0, 1, 0), c(-1, 0, 0), c(0, -1, 0))
  us <- lapply(us, function(u) c(u, f(u, ...)))
  list(
    list(us[[1]], us[[2]], us[[3]], us[[4]]),
    list(us[[1]], us[[5]], us[[6]], us[[2]]),
    list(us[[2]], us[[6]], us[[7]], us[[3]]),
    list(us[[3]], us[[7]], us[[8]], us[[4]]),
    list(us[[4]], us[[8]], us[[5]], us[[1]]))
}

# Given one box, returns a list of four boxes.
lineFunctionPlotRefinement <- function(box, f, ...) {
  mids <- lapply(1:4, function(i) lineInterpolation(box[[i]][1:3], box[[i %% 4 + 1]][1:3]))
  center <- lineInterpolation(lineInterpolation(mids[[1]], mids[[3]]), lineInterpolation(mids[[2]], mids[[4]]))
  mids <- lapply(mids, function(u) c(u, f(u, ...)))
  center <- c(center, f(center, ...))
  list(
    list(box[[1]], mids[[1]], center, mids[[4]]),
    list(box[[2]], mids[[2]], center, mids[[1]]),
    list(box[[3]], mids[[3]], center, mids[[2]]),
    list(box[[4]], mids[[4]], center, mids[[3]]))
}

# Given one box, returns a list of lines (pairs of unit vectors).
lineFunctionPlotCuts <- function(box, level) {
  cuts <- list()
  for (i in 1:4) {
    a <- box[[i]]
    b <- box[[i %% 4 + 1]]
    if (a[[4]] == level)
      cuts <- c(cuts, list(a[1:3]))
    else if (a[[4]] < level && b[[4]] > level)
      cuts <- c(cuts, list(lineInterpolation(a[1:3], b[1:3], (level - a[[4]]) / (b[[4]] - a[[4]]))))
    else if (b[[4]] < level && a[[4]] > level)
      cuts <- c(cuts, list(lineInterpolation(b[1:3], a[1:3], (level - b[[4]]) / (a[[4]] - b[[4]]))))
  }
  allPairs(cuts)
}

# Returns a list of length equal to that of levels. Each item in the list is a list of lines. Each line is a pair of unit vectors. Occasionally little 'stars' (complete graphs) may appear along a contour. If so, the non-adaptive refinement was not fine enough to resolve the contour in that box of the plot. So increase numNonAdapt.
lineFunctionPlotLineSets <- function(f, levels, numNonAdapt, ...) {
  # Non-adaptive refinement.
  boxes <- lineFunctionPlotBase(f, ...)
  if (numNonAdapt >= 1)
    for (i in 1:numNonAdapt) 
      boxes <- unlist(lapply(boxes, lineFunctionPlotRefinement, f, ...), recursive=FALSE)
  # Generate lines for this level of f.
  lineSets <- lapply(levels, function(level) unlist(lapply(boxes, lineFunctionPlotCuts, level), recursive=FALSE))
  lineSets
}

# lineFunctionPlotLineSets customized to the problem of Kamb contouring.
lineKambLineSets <- function(points, multiples=c(3, 6, 9, 12), k=3, degree=3, numNonAdapt=4) {
  # Compute the basic Kamb parameters.
  n <- length(points)
  p <- k^2 / (n + k^2)
  sigma <- squareRoot(n * p * (1 - p))
  r <- arcCos(1 - p)
  levels <- multiples * sigma
  # Compute the weighting function.
  if (degree == 0) {
    h <- function(p, u) {
      alpha <- lineDistance(p, u)
      if (alpha > r)
        0
      else
        1
    }
  } else {
    d0 <- 1 - cos(r)
    d1 <- -r * cos(r) + sin(r)
    d2 <- -2 + 2 * cos(r) - r^2 * cos(r) + 2 * r * sin(r)
    d3 <- 6 * r * cos(r) - r^3 * cos(r) - 6 * sin(r) + 3 * r^2 * sin(r)
    if (degree == 1)
      mat <- rbind(c(d0, d1, d2, d3), c(1, r, 0, 0), c(0, 0, 1, 0), c(0, 0, 0, 1))
    else
      # Default to degree-3 weighting polynomial.
      mat <- rbind(c(d0, d1, d2, d3), c(0, 1, 0, 0), c(0, 1, 2 * r, 3 * r^2), c(1, r, r^2, r^3))
    coeffs <- solve(mat, c(1 - cos(r), 0, 0, 0))
    h <- function(p, u) {
      alpha <- lineDistance(p, u)
      if (alpha > r)
        0
      else
        coeffs[[1]] + coeffs[[2]] * alpha + coeffs[[3]] * alpha^2 + coeffs[[4]] * alpha^3
    }
  }
  # Contour-plot the function.
  f <- function(u) sum(sapply(points, h, u))
  lineFunctionPlotLineSets(f, levels, numNonAdapt)
}



### PLOTTING ###

#' Equal-area, lower-hemisphere plot of lines.
#' 
#' This function essentially sends everything to the lower hemisphere and then invokes rayEqualAreaPlot. See rayEqualAreaPlot for details.
#' @param points A list of rays.
#' @param curves A list of lists of rays.
#' @param colors A character vector.
#' @param shapes A character vector. Traditionally lineations are squares, intermediates are triangles, and poles to foliation are circles.
#' @return NULL
lineEqualAreaPlot <- function(points=list(), curves=list(), colors=c("black"), shapes=c("c")) {
  # Convert the chosen shapes into their underlying shape codes.
  f <- function(i) {
    j <- (i - 1) %% length(shapes) + 1
    if (shapes[[j]] == "s")
      15
    else if (shapes[[j]] == "t")
      17
    else if (shapes[[j]] == "c")
      19
    else
      46
  }
  underShapes <- sapply(1:length(points), f)
  # Send the points to the lower hemisphere.
  underPoints <- lapply(points, lower)
  # Break the curves and send them to the lower hemisphere.
  underCurves <- list()
  for (curve in curves) {
    curvesSigns <- rayCurvesUpperLower(curve)
    for (i in 1:length(curvesSigns$curves)) {
      if (curvesSigns$signs[[i]] == 1)
        underCurves[[length(underCurves) + 1]] <- lapply(curvesSigns$curves[[i]], function(v) -v)
      else
        underCurves[[length(underCurves) + 1]] <- curvesSigns$curves[[i]]
    }
  }
  plotEqualArea(points=underPoints, curves=underCurves, colors=colors, shapes=underShapes)
}

#' Equal-area, lower-hemisphere plot of lines, with a circle about each point.
#' 
#' This function essentially sends everything to the lower hemisphere and then invokes rayEqualAreaRadiusPlot.
#' @param points A list of lines, to be plotted as points.
#' @param radii A vector of real numbers. The radii of circles to be drawn about the points. Radii are measured in radians, along the surface of the unit sphere.
#' @param colors A character vector. See rayEqualAreaPlot.
#' @param shapes A character vector. See rayEqualAreaPlot.
#' @return NULL
lineEqualAreaRadiusPlot <- function(points=list(), radii=c(), colors=c("black"), shapes=c("c")) {
  lineEqualAreaPlot(points, colors=colors, curves=thread(raySmallCircle, points, radii))
}

#' Equal-area, lower-hemisphere plot of lines, with contours representing level sets of an arbitrary given function.
#' 
#' @param f An R function. Its input is a single line u, and optionally additional parameters '...'. Its output is a real number f(u, ...).
#' @param levels A vector of real numbers. The numbers c such that the contours f(u, ...) = c are plotted.
#' @param numNonAdapt A real number (non-negative integer). The number of non-adaptive refinements to make. These control the quality of the plot. Each increment to numNonAdapt makes the contours smoother, and helps resolve 'star defects' that sometimes appear along the contours. On the other hand, each increment increases the time and memory requirements by a factor of four.
#' @param points A list of lines, to be plotted as points. Need not have anything to do with f.
#' @param colors A character vector. See rayEqualAreaPlot.
#' @param shapes A character vector. See rayEqualAreaPlot.
#' @param ... Additional arguments to be passed to f, other than u.
#' @return NULL
lineEqualAreaFunctionPlot <- function(f, levels, numNonAdapt=4, points=list(), colors=c("black"), shapes=c("c"), ...) {
  lineSets <- lineFunctionPlotLineSets(f, levels, numNonAdapt, ...)
  lineEqualAreaPlot(points=points, curves=unlist(lineSets, recursive=FALSE), colors=colors, shapes=shapes)
}

#' Equal-area, lower-hemisphere plot of lines, with Kamb contours representing density.
#' 
#' @param points A list of lines, to be plotted as points.
#' @param multiples A vector of real numbers (positive). The multiples of sigma to be plotted as contours.
#' @param k A real number (positive). The arbitrary smoothing factor of Kamb (1959). I see no practical reason to mess with this, other than curiosity.
#' @param degree A real number (0, 1, or 3). The degree of the weighting polynomial. Higher degrees tend to produce smoother plots. I see no practical reason to mess with this, other than curiosity.
#' @param numNonAdapt A real number (non-negative integer). The number of non-adaptive refinements to make. These control the quality of the plot. Each increment to numNonAdapt makes the contours smoother, and helps resolve 'star defects' that sometimes appear along the contours. On the other hand, each increment increases the time and memory requirements by a factor of four.
#' @param colors A character vector. See rayEqualAreaPlot.
#' @param shapes A character vector. See rayEqualAreaPlot.
#' @return NULL
lineKambPlot <- function(points, multiples=c(3, 6, 9, 12), k=3, degree=3, numNonAdapt=4, colors=c("black"), shapes=c("c")) {
  lineSets <- lineKambLineSets(points=points, multiples=multiples, k=k, degree=degree, numNonAdapt=numNonAdapt)
  lineEqualAreaPlot(points=points, curves=unlist(lineSets, recursive=FALSE), colors=colors, shapes=shapes)
}

#' Convenience shortcut for plotting two sets of lines.
#' @param pointsA A list of lines.
#' @param pointsB A list of lines.
#' @param colorA Character. A color, of the sort used in all R graphics routines.
#' @param colorB Character. A color, of the sort used in all R graphics routines.
#' @param curves A list of lists of lines. Curves to be plotted.
#' @return NULL.
lineEqualAreaPlotTwo <- function(pointsA, pointsB, colorA="red", colorB="cyan", curves=list()) {
  lineEqualAreaPlot(points=c(pointsA, pointsB), curves=curves,
                    colors=c(replicate(length(pointsA), colorA), replicate(length(pointsB), colorB)))
}

#' Convenience shortcut for plotting three sets of lines.
#' @param pointsA A list of lines.
#' @param pointsB A list of lines.
#' @param pointsC A list of lines.
#' @param colorA Character. A color, of the sort used in all R graphics routines.
#' @param colorB Character. A color, of the sort used in all R graphics routines.
#' @param colorC Character. A color, of the sort used in all R graphics routines.
#' @param curves A list of lists of lines. Curves to be plotted.
#' @return NULL.
lineEqualAreaPlotThree <- function(pointsA, pointsB, pointsC, colorA="red", colorB="green", colorC="blue", curves=list()) {
  lineEqualAreaPlot(points=c(pointsA, pointsB, pointsC), curves=curves,
                    colors=c(replicate(length(pointsA), colorA),
                             replicate(length(pointsB), colorB),
                             replicate(length(pointsC), colorC)))
}

#' Equal-area plot of lines, extruded into the third dimension arbitrarily by the user.
#' 
#' @param vs A list of lines.
#' @param zs A vector of real numbers, of the same length as vs. The third coordinate to be plotted.
#' @param ... Other plotting options to be passed to the underlying plot3D.
#' @return NULL.
lineEqualAreaScalarPlot <- function(vs, zs, ...) {
  # Make the points.
  xyzs <- thread(function(v, z) c(equalAreaProjection(lower(v)), z), vs, zs)
  # Make the boundary circle a few times.
  radius <- sqrt(2)
  xys <- lapply((0:72) * 2 * pi / 72, function(theta) {radius * c(cos(theta), sin(theta))})
  first <- lapply(xys, function(xy) c(xy, 0))
  second <- lapply(xys, function(xy) c(xy, min(zs)))
  third <- lapply(xys, function(xy) c(xy, max(zs)))
  plot3D(points=xyzs, curves=list(first, second, third), ...)
}



### TWO DIMENSIONS ###

# Some of the machinery above still works: lineMeanScatter, lineProjectedMean, lineDistance. Even lineWellner and lineWellnerInference.

#' Uniformly random 2D lines.
#' 
#' @param n A real number (positive integer) or NULL.
#' @return If n is NULL, then a single 2D line. If n is a positive integer, then a list of n 2D lines.
line2DUniform <- ray2DUniform

#' Random 2D lines, drawn from the wrapped normal distribution on the unit circle.
#' 
#' @param mean A 2D line.
#' @param sd A real number (positive). The standard deviation sigma of the underlying normal distribution, in radians.
#' @param n A real number (positive integer) or NULL.
#' @return If n is NULL, then a single 2D line. If n is a positive integer, then a list of n 2D lines.
line2DWrappedNormal <- ray2DWrappedNormal

#' Bootstrapped extrinsic mean for 2D lines.
#' 
#' Essentially ray2DBootstrapInference, adapted for lines.
#' @param ls A list of 2D lines.
#' @param numBoots A real number (positive integer). The number of bootstrapped means to compute. 10,000 might be a good value.
#' @return A list ($center, $bootstraps, $pvalue, $q000, $q025, $q050, $q075, $q095, $q099, $q100). $bootstraps are the bootstraps. $center is their mean. The other fields are quantiles of distance from the mean, in radians, among the bootstraps. For example, a 95% confidence region consists of all 2D rays within $q095 of $center. $pvalue is an R function from 2D rays to real numbers, assigning a p-value to any given null hypothesis for the mean.
line2DBootstrapInference <- function(ls, numBoots) {
  boots <- replicate(numBoots, lineProjectedMean(sample(ls, length(ls), replace=TRUE)), simplify=FALSE)
  bootMean <- lineProjectedMean(boots)
  dists <- sapply(boots, lineDistance, bootMean)
  empiricalCDF <- ecdf(dists)
  # Build the p-value function.
  f <- function(u) {
    1 - empiricalCDF(lineDistance(u, bootMean))
  }
  # Compute a few popular percentiles.
  qs <- quantile(dists, probs=c(0.00, 0.25, 0.50, 0.75, 0.95, 0.99, 1.00), names=FALSE)
  list(center=bootMean, bootstraps=boots, pvalue=f,
       q000=qs[[1]], q025=qs[[2]], q050=qs[[3]], q075=qs[[4]], q095=qs[[5]], q099=qs[[6]], q100=qs[[7]])
}

#' Rose plot for 2D lines.
#' 
#' Nearly identical to ray2DRosePlot, but with each line represented by two rays. In theory this rose diagram should be perfectly symmetric about the origin. In practice, data points often fall on the boundary between bins and are resolved arbitrarily, producing some leakage from bins into neighboring bins and hence asymmetry?
#' @param lines A list of 2D or 3D real vectors. In each one, only the 2D projection (the first two components) is used. It must be non-zero but need not be unit.
#' @param weights A vector of real numbers. The weights to attach to the lines.
#' @param ... Other options to be passed to the underlying rayRosePlot.
#' @return NULL.
line2DRosePlot <- function(lines, weights=replicate(length(lines), 1), ...) {
  rays <- c(lines, lapply(lines, function(v) -v))
  wghts <- c(weights / 2, weights / 2)
  ray2DRosePlot(rays, wghts, ...)
}


