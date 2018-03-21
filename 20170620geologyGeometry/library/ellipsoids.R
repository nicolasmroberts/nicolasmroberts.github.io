


# Copyright 2016-2017 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# An ellipsoid is a complete description of an ellipsoid in three-dimensional space: its size, shape, and orientation (but not its location). There are three degrees of freedom in orientation and three degrees of freedom in size-shape. Sometimes we normalize ellipsoids to have a particular volume, in which case only two degrees of freedom remain in shape. An ellipsoid can be described as a symmetric, positive-definite ellipsoid tensor E, in that the (boundary of the) ellipsoid is the set of points x such that x^T E x == 1. Or it can be described as an orientation and semi-axis lengths, or as a log-ellipsoid vector. In this R code, an 'ellipsoid' is a list of five elements: $vector, $tensor, $a, $logA, $rotation. $a are the semi-axis lengths, and $logA is their logarithms. $rotation is a special orthogonal matrix with the semi-axes of the ellipsoid along its rows. See also geoDataFromFile. Because they inhabit a vector space, we can throw all of multivariate statistics at the ellipsoid vectors. For ideas, see http://cran.r-project.org/web/views/Multivariate.html.



### CONVERSIONS AMONG REPRESENTATIONS ###

#' A log-ellipsoid tensor from an ellipsoid tensor, respecting normalization.
#' 
#' @param ell A real 3x3 matrix (symmetric, positive-definite).
#' @param A real 3x3 matrix (symmetric).
ellLog <- function(ell) {
  eig <- eigen(ell, symmetric=TRUE)
  eig$vectors %*% diag(log(eig$values)) %*% t(eig$vectors)
}

#' An ellipsoid tensor from a log-ellipsoid tensor, respecting normalization.
#' 
#' @param logEll A real 3x3 matrix (symmetric).
#' @param A real 3x3 matrix (symmetric, positive-definite).
ellExp <- function(logEll) {
  eig <- eigen(logEll, symmetric=TRUE)
  eig$vectors %*% diag(exp(eig$values)) %*% t(eig$vectors)
}

#' An ellipsoid 6-vector from an unnormalized log-ellipsoid tensor.
#' 
#' The conversion is tuned so that the Frobenius inner product on matrices maps to the dot product on vectors.
#' @param logEll A real 3x3 matrix (symmetric).
#' @return A 6-dimensional real vector.
ellVectorFromLog <- function(logEll) {
  v1 <- sqrt(2) * logEll[1, 2]
  v2 <- sqrt(2) * logEll[1, 3]
  v3 <- sqrt(2) * logEll[2, 3]
  v4 <- logEll[1, 1]
  v5 <- logEll[2, 2]
  v6 <- logEll[3, 3]
  c(v1, v2, v3, v4, v5, v6)
}

#' An unnormalized log-ellipsoid tensor from an ellipsoid 6-vector.
#' 
#' The conversion is tuned so that the Frobenius inner product on matrices maps to the dot product on vectors.
#' @param vec A 6-dimensional real vector.
#' @return A real 3x3 matrix (symmetric).
ellLogFromVector <- function(vec) {
  l12 <- vec[[1]] / sqrt(2)
  l13 <- vec[[2]] / sqrt(2)
  l23 <- vec[[3]] / sqrt(2)
  l11 <- vec[[4]]
  l22 <- vec[[5]]
  l33 <- vec[[6]]
  matrix(c(l11, l12, l13, l12, l22, l23, l13, l23, l33), 3, 3)
}

#' An ellipsoid 5-vector from a normalized log-ellipsoid tensor.
#' 
#' The conversion is tuned so that the Frobenius inner product on matrices maps to the dot product on vectors.
#' @param logEll A real 3x3 matrix (symmetric, trace-zero).
#' @return A 5-dimensional real vector.
ellNormalizedVectorFromLog <- function(normLogEll) {
  v1 <- sqrt(2) * normLogEll[1, 2]
  v2 <- sqrt(2) * normLogEll[1, 3]
  v3 <- sqrt(2) * normLogEll[2, 3]
  v4 <- (normLogEll[2, 2] + normLogEll[1, 1]) * sqrt(1.5)
  v5 <- (normLogEll[2, 2] - normLogEll[1, 1]) / sqrt(2)
  c(v1, v2, v3, v4, v5)
}

#' A normalized log-ellipsoid tensor from an ellipsoid 5-vector.
#' 
#' The conversion is tuned so that the Frobenius inner product on matrices maps to the dot product on vectors.
#' @param vec A 5-dimensional real vector.
#' @return A real 3x3 matrix (symmetric, trace-zero).
ellLogFromNormalizedVector <- function(vec) {
  l11 <- vec[[4]] / sqrt(6) - vec[[5]] / sqrt(2)
  l22 <- vec[[4]] / sqrt(6) + vec[[5]] / sqrt(2)
  l33 <- -sqrt(2 / 3) * vec[[4]]
  l12 <- vec[[1]] / sqrt(2)
  l13 <- vec[[2]] / sqrt(2)
  l23 <- vec[[3]] / sqrt(2)
  matrix(c(l11, l12, l13, l12, l22, l23, l13, l23, l33), 3, 3)
}

#' An ellipsoid from an ellipsoid vector.
#' 
#' @param v An ellipsoid vector, either 5-dimensional (if normalized) or 6-dimensional (if not).
#' @return An ellipsoid.
ellEllipsoidFromVector <- function(v) {
  # Diagonalize the log ellipsoid tensor.
  if (length(v) == 5)
    logEll <- ellLogFromNormalizedVector(v)
  else
    logEll <- ellLogFromVector(v)
  eig <- eigen(logEll, symmetric=TRUE)
  # Everything else follows from that diagonalization.
  logA <- -0.5 * eig$values
  a <- exp(logA)
  rotation <- t(eig$vectors)
  if (det(rotation) < 0)
    rotation[3,] <- -rotation[3,]
  tensor <- eig$vectors %*% diag(exp(eig$values)) %*% t(eig$vectors)
  list(vector=v, a=a, logA=logA, rotation=rotation, tensor=tensor)
}

#' An ellipsoid from its orientation and the logarithms of its semi-axis lengths.
#' 
#' @param r A real 3x3 matrix (special orthogonal), with the semi-axis directions along its rows.
#' @param logA A real 3-dimensional vector. The logarithms of the semi-axis lengths, in order corresponding to the rows of r.
#' @return An ellipsoid.
ellEllipsoidFromRotationLogA <- function(r, logA, doNormalize=FALSE) {
  if (doNormalize)
    logA <- logA - sum(logA) / 3
  a <- exp(logA)
  tensor <- t(r) %*% diag(a^-2) %*% r
  logEll <- t(r) %*% diag(-2 * logA) %*% r
  if (doNormalize)
    v <- ellNormalizedVectorFromLog(logEll)
  else
    v <- ellVectorFromLog(logEll)
  list(vector=v, a=a, logA=logA, rotation=r, tensor=tensor)
}

ellTensorFromRotationA <- function(r, a) {
  t(r) %*% diag(a^-2) %*% r
}

ellRotationAFromTensor <- function(e) {
  eig <- eigen(e, symmetric=TRUE)
  a <- eig$values^(-1 / 2)
  r <- t(eig$vectors)
  if (det(r) < 0)
    r[3,] <- -r[3,]
  list(rotation=r, a=a)
}

ellEllipsoidFromTensor <- function(e, doNormalize=FALSE) {
  ra <- ellRotationAFromTensor(e)
  logA <- log(ra$a)
  if (doNormalize)
    logA <- logA - sum(logA) / 3
  logEll <- t(ra$rotation) %*% diag(-2 * logA) %*% ra$rotation
  if (doNormalize)
    v <- ellNormalizedVectorFromLog(logEll)
  else
    v <- ellVectorFromLog(logEll)
  list(vector=v, logA=logA, a=exp(logA), rotation=ra$rotation, tensor=e)
}

#' Ellipsoid orientation with axes in ascending or descending order.
#' 
#' @param rot A real 3x3 matrix (special orthogonal), with the semi-axis directions along its rows.
#' @param aOrLogA A real 3-dimensional vector. The semi-axis lengths or their logarithms, in order corresponding to the rows of rot.
#' @param descending A Boolean. If TRUE, then sort axes in descending order. If FALSE, then sort in ascending order.
#' @return A real 3x3 matrix (special orthogonal). This matrix equals the input matrix, with its rows reordered and possibly one row negated to maintain determinant 1.
ellAscendingRotation <- function(rot, aOrLogA, descending=FALSE) {
  ord <- order(aOrLogA, decreasing=descending)
  first <- rot[ord[[1]],]
  second <- rot[ord[[2]],]
  third <- rot[ord[[3]],]
  if (dot(cross(first, second), third) > 0)
    rbind(first, second, third)
  else
    rbind(first, second, -third)
}

ellDescendingRotation <- function(rot, aOrLogA) {
  ellAscendingRotation(rot, aOrLogA, descending=TRUE)
}

#' Ellipsoid orientation with axes ordered short, then long, then intermediate.
#' 
#' This function is useful for comparing ellipsoid orientations to foliation-lineation orientations.
#' @param rot A real 3x3 matrix (special orthogonal), with the semi-axis directions along its rows.
#' @param aOrLogA A real 3-dimensional vector. The semi-axis lengths or their logarithms, in order corresponding to the rows of rot.
#' @return A real 3x3 matrix (special orthogonal). This matrix equals the input matrix, with its rows reordered and possibly one row negated to maintain determinant 1.
ellPoleDirectionRotation <- function(rot, aOrLogA) {
  ord <- order(aOrLogA)
  first <- rot[ord[[1]],]
  second <- rot[ord[[3]],]
  third <- rot[ord[[2]],]
  if (dot(cross(first, second), third) > 0)
    rbind(first, second, third)
  else
    rbind(first, second, -third)
}



### SIZE AND SHAPE ###

#' Size-related tensor invariant of ellipsoids. Tantamount to volume.
#' 
#' The first invariant (trace) of log E^(-1 / 2), where E is the ellipsoid tensor. The volume of the ellipsoid is (4 pi / 3) * exp(size). The size is positive for large ellipsoids, zero for normalized ellipsoids, and negative for small ellipsoids.
#' @param logs A real 3D vector. The logs of the ellipsoid's semi-axis lengths, in any order.
#' @return A real number (can achieve any real value).
ellSizeInvariant <- function(logs) {
  sum(logs)
}

#' Strain-related tensor invariant of ellipsoids. Tantamount to octahedral shear strain.
#'
#' The second invariant of log E^(-1 / 2), where E is the ellipsoid tensor. Equals zero for spheres. For normalized ellipsoids, this strain == -Es^2 / 2, where Es is the octahedral shear strain. 
#' @param logs A real 3D vector. The logs of the ellipsoid's semi-axis lengths, in any order.
#' @return A real number <= 0.
ellStrainInvariant <- function(logs) {
  logs[[1]] * logs[[2]] + logs[[2]] * logs[[3]] + logs[[3]] * logs[[1]]
}

#' Shape-related tensor invariant of ellipsoids. An analogue of Lode's parameter.
#'
#' The third invariant (determinant) of log E^(-1 / 2), where E is the ellipsoid tensor. For normalized ellipsoids, this shape is positive for oblate ellipsoids, zero for 'plane strain' ellipsoids, and negative for prolate ellipsoids. In this sense it is analogous to (but not equal to or even tantamount to) Lode's parameter.
#' @param logs A real 3D vector. The logs of the ellipsoid's semi-axis lengths, in any order.
#' @return A real number (can achieve any real value).
ellShapeInvariant <- function(logs) {
  logs[[1]] * logs[[2]] * logs[[3]]
}

#' Volume of an ellipsoid.
#' 
#' @param logs A real 3D vector. The logs of the ellipsoid's semi-axis lengths, in any order.
#' @return A real number > 0.
ellVolume <- function(logs) {
  exp(sum(logs)) * 4 * pi / 3
}

#' Octahedral shear strain e_s.
#' 
#' @param logs A real 3D vector. The logs of the ellipsoid's semi-axis lengths, in any order.
#' @return A real number >= 0.
ellOctahedralShearStrain <- function(logs) {
  mostOfIt <- (logs[[1]] - logs[[2]])^2 + (logs[[2]] - logs[[3]])^2 + (logs[[3]] - logs[[1]])^2
  sqrt(mostOfIt / 3)
}

# Lode's parameter nu.
#' 
#' nu is undefined for spheres, but we arbitrarily declare nu = 0 for them. Otherwise -1 <= nu <= 1. nu = -1 for prolate spheroids and nu = 1 for oblate spheroids.
#' @param logs A real 3D vector. The logs of the ellipsoid's semi-axis lengths, in any order.
#' @return A real number (in the interval [-1, 1]), unless it fails.
ellLodeNu <- function(logs) {
  # Sort the logs so that l1 >= l2 >= l3.
  l1 <- max(logs)
  l3 <- min(logs)
  l2 <- sum(logs) - l1 - l3
  if (l1 == l3)
    0
  else
    (2 * l2 - l1 - l3) / (l1 - l3)
}

#' The statistic P_j of Jelinek (1981).
#' 
#' @param logs A real 3D vector. The logs of the ellipsoid's semi-axis lengths, in any order.
#' @return A real number >= 1.
ellJelinekP <- function(logs) {
  # We used to have an '8' in the square root, where Jelinek had a '2', because we thought that he was working with eta_i = -2 l_i. But on 2016/03/10 Mike Jackson at the IRM told me that, according to all recent authors, AMS ellipsoids are 'magnitude ellipsoids'(e.g., Hrouda, 1982), whose semi-axis lengths are the principal susceptibilities, which are the eigenvalues of the susceptibility tensor. So it seems that eta_i = l_i now. And saying so reproduces the IRM's computed P_j.
  v <- logs - sum(logs) / 3
  exp(sqrt(2 * dot(v, v)))
}

#' Flinn's K measure of ellipsoid shape.
#' 
#' Fails in the case of a prolate spheroid or sphere. Zero in the case of an oblate spheroid (that is not a sphere).
#' @param logs A real 3D vector. The logs of the ellipsoid's semi-axis lengths, in any order.
#' @return A real number >= 0, unless it fails.
ellFlinnK <- function(logs) {
  # Sort the logs so that l1 >= l2 >= l3.
  l1 <- max(logs)
  l3 <- min(logs)
  l2 <- sum(logs) - l1 - l3
  a1 <- exp(l1)
  a2 <- exp(l2)
  a3 <- exp(l3)
  (a1 / a2 - 1) / (a2 / a3 - 1)
}

#' The logs of an ellipsoid's semi-axis lengths, from three other measures of shape.
#' 
#' @param vEsNu A real 3D vector consisting of volume, octahedral shear strain, and Lode's parameter nu.
#' @return A real 3D vector, consisting of the logs of the ellipsoid's semi-axis lengths.
ellLogAFromVEsNu <- function(vEsNu) {
  # Invert the volume-logs relationship.
  sumOfLogs <- log(vEsNu[[1]] * 3 / (4 * pi))
  # logs1 = alpha + beta logs3.
  alpha <- sumOfLogs * 2 / (vEsNu[[3]] + 3)
  beta <- (vEsNu[[3]] - 3) / (vEsNu[[3]] + 3)
  # logs2 = gamma + delta logs3.
  gamma <- (1 - 2 / (vEsNu[[3]] + 3)) * sumOfLogs
  delta <- -1 - beta
  # Compute the coefficients of the quadratic.
  aa <- 2 - 2 * beta + 2 * beta^2 - 2 * delta - 2 * beta * delta + 2 * delta^2
  bb <- -2 * alpha + 4 * alpha * beta - 2 * gamma - 2 * beta * gamma - 2 * alpha * delta + 4 * gamma * delta
  cc <- 2 * alpha^2 - 2 * alpha * gamma + 2 * gamma^2 - 3 * vEsNu[[2]]^2
  # Solve the quadratic aa logs3^2 + bb logs3 + cc == 0 and back out the other logs.
  logs3 <- realQuadraticSolutions(aa, bb, cc)
  logs1 <- sapply(logs3, function(l3) {alpha + beta * l3})
  logs2 <- sapply(logs3, function(l3) {gamma + delta * l3})
  sols <- lapply(1:length(logs3), function(i) c(logs1[[i]], logs2[[i]], logs3[[i]]))
  # Choose the solution such that logs1 >= logs2 >= logs3, as required by nu.
  sols <- Filter(function(sol) {sol[[1]] >= sol[[2]] && sol[[2]] >= sol[[3]]}, sols)
  if (length(sols) != 1) {
    print("warning: ellLogsFromVEsNu: did not find one and only one solution as expected")
    print(sols)
  }
  sols[[1]]
}



### DESCRIPTIVE STATISTICS ###

#' Geometric mean.
#' 
#' @param vectors A list of 5- or 6-dimensional real vectors.
#' @return An ellipsoid.
ellMean <- function(vectors) {
  ellEllipsoidFromVector(arithmeticMean(vectors))
}

#' Covariance matrix.
#' 
#' The vectors are automatically centered about their mean. The denominator is n - 1, not n.
#' @param vectors A list of 5- or 6-dimensional real vectors.
#' @return A 5x5 or 6x6 real matrix.
ellCovariance <- function(vectors) {
  var(t(simplify2array(vectors)))
}

#' Convenience shortcut to the eigenvalues of the covariance matrix.
#' 
#' These 5 or 6 numbers quantify the dispersion of the ellipsoids.
#' @param vectors A list of 5- or 6-dimensional real vectors.
#' @return A 5- or 6-dimensional vector, containing the eigenvalues of the variance.
ellCovarianceScalars <- function(vectors) {
  eigen(ellCovariance(vectors), only.values=TRUE)$values
}

#' Principal component analysis.
#' 
#' @param vectors A list of 5- or 6-dimensional real vectors.
#' @return See prcomp.
ellPrincipalComponentAnalysis <- function(vectors) {
  prcomp(t(simplify2array(vectors)))
}



### INFERENCE ###

#' One-sample inference based on Hotelling's T2 test.
#' 
#' @param vectors A list of 5- or 6-dimensional real vectors.
#' @param hypoth A 5- or 6-dimensional real vector, respectively.
#' @param fOrChi Character. Should be 'f' or 'chi'.
#' @return See HotellingsT2.
ellHotellingT2Inference <- function(vectors, hypoth, fOrChi="f") {
  HotellingsT2(X=t(simplify2array(vectors)), mu=hypoth, test=fOrChi)
}

#' One-sample inference based on Hotelling test using the MM estimator.
#' 
#' @param vectors A list of 5- or 6-dimensional real vectors.
#' @param hypoth A 5- or 6-dimensional real vector, respectively.
#' @param numBoots A real number (positive integer). The number of bootstrap samples.
#' @return See FRBhotellingMM.
ellBootstrapMMInference <- function(vectors, hypoth, numBoots=1000, ...) {
  FRBhotellingMM(X=t(simplify2array(vectors)), mu0=hypoth, R=numBoots, ...)
}

#' One-sample inference based on Hotelling test using the S estimator.
#' 
#' @param vectors A list of 5- or 6-dimensional real vectors.
#' @param hypoth A 5- or 6-dimensional real vector, respectively.
#' @param numBoots A real number (positive integer). The number of bootstrap samples.
#' @return See FRBhotellingS.
ellBootstrapSInference <- function(vectors, hypoth, numBoots=1000, ...) {
  FRBhotellingS(X=t(simplify2array(vectors)), mu0=hypoth, R=numBoots, ...)
}

# 32 or 64 corners of a crude confidence region.
ellCICombinationVectors <- function(ci) {
  if (!class(ci) == "matrix")
    list(ci[[1]], ci[[2]])
  else {
    recursive <- ellCICombinationVectors(ci[,(2:ncol(ci))])
    c(lapply(recursive, function(rec) c(ci[[1, 1]], rec)),
      lapply(recursive, function(rec) c(ci[[2, 1]], rec)))
  }
}

# Generates approximately 7^(dim - 1) points on the (dim - 1)-dimensional unit sphere in dim-dimensional Euclidean space.
ellHighSphereVectors <- function(ambientDimension, numSamples=7) {
  if (ambientDimension == 0)
    list()
  else if (ambientDimension == 1)
    list(1, -1)
  else if (ambientDimension == 2)
    lapply(0:(numSamples - 1),
           function(i) {a <- (i * 2 * pi + 1) / numSamples; c(sin(a), cos(a))})
  else {
    recursive <- ellHighSphereVectors(ambientDimension - 1, numSamples)
    unlist(lapply(
      recursive,
      function(rec) lapply(0:(numSamples - 1),
                           function(i) {a <- (i * 2 * pi + 1) / numSamples; c(sin(a) * rec, cos(a))})),
      recursive=FALSE, use.names=FALSE)
  }
}

#' A sampling of points on an ellipsoid in high-dimensional space.
#' 
#' d is the dimension of the ambient space. The arguments for this function typically come out of ellBootstrapInference, where d == 5 or d == 6.
#' @param covarInv A real dxd matrix (symmetric, positive-definite).
#' @param center A real dD vector.
#' @param level A real number. Typically q095^2 from ellBootstrapInference.
#' @param numSamples A real number (positive integer). Roughly, the number of samples per dimension.
#' @return A list of dD real vectors.
ellHighEllipsoidVectors <- function(covarInv, center, level, numSamples=7) {
  eig <- eigen(covarInv, symmetric=TRUE)
  q <- eig$vectors
  a <- sqrt(level) * eig$values^(-0.5)
  sphere <- ellHighSphereVectors(ncol(covarInv), numSamples)
  lapply(sphere, function(v) as.numeric((q %*% (a * v)) + center))
}

# Confidence region based on percentiles of Mahalanobis distance.
ellMahalanobisInference <- function(ss, sBar) {
  vs <- lapply(ss, function(s) {s - sBar})
  covar <- arithmeticMean(lapply(vs, function(v) {outer(v, v)}))
  covarInv <- solve(covar)
  norms <- sapply(vs, function(v) {sqrt(v %*% covarInv %*% v)})
  empiricalCDF <- ecdf(norms)
  # Build the p-value function.
  f <- function(s) {
    v <- s - sBar
    1 - empiricalCDF(sqrt(v %*% covarInv %*% v))
  }
  # Compute a few popular percentiles.
  qs <- quantile(norms, probs=c(0.00, 0.25, 0.50, 0.75, 0.95, 0.99, 1.00), names=FALSE)
  list(pvalue=f, center=sBar, covarInv=covarInv,
       q000=qs[[1]], q025=qs[[2]], q050=qs[[3]], q075=qs[[4]], q095=qs[[5]], q099=qs[[6]], q100=qs[[7]])
}

#' One-sample inference based on bootstrapping.
#' 
#' @param vectors A list of 5- or 6-dimensional real vectors.
#' @param numBoots A real number (positive integer). The number of bootstrap samples.
#' @return A list with members $pvalue, $center, $covarInv, $bootstraps, $q000, $q025, $q050, $q075, $q095, $q099, $q100. $pvalue is a function with input a 5D or 6D vector and output a real number, the p-value for the null hypothesis that the mean is that vector. $center is the mean of $bootstraps, which are the bootstraps. $covarInv is their inverse covariance matrix at the mean. The $qxxx values are percentiles of Mahalanobis distance among the bootstraps.
ellBootstrapInference <- function(vectors, numBoots=1000) {
  boots <- replicate(numBoots, arithmeticMean(sample(vectors, length(vectors), replace=TRUE)), simplify=FALSE)
  bootMean <- arithmeticMean(boots)
  infer <- ellMahalanobisInference(boots, bootMean)
  infer$bootstraps <- boots
  infer
}

#' Two-sample inference based on Hotelling's T2 test.
#' 
#' @param firsts A list of 5- or 6-dimensional real vectors.
#' @param seconds A list of 5- or 6-dimensional real vectors.
#' @param fOrChi Character. Should be 'f' or 'chi'.
#' @return See HotellingsT2.
ellTwoSampleHotellingT2Inference <- function(firsts, seconds, fOrChi="f") {
  HotellingsT2(X=t(simplify2array(firsts)), Y=t(simplify2array(seconds)), test=fOrChi)
}

#' Two-sample inference based on bootstrapping.
#' 
#' @param firsts A list of 5- or 6-dimensional real vectors.
#' @param seconds A list of 5- or 6-dimensional real vectors.
#' @param numBoots A real number (positive integer). The number of bootstrap samples.
#' @return A list with members $pvalue, $center, $covarInv, $bootstraps, $q000, $q025, $q050, $q075, $q095, $q099, $q100. $pvalue is a function with input a 5D or 6D vector and output a real number, the p-value for the null hypothesis that second-mean - first-mean is that vector. $center is the mean of $bootstraps, which are the bootstraps. $covarInv is their inverse covariance matrix at the mean. The $qxxx values are percentiles of Mahalanobis distance among the bootstraps.
ellTwoSampleBootstrapInference <- function(firsts, seconds, numBoots=1000) {
  f <- function() {
    firstMean <- arithmeticMean(sample(firsts, length(firsts), replace=TRUE))
    secondMean <- arithmeticMean(sample(seconds, length(seconds), replace=TRUE))
    secondMean - firstMean
  }
  boots <- replicate(numBoots, f(), simplify=FALSE)
  bootMean <- arithmeticMean(boots)
  infer <- ellMahalanobisInference(boots, bootMean)
  infer$bootstraps <- boots
  infer
}



### FITTING ELLIPSOIDS TO ELLIPTICAL SECTIONS (ROBIN, 2002) ###

# poleRakeOther is a rotation matrix with rows pointing along pole, rake, and other direction.
# Returns coefficients of B11, B12, B13, B23, B22, 1, and extra variable C in three equations.
ellRobinCoefficients <- function(poleRakeOther, rakeSemiaxisLength, otherSemiaxisLength) {
  # l is the matrix with columns pole, rake, other. Only its second and third columns will be used.
  l <- t(poleRakeOther)
  # First equation, based on (1, 1) entry.
  first <- c(
    l[[1, 2]]^2 - l[[3, 2]]^2,
    2 * l[[1, 2]] * l[[2, 2]],
    2 * l[[1, 2]] * l[[3, 2]],
    2 * l[[2, 2]] * l[[3, 2]],
    l[[2, 2]]^2 - l[[3, 2]]^2,
    0,
    -rakeSemiaxisLength^-2)
  # Second equation, based on (2, 2) entry.
  second <- c(
    l[[1, 3]]^2 - l[[3, 3]]^2,
    2 * l[[1, 3]] * l[[2, 3]],
    2 * l[[1, 3]] * l[[3, 3]],
    2 * l[[2, 3]] * l[[3, 3]],
    l[[2, 3]]^2 - l[[3, 3]]^2,
    0,
    -otherSemiaxisLength^-2)
  # Third equation, based on (1, 2) or (2, 1) entry.
  third <- c(
    l[[1, 2]] * l[[1, 3]] - l[[3, 2]] * l[[3, 3]],
    l[[1, 3]] * l[[2, 2]] + l[[1, 2]] * l[[2, 3]],
    l[[1, 3]] * l[[3, 2]] + l[[1, 2]] * l[[3, 3]],
    l[[2, 3]] * l[[3, 2]] + l[[2, 2]] * l[[3, 3]],
    l[[2, 2]] * l[[2, 3]] - l[[3, 2]] * l[[3, 3]],
    3 * l[[3, 2]] * l[[3, 3]],
    0)
  list(first, second, third)
}

#' Fit an ellipsoid to elliptical sections, using their shape but not size.
#' 
#' Warning: This function is not well tested. Actually I have reason to believe that it is quite wrong. Anyway, this is the second case treated by Robin (2002). The output ellipsoid tensor is normalized to have trace 3, and might not actually be positive-definite at all.
#' @param poleRakeOthers A list of real 3x3 matrices (special orthogonal). Each matrix describes the orientation of an ellipse in space. The first row is the pole to the plane. The second row is the rake of one of the ellipse's axes in that plane. The third row is the cross product of the first two.
#' @param rakeSemiaxisLengths A vector of real numbers. The length of the semi-axis indicated by the rake in the first argument.
#' @param otherSemiaxisLengths A vector of real numbers. The length of the semi-axis perpendicular to the rake.
#' @return A real 3x3 matrix (symmetric, trace-3). The putative ellipsoid tensor.
ellRobin <- function(poleRakeOthers, rakeSemiaxisLengths, otherSemiaxisLengths) {
  # Construct a system X B = Y of linear equations.
  n <- length(poleRakeOthers)
  x <- matrix(0, 3 * n, 5 + n)
  y <- replicate(3 * n, 0)
  for (i in 1:n) {
    eqns <- ellRobinCoefficients(poleRakeOthers[[i]], rakeSemiaxisLengths[[i]], otherSemiaxisLengths[[i]])
    x[(i * 3 - 2),1:5] <- eqns[[1]][1:5]
    x[(i * 3 - 2),(5 + i)] <- eqns[[1]][7]
    x[(i * 3 - 1),1:5] <- eqns[[2]][1:5]
    x[(i * 3 - 1),(5 + i)] <- eqns[[2]][7]
    x[(i * 3),1:5] <- eqns[[3]][1:5]
    x[(i * 3),(5 + i)] <- eqns[[3]][7]
    y[[i * 3 - 2]] <- -eqns[[1]][[6]]
    y[[i * 3 - 1]] <- -eqns[[2]][[6]]
    y[[i * 3]] <- -eqns[[3]][[6]]
  }
  # Solve for B.
  fit <- lm.fit(x, y)
  # For now, just rebuild the trace-3 ellipsoid tensor.
  es <- fit$coefficients
  rbind(c(es[[1]], es[[2]], es[[3]]),
        c(es[[2]], es[[4]], es[[5]]),
        c(es[[3]], es[[5]], 3 - es[[1]] - es[[4]]))
}

#' Fit an ellipsoid to elliptical sections, using their shape but not size.
#' 
#' This is my custom method for fitting SPO. Unlike the method of Robin (2002, Case 2), this method is guaranteed to produce a positive-definite ellipsoid tensor E. Currently I force volume normalization (det E = 1) as well. Works well, except when minEigenvalue is negative. Works less well if BFGS is replaced with the default (Nelder-Mead).
#' @param poleRakeOthers A list of real 3x3 matrices (special orthogonal). Each matrix describes the orientation of an ellipse in space. The first row is the pole to the plane. The second row is the rake of one of the ellipse's axes in that plane. The third row is the cross product of the first two.
#' @param rakeSemiaxisLengths A vector of real numbers. The length of the semi-axis indicated by the rake in the first argument.
#' @param otherSemiaxisLengths A vector of real numbers. The length of the semi-axis perpendicular to the rake.
#' @param numSteps A real number (positive integer). The number of steps to use in the optimization algorithm.
#' @return A list with members $ellipsoid, $error, $minEigenvalue, and $value. $ellipsoid is an ellipsoid. $error is an error code; if it is non-zero, then an error occurred; try increasing numSteps. $minEigenvalue is the least eigenvalue of the Hessian at the putative optimum; if it is non-positive, then an error occurred. $value is the value of the misfit function at the optimum.
ellSPO <- function(poleRakeOthers, rakeSemiaxisLengths, otherSemiaxisLengths, numSteps=10000) {
  # Pre-process the data.
  n <- length(poleRakeOthers)
  ls <- lapply(poleRakeOthers, function(r) r[2:3,])
  bs <- thread(function(a1, a2) diag(c(a1, a2)^-2), rakeSemiaxisLengths, otherSemiaxisLengths)
  # Build the misfit function to be minimized.
  misfit <- function(pars) {
    s <- rbind(
      c(pars[[1]], pars[[2]], pars[[3]]), 
      c(pars[[2]], pars[[4]], pars[[5]]), 
      c(pars[[3]], pars[[5]], -pars[[1]] - pars[[4]]))
    e <- ellExp(s)
    diffs <- lapply(1:n, function(i) {exp(pars[[5 + i]]) * bs[[i]] - ls[[i]] %*% e %*% t(ls[[i]])})
    normSqs <- sapply(diffs, function(diff) tr(t(diff) %*% diff))
    sum(normSqs)
  }
  # Seed the minimization from the unit sphere and all k_i = 0.5 arbitrarily.
  seed <- c(0, 0, 0, 0, 0, replicate(n, 0.5))
  solution <- optim(seed, misfit, hessian=TRUE, method="BFGS", control=list(maxit=numSteps))
  # Report the answer and diagnostic information.
  eigvals <- eigen(solution$hessian, symmetric=TRUE, only.values=TRUE)$values
  s <- rbind(
    c(solution$par[[1]], solution$par[[2]], solution$par[[3]]), 
    c(solution$par[[2]], solution$par[[4]], solution$par[[5]]), 
    c(solution$par[[3]], solution$par[[5]], -solution$par[[1]] - solution$par[[4]]))
  ell <- ellEllipsoidFromTensor(ellExp(s), doNormalize=TRUE)
  list(ellipsoid=ell, error=solution$convergence, minEigenvalue=min(eigvals), value=solution$value)
}

# Testing for ellSPO. Noiseless. Test exactly mirrors the optimization procedure.
ellSPOTest <- function(n) {
  # Make a random ellipsoid.
  q <- rotUniform()
  a <- exp(rnorm(3))
  e <- t(q) %*% diag(a^-2) %*% q
  # Make n random sections.
  f <- function() {
    l <- rotUniform()[2:3,]
    b <- l %*% e %*% t(l)
    eig <- eigen(b, symmetric=TRUE)
    l <- t(eig$vectors) %*% l
    b <- l %*% e %*% t(l)
    rake <- b[[1, 1]]^(-1 / 2)
    other <- b[[2, 2]]^(-1 / 2)
    list(l=l, rake=rake, other=other)
  }
  sections <- replicate(n, f(), simplify=FALSE)
  # Dissect the sections into the format desired by ellRobin and ellSPO.
  poleRakeOthers <- lapply(sections, 
                           function(s) rbind(cross(s$l[1,], s$l[2,]), s$l[1,], s$l[2,]))
  rakeSemiaxisLengths <- sapply(sections, function(s) s$rake)
  otherSemiaxisLengths <- sapply(sections, function(s) s$other)
  # Compare the true answer to the deduced answer.
  print("true ellipsoid tensor, volume-normalized:")
  print(e * det(e)^(-1 / 3))
  print("ellSPO result:")
  pred <- ellSPO(poleRakeOthers, rakeSemiaxisLengths, otherSemiaxisLengths)
  print(pred$ellipsoid$tensor)
  print(c(pred$error, pred$minEigenvalue))
  print("ellRobin result, volume-normalized:")
  pred <- ellRobin(poleRakeOthers, rakeSemiaxisLengths, otherSemiaxisLengths)
  print(pred * det(pred)^(-1 / 3))
}

# This is an example hand-ported from Mathematica. Again noiseless, but the sections are not being generated by the same code that does the optimization. We should get rbind(c(0.879642, -0.0768036, -0.0419311), c(-0.0768036, 1.06686, 0.0109123), c(-0.0419311, 0.0109123, 1.07437)).
ellSPOTestSpecific <- function() {
  rakeOthers <- list(
    rbind(
      c(-0.336673, 0.941231, 0.0271225),
      c(-0.941158, -0.335463, -0.0410727)),
    rbind(
      c(-0.251698, 0.938829, 0.23505),
      c(-0.506325, 0.0792426, -0.858694)),
    rbind(
      c(-0.263783, 0.947118, -0.182718),
      c(0.865012, 0.31609, 0.389668)),
    rbind(
      c(0.065659, -0.526426, 0.847682),
      c(-0.324202, -0.814681, -0.48082)))
  rakeSemiaxisLengths <- c(0.955355, 0.952817, 0.960189, 0.970211)
  otherSemiaxisLengths <- c(1.08491, 1.00371, 1.07812, 0.998093)
  poleRakeOthers <- lapply(
    rakeOthers, function(ro) rotProjectedMatrix(rbind(cross(ro[1,], ro[2,]), ro[1,], ro[2,])))
  ellSPO(poleRakeOthers, rakeSemiaxisLengths, otherSemiaxisLengths)
}



### PLOTTING ###

#' Visualization of ellipsoids.
#' 
#' @param rots A list of 3x3 real matrices (special orthogonal). The ellipsoid orientations, as given by the $rotation field of an ellipsoid.
#' @param as A list of 3D real vectors. The ellipsoid semi-axis lengths, in order corresponding to the rows of the rots, as given by the $a field of an ellipsoid.
#' @param centers A list of 3D real vectors. The locations of the ellipsoid centers.
#' @param numNonAdapt A real number (non-negative integer). The number of refinements to use. Each refinement makes the ellipsoids smoother, but increases time and memory requirements by a factor of four.
#' @return NULL.
ellEllipsoidPlot <- function(rots, as, centers=replicate(length(rots), c(0, 0, 0), simplify=FALSE), numNonAdapt=4, ...) {
  sphere <- rayTetrahedralSphere(numNonAdapt)
  triangles <- unlist(
    lapply(1:length(as),
           function(i) lapply(sphere, function(tri) 
             lapply(tri, function(v) {centers[[i]] + as.numeric(t(rots[[i]]) %*% (as[[i]] * v))}))),
    recursive=FALSE, use.names=FALSE)
  plot3D(triangles=triangles, ...)
}

#' Equal-area plot of ellipsoid axes.
#' 
#' Short axes are shown as circles, intermediate as triangles, long as squares. Warning: Curves are not well tested.
#' @param rots A list of 3x3 real matrices (special orthogonal). The ellipsoid orientations, as given by the $rotation field of an ellipsoid.
#' @param as A list of 3D real vectors. The ellipsoid semi-axis lengths, in order corresponding to the rows of the rots, as given by the $a field of an ellipsoid. Alternatively, you can pass logA; this information is used only to determine the order of the axes.
#' @param rotCurves A list of lists of 3x3 real matrices (special orthogonal). Like rots, but curves rather than points.
#' @param aCurves A list of lists of 3D real vectors. Like as, but curves rather than points.
#' @param colors Character. A vector of colors for coloring the points, as in all R graphics functions.
#' @return NULL.
ellEqualAreaPlot <- function(rots, as, rotCurves=list(), aCurves=list(), colors=c("black")) {
  # Prepare to plot points based on rots and as.
  f <- function(i, rs, as) {
    ord <- order(as[[i]])
    list(rs[[i]][ord[[1]],], rs[[i]][ord[[2]],], rs[[i]][ord[[3]],])
  }
  points <- unlist(lapply(1:length(rots), f, rots, as), recursive=FALSE, use.names=FALSE)
  # Prepare to plot curves based on curvesRots and curvesAs.
  if (length(rotCurves) >= 1 && length(rotCurves) == length(aCurves)) {
    curves <- lapply(1:length(rotCurves), function(j) lapply(1:length(rotCurves[[j]]), f, rotCurves[[j]], aCurves[[j]]))
    curves1 <- lapply(curves, function(curve) lapply(curve, function(tri) tri[[1]]))
    curves2 <- lapply(curves, function(curve) lapply(curve, function(tri) tri[[2]]))
    curves3 <- lapply(curves, function(curve) lapply(curve, function(tri) tri[[3]]))
    curves <- c(curves1, curves2, curves3)
  } else
    curves <- list()
  # Plot.
  newColors <- as.character(sapply(colors, function(s) c(s, s, s)))
  lineEqualAreaPlot(points, curves=curves, colors=newColors, shapes=c("c", "t", "s"))
}

#' Equal-volume plot of ellipsoid orientations.
#' 
#' @param rots A list of 3x3 real matrices (special orthogonal). The ellipsoid orientations, as given by the $rotation field of an ellipsoid.
#' @param as A list of 3D real vectors. The ellipsoid semi-axis lengths, in order corresponding to the rows of the rots, as given by the $a field of an ellipsoid. Alternatively, you can pass logA; this information is used only to determine the order of the axes.
#' @param rotCurves A list of lists of 3x3 real matrices (special orthogonal). Like rots, but curves rather than points.
#' @param aCurves A list of lists of 3D real vectors. Like as, but curves rather than points.
#' @param colors Character. A vector of colors for coloring the points, as in all R graphics functions.
#' @return NULL.
ellEqualVolumePlot <- function(rots, as, rotCurves=list(), aCurves=list(), colors=c("white"), ...) {
  # The rotations are permuted into this row-order: short, long, intermediate. To match other plane-line stuff in our library.
  f <- function(i, rs, as) {
    ord <- order(as[[i]])
    short <- rs[[i]][ord[[1]],]
    long <- rs[[i]][ord[[3]],]
    rbind(short, long, cross(short, long))
  }
  points <- lapply(1:length(rots), f, rots, as)
  if (length(rotCurves) >= 1 && length(rotCurves) == length(aCurves)) {
    curves <- lapply(1:length(rotCurves), function(j) lapply(1:length(rotCurves[[j]]), f, rotCurves[[j]], aCurves[[j]]))
    f <- function(curve) {
      cur <- list(curve[[1]])
      for (r in curve[2:length(curve)])
        cur[[length(cur) + 1]] <- oriNearestRepresentative(r, cur[[length(cur)]], oriLineInPlaneGroup)
      cur
    }
    curves <- lapply(curves, f)
  } else
    curves <- list()
  oriEqualVolumePlot(points=points, curves=curves, group=oriLineInPlaneGroup, colors=colors, ...)
}

#' A distorted version of the Hsu-Nadai plot of ellipsoid shapes.
#' 
#' This is a polar plot, in which the radial coordinate is octahedral shear strain and the angular coordinate is Lode's parameter. This plot is similar to, but not identical to, the Hsu-Nadai plot. See ellHsuNadaiPlot for the real thing.
#' @param logAs A list of 3D real vectors. The ellipsoid semi-axis log-lengths.
#' @param curves A list of lists of 3D real vectors. Like logAs, but curves rather than points.
#' @param es A real number (positive) or NULL. If a number, then that is the radius of the plot in the E_s direction. If NULL, then the radius of the plot is inferred from the points (not the curves, currently).
#' @param colors Character. A vector of colors for coloring the points, as in all R graphics functions.
#' @return NULL.
ellWrongHsuNadaiPlot <- function(logAs, curves=list(), es=NULL, colors=c("black")) {
  ess <- sapply(logAs, ellOctahedralShearStrain)
  nus <- sapply(logAs, ellLodeNu)
  esCurves <- lapply(curves, function(curve) sapply(curve, ellOctahedralShearStrain))
  nuCurves <- lapply(curves, function(curve) sapply(curve, ellLodeNu))
  # Make the plot window.
  if (is.null(es))
    es <- max(c(ess, 1))
  plot.new()
  plot.window(xlim=c(-0.55 * es, 0.55 * es), ylim=c(-0.05 * es, 1.05 * es))
  # Plot the points.
  if (length(logAs) >= 1) {
    xs <- sapply(1:length(logAs), function(i) ess[[i]] * cos(pi / 2 - nus[[i]] * pi / 6))
    ys <- sapply(1:length(logAs), function(i) ess[[i]] * sin(pi / 2 - nus[[i]] * pi / 6))
    points(xs, ys, col=colors, pch=c(19))
  }
  # Plot the curves.
  if (length(curves) >= 1)
    for (j in 1:length(curves)) {
      xs <- sapply(1:length(curves[[j]]), function(i) esCurves[[j]][[i]] * cos(pi / 2 - nuCurves[[j]][[i]] * pi / 6))
      ys <- sapply(1:length(curves[[j]]), function(i) esCurves[[j]][[i]] * sin(pi / 2 - nuCurves[[j]][[i]] * pi / 6))
      lines(xs, ys)
    }
  # Plot the boundary.
  xys <- sapply(0:30, function(t) {
    theta <- (t / 30) * (pi / 3) + (pi / 3)
    es * c(cos(theta), sin(theta))
  })
  lines(xys[1,], xys[2,])
  lines(c(-0.5 * es, 0, 0.5 * es), c(sqrt(3) / 2 * es, 0, sqrt(3) / 2 * es))
  # Plot some tick marks.
  if (es >= 1)
    for (i in 1:es) {
      xs <- c(i * 0.5, i * 0.5 + sqrt(3) * 0.5 * 0.05)
      ys <- c(i * sqrt(3) * 0.5, i * sqrt(3) * 0.5 - 0.5 * 0.05)
      lines(xs, ys)
      lines(-xs, ys)
    }
}

#' Hsu-Nadai plot of ellipsoid shapes.
#' 
#' @param logAs A list of 3D real vectors. The ellipsoid semi-axis log-lengths.
#' @param curves A list of lists of 3D real vectors. Like logAs, but curves rather than points.
#' @param es A real number (positive) or NULL. If a number, then that is the radius of the plot in the E_s direction. If NULL, then the radius of the plot is inferred from the points (not the curves, currently).
#' @param colors Character. A vector of colors for coloring the points, as in all R graphics functions.
#' @return NULL.
ellHsuNadaiPlot <- function(logAs, curves=list(), es=NULL, colors=c("black")) {
  x <- function(logA) {-sqrt(1.5) * max(logA - sum(logA) / 3) - sqrt(1.5) * min(logA - sum(logA) / 3)}
  y <- function(logA) {sqrt(0.5) * max(logA - sum(logA) / 3) - sqrt(0.5) * min(logA - sum(logA) / 3)}
  # Make the plot window.
  if (is.null(es))
      es <- max(c(1, sapply(logAs, ellOctahedralShearStrain)))
  plot.new()
  plot.window(xlim=c(-0.55 * es, 0.55 * es), ylim=c(-0.05 * es, 1.05 * es))
  # Plot the points.
  if (length(logAs) >= 1)
    points(sapply(logAs, x), sapply(logAs, y), col=colors, pch=c(19))
  # Plot the curves.
  if (length(curves) >= 1)
    for (j in 1:length(curves))
      lines(sapply(curves[[j]], x), sapply(curves[[j]], y))
  # Plot the boundary.
  xys <- sapply(0:30, function(t) {
    theta <- (t / 30) * (pi / 3) + (pi / 3)
    es * c(cos(theta), sin(theta))
  })
  lines(xys[1,], xys[2,])
  lines(c(-0.5 * es, 0, 0.5 * es), c(sqrt(3) / 2 * es, 0, sqrt(3) / 2 * es))
  # Plot some tick marks.
  if (es >= 1)
    for (i in 1:es) {
      xs <- c(i * 0.5, i * 0.5 + sqrt(3) * 0.5 * 0.05)
      ys <- c(i * sqrt(3) * 0.5, i * sqrt(3) * 0.5 - 0.5 * 0.05)
      lines(xs, ys)
      lines(-xs, ys)
    }
}

#' Hsu-Nadai plot of ellipsoid shapes, with a third dimension specified by the user.
#' 
#' @param logAs A list of 3D real vectors. The ellipsoid semi-axis log-lengths.
#' @param zs A vector of real numbers. The coordinates of the points in the direction perpendicular to the Hsu-Nadai plot.
#' @param es A real number (positive) or NULL. If a number, then that is the radius of the plot in the E_s direction. If NULL, then the radius of the plot is inferred from the points.
#' @param colors Character. A vector of colors for coloring the points, as in all R graphics functions.
#' @param ... Other arguments to pass to the underlying plot3D.
#' @return NULL.
ellHsuNadaiScalarPlot <- function(logAs, zs, es=NULL, colors=c("white"), ...) {
  x <- function(logA) {-sqrt(1.5) * max(logA - sum(logA) / 3) - sqrt(1.5) * min(logA - sum(logA) / 3)}
  y <- function(logA) {sqrt(0.5) * max(logA - sum(logA) / 3) - sqrt(0.5) * min(logA - sum(logA) / 3)}
  # Determine how big the plot will be.
  if (is.null(es))
    es <- max(c(1, sapply(logAs, ellOctahedralShearStrain)))
  z <- max(abs(zs))
  radius <- max(es, z)
  # Build the points.
  points <- lapply(1:length(logAs), function(i) c(x(logAs[[i]]), y(logAs[[i]]), zs[[i]]))
  # Build the curves.
  f <- function(t, zz) {
    theta <- (t / 30) * (pi / 3) + (pi / 3)
    c(es * c(cos(theta), sin(theta)), zz)
  }
  bottom <- lapply(0:30, f, -z)
  bottom <- c(bottom, list(c(-0.5 * es, sqrt(3) / 2 * es, -z), c(0, 0, -z), c(0.5 * es, sqrt(3) / 2 * es, -z)))
  middle <- lapply(0:30, f, 0)
  middle <- c(middle, list(c(-0.5 * es, sqrt(3) / 2 * es, 0), c(0, 0, 0), c(0.5 * es, sqrt(3) / 2 * es, 0)))
  top <- lapply(0:30, f, z)
  top <- c(top, list(c(-0.5 * es, sqrt(3) / 2 * es, z), c(0, 0, z), c(0.5 * es, sqrt(3) / 2 * es, z)))
  plot3D(radius=radius, points=points, curves=list(bottom, middle, top), colors=colors, ...)
}

#' Flinn plot of ellipsoid shapes.
#' 
#' @param as A list of 3D real vectors, with all entries positive. The ellipsoid semi-axis lengths.
#' @param colors Character. A vector of colors for coloring the points, as in all R graphics functions.
#' @return NULL.
ellFlinnPlot <- function(as, colors=c("black")) {
  xs <- sapply(as, function(a) {(sum(a) - min(a) - max(a)) / min(a)})
  ys <- sapply(as, function(a) {max(a) / (sum(a) - min(a) - max(a))})
  plot(x=xs, y=ys, xlim=c(1, max(xs)), ylim=c(1, max(ys)),
       xlab="intermediate / short", ylab="long / intermediate")
}

#' Logarithmic Flinn plot (Ramsay plot) of ellipsoid shapes.
#' 
#' @param logAs A list of 3D real vectors. The ellipsoid semi-axis log-lengths.
#' @param colors Character. A vector of colors for coloring the points, as in all R graphics functions.
#' @return NULL.
ellLogFlinnPlot <- function(logAs, colors=c("black")) {
  xs <- sapply(logAs, function(logA) {-max(logA - sum(logA) / 3) - 2 * min(logA - sum(logA) / 3)})
  ys <- sapply(logAs, function(logA) {2 * max(logA - sum(logA) / 3) + min(logA - sum(logA) / 3)})
  plot(x=xs, y=ys, xlim=c(0, max(xs)), ylim=c(0, max(ys)),
       xlab="log(intermediate / short)", ylab="log(long / intermediate)")
}

#' Jelinek plot of ellipsoid shapes.
#' 
#' This is a rectangular plot of Jelinek's Pj vs. Lode's nu. Tauxe (2010) called it the Jelinek plot, after Jelinek (1981).
#' @param logAs A list of 3D real vectors. The ellipsoid semi-axis log-lengths.
#' @param colors Character. A vector of colors for coloring the points, as in all R graphics functions.
#' @return NULL.
ellJelinekPlot <- function(logAs, colors=c("black")) {
  plot(x=sapply(logAs, ellJelinekP), 
       y=sapply(logAs, ellLodeNu), 
       col=colors, xlab="P_j", ylab="nu", ylim=c(-1, 1))
}

#' Pair plot of ellipsoid vectors.
#' 
#' @param points A list of 5D or 6D real vectors. Ellipsoid vectors, as in the $vector field of an ellipsoid.
#' @param colors Character. A vector of colors for coloring the points, as in all R graphics functions.
#' @param ... Other parameters to be passed to the underlying pairs function.
#' @param NULL.
ellPairsPlot <- function(points, colors=c("black"), ...) {
  pairs(t(simplify2array(points)), labels=c("v_1", "v_2", "v_3", "v_4", "v_5"), col=colors, ...)
}

#' 2D or 3D plot of ellipsoid vectors.
#' 
#' The 2D case is like a single panel of ellPairsPlot. Warning: In 2D, curves are not well tested.
#' @param ijk A 2D or 3D vector of real numbers (positive integers). These should be in 1, ..., d, where d is the dimension of the vectors. They select out which coordinates of the vectors to display. For example, c(1, 2, 3) indicates to make a 3D plot of the first three vector coordinates.
#' @param points A list of 5D or 6D real vectors. Ellipsoid vectors, as in the $vector field of an ellipsoid.
#' @param curves A list of lists of 5D or 6D real vectors. Like points, but curves.
#' @param colors Character or NULL. A vector of colors for coloring the points, as in all R graphics functions. If NULL, then defaults to black in 2D or white in 3D.
#' @param ... Other parameters to be passed to the underlying plot3D function. Ignored in the 2D case.
#' @param NULL.
ellVectorPlot <- function(ijk, points=list(), curves=list(), colors=NULL, ...) {
  pointsNew <- lapply(points, function(v) v[ijk])
  curvesNew <- lapply(curves, function(curve) lapply(curve, function(v) v[ijk]))
  if (length(ijk) == 3) {
    if (is.null(colors))
      colors="white"
    plot3D(points=pointsNew, curves=curvesNew, colors=colors, ...)
  } else {
    if (is.null(colors))
      colors="black"
    plot(t(simplify2array(pointsNew)), col=colors,
         xlab=paste0("ellipsoid v_", as.character(ijk[[1]])),
         ylab=paste0("ellipsoid v_", as.character(ijk[[2]])))
    for (curve in curvesNew)
      lines(t(simplify2array(curve)))
  }
}


