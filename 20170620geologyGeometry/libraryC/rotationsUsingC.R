


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# In this file, 'rotation matrix' means 'special orthogonal 3x3 real matrix'. These functions complement and replace some functions in rotations.R, with versions written in C and called from R. These functions require compilation of the shared library rotationsForR.c. So you have to load it, something like this:
#   dyn.load("rotationsForR.so")



### INFERENCE ###

#' MCMC of posterior distribution rho(S, eta | D) for wrapped trivariate normal distribution parameters.
#'
#' Implemented in C for speed. See Qiu et al. (2014). Requires compilation of the C shared library rotationsForR.c.
#' @param rs A list of rotation matrices. The data set D.
#' @param numTerms A real number (non-negative integer). Controls how many terms are used in the asymptotic expansions in the Jeffreys prior for kappa.
#' @param numTuning A real number (non-negative integer). The tuning parameters are re-tuned every numTuning MCMC iterations, based on the acceptance rate since the last tuning.
#' @param numBurnin A real number (non-negative integer). A bound on the number of tuning cycles in the burn-in phase. If the system is well-tuned before numBurnin cycles, then burn-in is stopped early. If the system is not well-tuned after numBurnin cycles, then the meta-data reveal it.
#' @param numCollection A real number (non-negative integer). The number of tuning cycles in the collection phase. numTuning * numCollection should not exceed 100,000,000, or else we might overflow 32-bit integers along the way.
#' @param numReport A real number (non-negative integer). The number of MCMC samples (M, kappa) to report. If 0, then none are reported. If numTuning * numCollection (or greater), then all are reported.
#' @return When burn-in fails, a list of metadata after burn-in: nu, nuRate, gamma, gammaRate. Otherwise, a list with elements pvalue (R function from rotations to {NA} union [0, 0.1]), ms (the reported Ms), kappas (the reported kappas), mBar (rotation matrix, the mean of the collected Ms), leftCovarInv (the inverse covariance matrix of the collected Ms in the tangent space at mBar), q090, q095, q099, q100 (percentiles), metadata after burnin, and metadata after collection.
rotcWrappedTrivariateNormalMCMCInference <- function(rs, numTerms=10, numTuning=10000, numBurnin=100, numCollection=1000, numReport=10000) {
  # Check that the inputs are of the correct types.
  # !!nums are integer; rs are real
  if (numTuning * numCollection < numReport) {
    numReport <- numTuning * numCollection
    print("warning: rotcWrappedTrivariateNormalMCMCInference: clamping numReport to numTuning * numCollection")
  }
  # Flatten rs into an array of n * 9 numbers, by column-major order.
  flat <- c(simplify2array(rs))
  # Get a huge vector of numbers from C, to be parsed in a moment.
  raw <- .Call("mcmcRotationWrappedTrivariateNormalC",
               flat, numTerms, numTuning, numBurnin, numCollection, numReport)
  if (length(raw) == 4) {
    # Burn-in failed. Just report the burn-in metadata.
    list(nu=raw[[1]], nuRate=raw[[2]], gamma=exp(raw[[3]]), gammaRate=raw[[4]])
  } else {
    # raw consists of means-mean, means-covarInv, 101 percentiles, 8 meta-data, numReport means, and numReport etas.
    mBar <- matrix(raw[1:9], 3, 3)
    covarInv <- matrix(raw[10:18], 3, 3)
    percs <- raw[19:119]
    burnin <- list(nu=raw[[120]], nuRate=raw[[121]], gamma=exp(raw[[122]]), gammaRate=raw[[123]])
    collection <- list(nu=raw[[124]], nuRate=raw[[125]], gamma=exp(raw[[126]]), gammaRate=raw[[127]])
    ms <- list()
    if (numReport <= 0) 
      kappas <- c()
    else {
      for (j in 0:(numReport - 1))
        ms[[length(ms) + 1]] <- matrix(raw[(128 + j * 9):(128 + j * 9 + 8)], 3, 3)
      kappas <- exp(-raw[(128 + numReport * 9):(128 + numReport * 9 + numReport - 1)])
    }
    # p-value function based on Mahalanobis distance.
    f <- approxfun(x=percs, y=((0:100) / 1000 + 0.9), yleft=NA, yright=1)
    pvalue <- function(r) {
      v <- rotLeftTangentFromMatrix(r, mBar)
      1 - f(sqrt(v %*% covarInv %*% v))
    }
    list(pvalue=pvalue, ms=ms, kappas=kappas, mBar=mBar, leftCovarInv=covarInv,
      q090=percs[[1]], q095=percs[[51]], q099=percs[[91]], q100=percs[[101]],
      burnin=burnin, collection=collection)
  }
}

#' Bootstrapping of the sample mean.
#'
#' Similar to rotBootstrapInference with rotFrechetMean, but implemented in C for speed, and offers different percentiles of Mahalanobis norm. Requires compilation of the C shared library rotationsForR.c.
#' @param rs A list of rotation matrices.
#' @param numBoots A real number (positive integer). The number of bootstrap samples. Affects the memory requirements of the function.
#' @return A list with elements pvalue (R function from rotations to {NA} union [0, 0.1]), bootstraps (the bootstrapped means), rBar (rotation matrix, the mean of the bootstrapped means), leftCovarInv (the inverse covariance matrix at rBar), q090, q095, q099, q100 (percentiles of Mahalanobis norm).
rotcBootstrapInference <- function(rs, numBoots=10000) {
  # Check that the inputs are of the correct types.
  # !!numBoots is integer; rs are real
  # Flatten rs into an array of n * 9 numbers, by column-major order.
  flat <- c(simplify2array(rs))
  # Get a huge vector of numbers from C, to be parsed in a moment.
  raw <- .Call("pvalueRotationBootstrappingC", flat, numBoots)
  # raw consists of means-mean, means-covarInv, 101 percentiles, numBoots means.
  mBar <- matrix(raw[1:9], 3, 3)
  covarInv <- matrix(raw[10:18], 3, 3)
  percs <- raw[19:119]
  ms <- list()
  for (j in 0:(numBoots - 1))
    ms[[length(ms) + 1]] <- matrix(raw[(120 + j * 9):(120 + j * 9 + 8)], 3, 3)
  # p-value function based on Mahalanobis distance.
  f <- approxfun(x=percs, y=((0:100) / 1000 + 0.9), yleft=NA, yright=1)
  pvalue <- function(r) {
    v <- rotLeftTangentFromMatrix(r, mBar)
    1 - f(sqrt(v %*% covarInv %*% v))
  }
  list(pvalue=pvalue, bootstraps=ms, center=mBar, leftCovarInv=covarInv,
    q090=percs[[1]], q095=percs[[51]], q099=percs[[91]], q100=percs[[101]])
}



### PLOTTING ###

#' Kamb level surface for density, analogous to Kamb contouring.
#'
#' You probably do not want to call this function directly. You probably want to use rotcKambPlot. This function returns its triangles in "raw" format (a list of x-coordinates, then y-coordinates, then z-coordinates), suitable for passing to rotBallPlot as trianglesRaw. Requires compilation of the C shared library rotationsForR.c.
#' @param rs A list of rotation matrices.
#' @param multiple A real number (positive). Indicates which multiple of the standard deviation to plot; for example, mult is 12 in a 12-sigma plot.
#' @param k A real number (positive). A smoothing factor, which equaled 3 in the original paper of Kamb (1959).
#' @param degree A real number (0, 1, or 3). The degree of the weighting polynomial; higher degrees generally produce smoother plots.
#' @param numNonAdapt A real number (non-negative integer). The number of non-adaptive refinements. Time and space required are proportional to 8^(numNonAdapt + numAdapt), so don't make it too big.
#' @param numAdapt A real number (non-negative integer). The number of adaptive refinements after the non-adaptive ones. Time and space required are proportional to 8^(numNonAdapt + numAdapt), so don't make it too big.
#' @param groupSize A real number (positive integer). The size |G| of the symmetry group G acting on the rotations. If the rotations are truly rotations, without symmetry, then |G| == 1. If |G| > 1, then rs has been pre-symmetrized to hold m * |G| rotations representing m distinct orientations.
#' @return A set of raw triangles: x-coordinates, then y- and z-coordinates.
rotcKambTriangles <- function(rs, multiple, k, degree, numNonAdapt, numAdapt, groupSize=1) {
  # Check that the inputs are of the correct types.
  # !!degree and depths are integer; others are real
  # Flatten rs into an array of n * 9 numbers, by column-major order.
  flat <- c(simplify2array(rs))
  # Get a list of m * 9 numbers from C, representing 3 * m xs, 3 * m ys, and 3 * m zs.
  .Call("volumetricKambTrianglesRawC", flat, multiple, k, degree, numNonAdapt, numAdapt, groupSize)
}

#' Equal-volume plot of rotations with an accompanying Kamb density level surface.
#' 
#' This function requires compilation of the C shared library rotationsForR.c.
#' @param points A list of rotation matrices.
#' @param multiple A real number (positive). Indicates which multiple of the standard deviation to plot; for example, mult is 12 in a 12-sigma plot.
#' @param k A real number (positive). A smoothing factor, which equaled 3 in the original paper of Kamb (1959).
#' @param degree A real number (0, 1, or 3). The degree of the weighting polynomial; higher degrees generally produce smoother plots.
#' @param numNonAdapt A real number (non-negative integer). The number of non-adaptive refinements. Time and space required are proportional to 8^(numNonAdapt + numAdapt), so don't make it too big.
#' @param numAdapt A real number (non-negative integer). The number of adaptive refinements after the non-adaptive ones. Time and space required are proportional to 8^(numNonAdapt + numAdapt), so don't make it too big.
#' @param colors A list of strings (colors). Used to color the points only.
#' @param ... Plotting options: boundaryAlpha, simplePoints, etc. See rotBallPlot for details. Options about curves are ignored.
#' @return NULL.
rotcKambPlot <- function(points, multiple=6, k=3, degree=3, numNonAdapt=3, numAdapt=3, colors=c("white"), ...) {
  raws <- rotcKambTriangles(points, multiple, k, degree, numNonAdapt, numAdapt, 1)
  vs <- lapply(points, rotEqualVolumeFromMatrix)
  rotNativeEqualVolumePlot(points=vs, colors=colors, trianglesRaw=raws, ...)
}

#' Raw triangles for an ellipsoid in equal-volume coordinates.
#'
#' You probably don't want to call this function directly. You probably want to call rotcEllipsoidPlot. Returns its triangles in 'raw' format (a list of x-coordinates, then y-coordinates, then z-coordinates), suitable for passing to rotBallPlot as trianglesRaw. Requires compilation of the shared library rotationsForR.c.
#' @param center A rotation matrix.
#' @param leftCovarInv A 3x3 real matrix (symmetric).
#' @param level A real number. The ellipsoid consists of rotations whose tangent vectors v at center satisfy v %*% covarInv %*% v == level.
#' @param numNonAdapt A real number (non-negative integer). The number of non-adaptive refinements. Time and space required are proportional to 8^(numNonAdapt + numAdapt), so don't make it too big.
#' @param numAdapt A real number (non-negative integer). The number of adaptive refinements after the non-adaptive ones. Time and space required are proportional to 8^(numNonAdapt + numAdapt), so don't make it too big.
#' @return A set of raw triangles: x-coordinates, then y- and z-coordinates.
rotcEllipsoidTriangles <- function(center, leftCovarInv, level, numNonAdapt, numAdapt) {
  # Check that the inputs are of the correct types.
  # !!depths are integer; others are real
  # Get a list of m * 9 numbers from C, representing 3 * m xs, 3 * m ys, and 3 * m zs.
  .Call("volumetricEllipsoidTrianglesRawC", center, leftCovarInv, level, numNonAdapt, numAdapt)
}

#' Equal-volume plot of an ellipsoid.
#' 
#' Plots an ellipsoid based on a covariance in the left-invariant tangent space at the center. For example, to plot a 95% region produced by rotMahalanobisInference, pass leftCovarInv and q095^2. Triangles are not clipped to the edges of the plot. This function requires compilation of the C shared library rotationsForR.c.
#' @param points A list of rotation matrices. Points to plot, possibly unrelated to the ellipsoid.
#' @param center A rotation matrix.
#' @param leftCovarInv A 3x3 real matrix (symmetric).
#' @param level A real number. The ellipsoid consists of rotations whose left-invariant tangent vectors v at center satisfy v %*% leftCovarInv %*% v == level.
#' @param numNonAdapt A real number (non-negative integer). The number of non-adaptive refinements. Time and space required are proportional to 8^(numNonAdapt + numAdapt), so don't make it too big.
#' @param numAdapt A real number (non-negative integer). The number of adaptive refinements after the non-adaptive ones. Time and space required are proportional to 8^(numNonAdapt + numAdapt), so don't make it too big.
#' @param colors A list of strings (colors). Used to color the points only.
#' @param ... Plotting options: boundaryAlpha, simplePoints, etc. See rotBallPlot for details.
#' @return NULL.
rotcEllipsoidPlot <- function(points, center, leftCovarInv, level, numNonAdapt=3, numAdapt=3, colors=c("white"), ...) {
  vs <- lapply(points, rotEqualVolumeFromMatrix)
  raws <- rotcEllipsoidTriangles(center, leftCovarInv, level, numNonAdapt, numAdapt)
  if (length(raws) == 0) {
    print("warning: rotcEllipsoidPlot: No triangles generated. Try increasing numNonAdapt.")
    rotNativeEqualVolumePlot(points=vs, colors=colors, ...)
  } else
    rotNativeEqualVolumePlot(points=vs, colors=colors, trianglesRaw=raws, ...)
}


