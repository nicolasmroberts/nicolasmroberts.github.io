


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# In this file, 'rotation matrix' means 'special orthogonal 3x3 real matrix'. These functions complement and replace some functions in orientations.R, with versions written in C and called from R. These functions require compilation of the shared library orientationsForR.c. After that, you have to issue a command like this:
#   dyn.load("orientationsForR.so")



### INFERENCE ###

# Helper function.
oricMCMCParsing <- function(raw, group, numReport) {
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
      vs <- lapply(group, function(g) rotLeftTangentFromMatrix(g %*% r, mBar))
      ps <- sapply(vs, function(v) {1 - f(sqrt(v %*% covarInv %*% v))})
      # If any of the ps are NA, then max will return NA.
      max(ps)
    }
    list(pvalue=pvalue, ms=ms, kappas=kappas, mBar=mBar, leftCovarInv=covarInv,
         q090=percs[[1]], q095=percs[[51]], q099=percs[[91]], q100=percs[[101]],
         burnin=burnin, collection=collection)
  }
}

#' MCMC of posterior distribution rho(S, eta | D) for wrapped trivariate normal distribution parameters.
#'
#' Implemented in C for speed. See Qiu et al. (2014). This function requires compilation of the C shared library orientationsForR.c.
#' @param rs A list of rotation matrices. The data set D.
#' @param group A list of rotation matrices. The symmetry group G.
#' @param numTerms A real number (non-negative integer). Controls how many terms are used in the asymptotic expansions in the Jeffreys prior for kappa.
#' @param numTuning A real number (non-negative integer). The tuning parameters are re-tuned every numTuning MCMC iterations, based on the acceptance rate since the last tuning.
#' @param numBurnin A real number (non-negative integer). The number of MCMC iterations in the burn-in phase.
#' @param numCollection A real number (non-negative integer). The number of MCMC iterations in the collection phase. Should not exceed 100,000,000, or else we might overflow 32-bit integers along the way.
#' @param numReport A real number (non-negative integer). The number of MCMC samples (M, kappa) to report. If 0, then none are reported. If numCollection (or greater), then all are reported.
#' @return A list with elements pvalue (R function from rotations to {NA} union [0, 0.1]), ms (the reported Ms), kappas (the reported kappas), mBar (rotation matrix, the mean of the collected Ms), leftCovarInv (the inverse covariance matrix of the collected Ms in the left-invariant tangent space at mBar).
oricWrappedTrivariateNormalMCMCInference <- function(rs, group, numTerms=10, numTuning=10000, numBurnin=100, numCollection=1000, numReport=10000) {
  # Check that the inputs are of the correct types.
  # !!nums are integer; rs are real
  if (numTuning * numCollection < numReport) {
    numReport <- numTuning * numCollection
    print("warning: oricWrappedTrivariateNormalMCMCInference: clamping numReport to numTuning * numCollection")
  }
  # Flatten rs into an array of n * 9 numbers, by column-major order. Same with group.
  flat <- c(simplify2array(rs))
  flatGroup <- c(simplify2array(group))
  # Get a huge vector of numbers from C, to be parsed in a moment.
  raw <- .Call("mcmcOrientationWrappedTrivariateNormalC",
               flat, flatGroup, numTerms, numTuning, numBurnin, numCollection, numReport)
  oricMCMCParsing(raw, group, numReport)
}

#' Bootstrapping of the Frechet mean.
#'
#' Similar to oriBootstrapInference, but implemented in C for speed, and offers different percentiles of Mahalanobis norm. This function requires compilation of the C shared library orientationsForR.c.
#' @param rs A list of rotation matrices. The data set.
#' @param group A list of rotation matrices. The symmetry group G.
#' @param numBoots A real number (non-negative integer). The number of bootstrap samples. Affects the memory requirements of the function.
#' @return A list with elements pvalue (R function from rotations to {NA} union [0, 0.1]), bootstraps (the bootstrapped means), center (rotation matrix, the mean of the bootstrapped means), leftCovarInv (the inverse covariance matrix at rBar), q090, q095, q099, q100 (percentiles of Mahalanobis norm).
oricBootstrapInference <- function(rs, group, numBoots=10000) {
  # Check that the inputs are of the correct types.
  # !!numBoots is integer; rs are real
  # Flatten rs into an array of n * 9 numbers, by column-major order. Same for group.
  flat <- c(simplify2array(rs))
  flatGroup <- c(simplify2array(group))
  # Get a huge vector of numbers from C, to be parsed in a moment.
  raw <- .Call("pvalueOrientationBootstrappingC", flat, flatGroup, numBoots)
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
    vs <- lapply(group, function(g) rotLeftTangentFromMatrix(g %*% r, mBar))
    ps <- sapply(vs, function(v) {1 - f(sqrt(v %*% covarInv %*% v))})
    # If any of the ps are NA, then max will return NA.
    max(ps)
  }
  list(pvalue=pvalue, bootstraps=ms, center=mBar, leftCovarInv=covarInv,
       q090=percs[[1]], q095=percs[[51]], q099=percs[[91]], q100=percs[[101]])
}



### PLOTTING ###

#' Equal-volume plot of orientations with an accompanying Kamb density level surface.
#' 
#' This function requires compilation of the C shared library orientationsForR.c.
#' @param points A list of rotation matrices.
#' @param group A list of rotation matrices. The symmetry group G.
#' @param multiple A real number (positive). Indicates which multiple of the standard deviation to plot; for example, mult is 12 in a 12-sigma plot.
#' @param k A real number (positive). A smoothing factor, which equaled 3 in the original paper of Kamb (1959).
#' @param degree A real number (0, 1, or 3). The degree of the weighting polynomial; higher degrees generally produce smoother plots.
#' @param numNonAdapt A real number (non-negative integer). The number of non-adaptive refinements. Time and space required are proportional to 8^(numNonAdapt + numAdapt), so don't make it too big.
#' @param numAdapt A real number (non-negative integer). The number of adaptive refinements after the non-adaptive ones. Time and space required are proportional to 8^(numNonAdapt + numAdapt), so don't make it too big.
#' @param colors A list of strings (colors). Used to color the points only.
#' @param ... Plotting options: boundaryAlpha, simplePoints, etc. See rotPlotBall for details. Options about curves are ignored.
#' @return NULL.
oricKambPlot <- function(points, group, multiple=6, k=3, degree=3, numNonAdapt=3, numAdapt=3, colors=c("white"), ...) {
  pointss <- Reduce(c, lapply(group, function(g) lapply(points, function(p) g %*% p)))
  raws <- rotcKambTriangles(pointss, multiple, k, degree, numNonAdapt, numAdapt, length(group))
  vs <- lapply(pointss, rotEqualVolumeFromMatrix)
  colorss <- Reduce(c, replicate(length(group), colors))
  rotNativeEqualVolumePlot(points=vs, colors=colorss, trianglesRaw=raws, ...)
}


