


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# This exercise describes how to measure the dispersion of an ellipsoid data set using the covariance of the ellipsoid vectors. It then works through a detailed application: Understanding how fabric develops in a rock that contains deformable ellipsoidal clasts. After some examples, we use dispersion measures to 'automate' the process of detecting fabric development.

# For data we will use a combined data set from New Caledonia: orthopyroxene shape preferred orientation (SPO) ellipsoids (Titus et al.) and spinel X-ray computed tomography ellipsoids (Chatzaras). We also do numerical simulations of how deformation affects deformable ellipsoidal clasts (Eshelby, 1957; Bilby et al., 1975).



### PRELIMINARY WORK ###

# If you have just restarted R, then set the working directory to the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc. Execute the following line of code to load our custom library and its dependencies.
source("library/all.R")

# Load data about shape preferred orientation (SPO) in New Caledonia. There are 34 SPO ellipsoids derived from field observations of orthopyroxene (Titus et al.). Then there are 31 SPO ellipsoids derived from X-ray computed tomography of spinel (Vasilis Chatzaras). The spinel clast ellipsoids at each station have been averaged into a single SPO using ellMean, as in ellipsoid exercise 3mean.R.
source("data/newcalOPXSpinelSPO.R")
length(ncOrthos)
length(ncSpinels)

# Examine some plots, wih the orthopyroxene SPOs in red and the spinel SPOs in cyan. Which do you think is more dispersed?
ellHsuNadaiPlot(
  c(lapply(ncOrthos, function(ell) ell$ortho$logA), lapply(ncSpinels, function(ell) ell$spinel$logA)),
  colors=c(replicate(length(ncOrthos), "red"), replicate(length(ncSpinels), "cyan")))
ellEqualVolumePlot(
  c(lapply(ncOrthos, function(ell) ell$ortho$rotation), lapply(ncSpinels, function(ell) ell$spinel$rotation)),
  c(lapply(ncOrthos, function(ell) ell$ortho$logA), lapply(ncSpinels, function(ell) ell$spinel$logA)),
  colors=c(replicate(length(ncOrthos), "red"), replicate(length(ncSpinels), "cyan")))

# Does this confirm your answer to the preceding question?
ellPairsPlot(
  c(lapply(ncOrthos, function(ell) ell$ortho$vector), lapply(ncSpinels, function(ell) ell$spinel$vector)),
  colors=c(replicate(length(ncOrthos), "red"), replicate(length(ncSpinels), "cyan")))



### MEASURING DISPERSION ###

# When ellipsoids are converted into vectors, you can do all kinds of multivariate statstics on them. If you wanted to measure the dispersion of a vector data set, a common approach would be to compute the sample variance matrix. This is a symmetric matrix whose eigenvectors describe perpendicular directions of dispersion (as columns of the eigenvector matrix) and whose associated eigenvalues describe the corresponding amounts of dispersion. Here's the variance of the orthopyroxene SPOs.
ncOrthoVariance <- ellVariance(lapply(ncOrthos, function(ell) ell$ortho$vector))
ncOrthoVariance
eigen(ncOrthoVariance)

# We could do the same computation for the spinel SPOs. But let's just focus on the eigenvalues. Does the following comparison confirm or refute your earlier opinion about which set of SPOs was more dispersed?
ellVarianceScalars(lapply(ncOrthos, function(ell) ell$ortho$vector))
ellVarianceScalars(lapply(ncSpinels, function(ell) ell$spinel$vector))

# In orientation statistics, we did a technique called principal component analysis, which produced principal directions and magnitudes of dispersion in a data set. Well, the same idea works for ellipsoid vectors (even more naturally, actually, because they form a vector space). You can access it through ellPrincipalComponentAnalysis if you like. But let's get to something more interesting.



### APPLICATION: FABRIC DEVELOPMENT ###

# Suppose that you're studying a rock unit that contains some clasts of a more competent material. The rock is deformed, and the clasts within it are deformed. You wish to understand how the deformed clasts contribute to the fabric of the rock. So you set up a numerical experiment, making several simplifying assumptions:
# * Before deformation, the clast ellipsoids all begin with the same size-shape. Their orientations are uniformly random.
# * The clasts and the surrounding rock are linear-viscous. All of the clasts have the same viscosity, but the host rock may have another viscosity.
# * The clasts are not tightly packed. They don't interact with each other during deformation.
# * The ambient deformation is steady, homogeneous, and volume-preserving.
# You simulate the deformation of the clasts using the dynamical theory of Eshelby (1957) and Bilby et al. (1975). Then you inspect the results, to understand the fabric.

# To implement this plan, let's make some specific choices: monoclinic transpression with gamma = 5 and logK = -1. Viscosity ratio 100 between the clasts and the host. Initial axial ratios 18.33 : 8.96 : 5.21, based on averaging the axial ratios of the unit cells of enstatite and ferrosillite. (We're going for orthopyroxene forming in a magma chamber and then deforming.)
fabricL <- defMonoclinicVGT(gamma=5, logK=-1)
fabricViscRatio <- 100
fabricLogA <- log(c(5.21, 8.96, 18.33))

# Generate a synthetic data set of 1,000 normalized ellipsoids, uniformly randomly oriented.
fabricInitials <- replicate(1000, ellEllipsoidFromRotationLogA(rotUniform(), fabricLogA, doNormalize=TRUE), simplify=FALSE)

# Inspect them in various plots. Compute the dispersion.
ellEqualVolumePlot(lapply(fabricInitials, function(ell) ell$rotation),
                   lapply(fabricInitials, function(ell) ell$a), simplePoints=TRUE)
ellEqualAreaPlot(lapply(fabricInitials, function(ell) ell$rotation),
                 lapply(fabricInitials, function(ell) ell$a))
ellHsuNadaiPlot(lapply(fabricInitials, function(ell) ell$logA))
ellPairsPlot(lapply(fabricInitials, function(ell) ell$vector))
ellVarianceScalars(lapply(fabricInitials, function(ell) ell$vector))

# Simulate the deformation of the ellipsoids forward in time, from t = 0 to t = 1, using 10 time steps. (Using 20 time steps produces very similar results.) On my laptop this takes about 40 seconds. Inspect the final ellipsoids in various plots.
fabricL <- defMonoclinicVGT(gamma=5, logK=-1)
fabricFinals050010 <- lapply(
  fabricInitials,
  function(ell) ellEllipsoidFromTensor(defLeftEshelby(ell$tensor, fabricL, fabricViscRatio, 10), doNormalize=TRUE))
ellEqualVolumePlot(lapply(fabricFinals050010, function(ell) ell$rotation),
                   lapply(fabricFinals050010, function(ell) ell$a), simplePoints=TRUE)
ellEqualAreaPlot(lapply(fabricFinals050010, function(ell) ell$rotation),
                 lapply(fabricFinals050010, function(ell) ell$a))
ellHsuNadaiPlot(lapply(fabricFinals050010, function(ell) ell$logA))
ellPairsPlot(lapply(fabricFinals050010, function(ell) ell$vector))
ellVectorPlot(c(2, 4, 5), lapply(fabricFinals050010, function(ell) ell$vector))

# Compute the dispersion. Notice that the two higher dispersions have increased slightly, while the three lower dispersions have decreased greatly.
ellVarianceScalars(lapply(fabricFinals050010, function(ell) ell$vector))

# Here's a more extreme example. What happened to the dispersion? Does it match with what's happening in the plots?
fabricL <- defMonoclinicVGT(gamma=14, logK=-4)
fabricFinals140040 <- lapply(
  fabricInitials,
  function(ell) ellEllipsoidFromTensor(defLeftEshelby(ell$tensor, fabricL, fabricViscRatio, 10), doNormalize=TRUE))
ellEqualVolumePlot(lapply(fabricFinals140040, function(ell) ell$rotation),
                   lapply(fabricFinals140040, function(ell) ell$a), simplePoints=TRUE)
ellEqualAreaPlot(lapply(fabricFinals140040, function(ell) ell$rotation),
                 lapply(fabricFinals140040, function(ell) ell$a))
ellHsuNadaiPlot(lapply(fabricFinals140040, function(ell) ell$logA))
ellPairsPlot(lapply(fabricFinals140040, function(ell) ell$vector))
ellVectorPlot(c(1, 4, 5), lapply(fabricFinals140040, function(ell) ell$vector))
ellVarianceScalars(lapply(fabricFinals140040, function(ell) ell$vector))



### APPLICATION CONTINUED: AUTOMATION ###

# We could investigate isolated examples of monoclinic transpression forever. But we can glean broader insights by automating the process: Deform the ellipsoids for many combinations of (gamma, log k), measure the final dispersion at each one, and make contour plots of those dispersions.

# This code computes about 700 examples. It takes my laptop 6 or 7 hours. (So I run it overnight, while I'm sleeping.) If you accidentally start it, then press the Stop button in RStudio's Console pane to stop it. We'll view precomputed results below...
fabricScalars <- function(gamma, logK, r, n) {
  fabricL <- defMonoclinicVGT(gamma=gamma, logK=logK)
  tryCatch(
    ellVarianceScalars(lapply(fabricInitials,
      function(ell) ellEllipsoidFromTensor(defLeftEshelby(ell$tensor, fabricL, r, n), doNormalize=TRUE)$vector)),
    error=function(e) c(NA, NA, NA, NA, NA))
}
fabricGammas <- seq(from=0, to=20, by=0.5)
fabricLogKs <- seq(from=0, to=-15, by=-1)
for (j in 1:length(fabricLogKs)) {
  firsts <- replicate(length(fabricGammas), 0)
  seconds <- replicate(length(fabricGammas), 0)
  thirds <- replicate(length(fabricGammas), 0)
  fourths <- replicate(length(fabricGammas), 0)
  fifths <- replicate(length(fabricGammas), 0)
  for (i in 1:length(fabricGammas)) {
    print(c(as.character(fabricGammas[[i]]), as.character(fabricLogKs[[j]]), date()))
    tuple <- fabricScalars(fabricGammas[[i]], fabricLogKs[[j]], 100, 10)
    firsts[[i]] <- tuple[[1]]
    seconds[[i]] <- tuple[[2]]
    thirds[[i]] <- tuple[[3]]
    fourths[[i]] <- tuple[[4]]
    fifths[[i]] <- tuple[[5]]
  }
  write.table(t(simplify2array(c(fabricLogKs[[j]], firsts))), file="fabricFirsts.tsv",
              append=TRUE, sep="\t", row.names=FALSE, col.names=FALSE)
  write.table(t(simplify2array(c(fabricLogKs[[j]], seconds))), file="fabricSeconds.tsv",
              append=TRUE, sep="\t", row.names=FALSE, col.names=FALSE)
  write.table(t(simplify2array(c(fabricLogKs[[j]], thirds))), file="fabricThirds.tsv",
              append=TRUE, sep="\t", row.names=FALSE, col.names=FALSE)
  write.table(t(simplify2array(c(fabricLogKs[[j]], fourths))), file="fabricFourths.tsv",
              append=TRUE, sep="\t", row.names=FALSE, col.names=FALSE)
  write.table(t(simplify2array(c(fabricLogKs[[j]], fifths))), file="fabricFifths.tsv",
              append=TRUE, sep="\t", row.names=FALSE, col.names=FALSE)
}

# This code loads the results from file and plots them.
fabricGammas <- seq(from=0, to=20, by=0.5)
fabricLogKs <- seq(from=-15, to=0, by=1)
fabricContourPlot <- function(fileName, levels=NULL) {
  a <- simplify2array(read.table(fileName, sep="\t", header=FALSE))
  if (is.null(levels))
    contour(x=fabricGammas, y=fabricLogKs,
            z=t(a[(from=seq(nrow(a), to=1, by=-1)),(2:ncol(a))]),
            xlab="gamma", ylab="log k")
  else
    contour(x=fabricGammas, y=fabricLogKs,
            z=t(a[(from=seq(nrow(a), to=1, by=-1)),(2:ncol(a))]),
            xlab="gamma", ylab="log k", levels=levels)
}
fabricContourPlot("data/fabricFirsts.tsv")
fabricContourPlot("data/fabricSeconds.tsv")
fabricContourPlot("data/fabricThirds.tsv")
fabricContourPlot("data/fabricFourths.tsv")
fabricContourPlot("data/fabricFifths.tsv")

# The precise results depend somewhat on the random initial ellipsoids (quantitatively, but not qualitatively). So that you can reconcile your future calculations with my contour plots, load the exact set of initial ellipsoids that was used to make them.
fabricInitials <- listFromDataFrame(geoDataFromFile("data/fabricInitials.csv"))



### APPLICATION CONTINUED: SIMPLE SHEARING ###

# Simple shears plot along the top edge of the plot (log k = 0). Let's check out the simple shear at (gamma, log k) = (12, 0). Its final ellipsoids are pretty close to the initial ellipsoids, with which we started. That's odd (?).
fabricL <- defMonoclinicVGT(gamma=12, logK=0)
fabric120000s <- lapply(
  fabricInitials,
  function(ell) ellEllipsoidFromTensor(defLeftEshelby(ell$tensor, fabricL, fabricViscRatio, 10), doNormalize=TRUE))
ellEqualVolumePlot(lapply(fabric120000s, function(ell) ell$rotation),
                   lapply(fabric120000s, function(ell) ell$a), simplePoints=TRUE)
ellEqualAreaPlot(lapply(fabric120000s, function(ell) ell$rotation),
                 lapply(fabric120000s, function(ell) ell$a))
ellHsuNadaiPlot(lapply(fabric120000s, function(ell) ell$logA))
ellPairsPlot(lapply(fabric120000s, function(ell) ell$vector))
ellVectorPlot(c(1, 4, 5), lapply(fabric120000s, function(ell) ell$vector))
ellVarianceScalars(lapply(fabric120000s, function(ell) ell$vector))

# In steady homogeneous deformations such as these, there is a time-intensity equivalence. For example, running a deformation with velocity gradient tensor L for a certain amount of time produces the same finite deformation as running 2 L for half as much time. Therefore, as a steady monoclinic transpression progresses, its cumulative effect traces out a line in the landscape, starting at the origin and proceeding on some constant heading through the landscape.

# In the case of steady simple shearing, you should imagine the deformation proceeding from left to right across the top of the plot. The ellipsoids start out uniform. Then they rotate and deform a bit. By gamma = 12 they have returned to something like their original state. There is a periodic, 'pulsating' behavior, which has been noted by every author who has studied problems like this.



### APPLICATION CONTINUED: COAXIAL DEFORMING ###

# Similarly, a steady coaxial deformation proceeds down the left side of the plot. Let's explore the valley around (gamma, log k) = (0, -10). It looks really different from the example above. Do the plots reflect the very low dispersion values?
fabricL <- defMonoclinicVGT(gamma=0, logK=-10)
fabric000100s <- lapply(
  fabricInitials,
  function(ell) ellEllipsoidFromTensor(defLeftEshelby(ell$tensor, fabricL, fabricViscRatio, 10), doNormalize=TRUE))
ellEqualVolumePlot(lapply(fabric000100s, function(ell) ell$rotation),
                   lapply(fabric000100s, function(ell) ell$a), simplePoints=TRUE)
ellEqualAreaPlot(lapply(fabric000100s, function(ell) ell$rotation),
                 lapply(fabric000100s, function(ell) ell$a))
ellHsuNadaiPlot(lapply(fabric000100s, function(ell) ell$logA))
ellPairsPlot(lapply(fabric000100s, function(ell) ell$vector))
ellVectorPlot(c(1, 4, 5), lapply(fabric000100s, function(ell) ell$vector))
ellVarianceScalars(lapply(fabric000100s, function(ell) ell$vector))

# By the way, log k = -3 means 95% shortening, log k = -4 means 98%, log k = -5 means 99%, log k = -10 means 99.99546%, and log k = -15 means 99.99997%. Studies of fabric development often go out to such extremes to approximate asymptotic behavior, which is easier to analyze mathematically than non-asymptotic behavior. If the viscosity ratio were smaller --- say, 5 instead of 100 --- then the clasts would deform more rapidly and we would not have to impose such extreme deformations.



### APPLICATION CONTINUED: THE REST OF THE LANDSCAPE ###

# Notice that the fourth and fifth dispersions drop to nearly nothing, as soon as there's a lot of shortening in the system (log k < -5).

# Based on the other three contour plots, I'm most interested in the progressive deformation that passes through (13, -15). Do you see why? Let's investigate. Yep, it's really different from the two previous examples.
fabricL <- defMonoclinicVGT(gamma=13, logK=-15)
fabric130150s <- lapply(
  fabricInitials,
  function(ell) ellEllipsoidFromTensor(defLeftEshelby(ell$tensor, fabricL, fabricViscRatio, 10), doNormalize=TRUE))
ellEqualVolumePlot(lapply(fabric130150s, function(ell) ell$rotation),
                   lapply(fabric130150s, function(ell) ell$a), simplePoints=TRUE)
ellEqualAreaPlot(lapply(fabric130150s, function(ell) ell$rotation),
                 lapply(fabric130150s, function(ell) ell$a))
ellHsuNadaiPlot(lapply(fabric130150s, function(ell) ell$logA))
ellPairsPlot(lapply(fabric130150s, function(ell) ell$vector))
ellVectorPlot(c(1, 4, 5), lapply(fabric130150s, function(ell) ell$vector), simplePoints=TRUE)
ellVarianceScalars(lapply(fabric130150s, function(ell) ell$vector))



### WHAT HAVE WE LEARNED? ###

# Statistics helps you measure dispersion of ellipsoids, which can help you quantify fabrics defined by ellipsoids.



### SEE ALSO ###

# In tutorialsBonus, oriRigidFabric.R does a similar study for rigid ellipsoids rotating under the dynamics of Jeffery (1922).


