


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# In this exercise we use orientation dispersion measures to study how fabric develops in rock. We assume that the fabric reflects the orientations of rigid triaxial ellipsoidal clasts, which rotate under the dynamical theory of Jeffery (1922). The orientations form fantastical patterns, which can be analyzed using orientation dispersion measures. Our ellipsoid tutorial 5dispersion.R does a similar analysis for deformable ellipsoidal clasts.

# The data here are all synthetic, generated to study a mathematical model of deformation.



### PRELIMINARY WORK ###

# If you have just restarted R, then set the working directory to the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc. Execute the following line of code to load our custom library and its dependencies.
source("library/all.R")

# Suppose that we're studying rocks that contain clasts of a highly competent material. We want to understand how fabric develops in such rocks. We make the following simplifying assumptions:
# * The large-scale deformation of the rock is homogeneous and steady.
# * The clasts are so competent, compared to the surrounding rock, that we can treat them as rigid. As the rock around them deforms, they rotate but do not deform.
# * The clasts are ellipsoids. The ellipsoids are triaxial; no two of their axes are equal in length. All of the ellipsoids have the same axial ratios.
# * The clasts are not crammed together in the rock. They are separated by at least four radii or so, so that they do not interact.
# * The initial (pre-deformation) orientations of the clasts were uniformly random.
# Under these assumptions, each clast obeys the dynamics of Jeffery (1922). So we can run numerical simulations, starting from uniformly random orientations, under various transpressions, to observe how the clasts align or not. By comparing the results of such experiments to natural fabrics, a geologist can gain insight into the deformations that produced those fabrics.

# Actually we're going to assume that the deformation is monoclinic transpression (Fossen and Tikoff, 1993). And we're going to use axial ratios 1 : 2 : 3, for no good reason.



### SIMULATING ISOLATED DEFORMATIONS ###

# This function takes as input a set of initial clast orientations, a gamma (magnitude of simple shear), and a logK (magnitude of coaxial deformation). All clasts are assumed to have semi-axis lengths 1, 2, 3. The function produces as output a set of final clast orientations, after rotation according to the Jeffery (1922) theory.
fabricFinalOris <- function(initials, gamma, logK) {
  l <- defMonoclinicVGT(gamma, logK)
  finals <- lapply(initials, function(r) defLeftJeffery(r, c(1, 2, 3), l, 20))
  rBar <- oriFrechetMean(finals, group=oriLineInPlaneGroup, numSeeds=5)
  finals <- oriNearestRepresentatives(finals, rBar, group=oriLineInPlaneGroup)
  finals
}

# Generate 1000 ellipsoid orientations uniformly. They seem to be evenly dense in the equal-volume plot, because of that plot's equal-volume property.
fabric000000 <- rotUniform(1000)
oriEqualVolumePlot(fabric000000, group=oriLineInPlaneGroup, simplePoints=TRUE)

# Subject the initially uniform orientations to a little bit of deformation --- gamma == 0.5 and logK == -0.5. This may take a minute. The orientations are still quite spread-out, but they're starting to show some clumping.
fabric005005 <- fabricFinalOris(fabric000000, 0.5, -0.5)
oriEqualVolumePlot(fabric005005, group=oriLineInPlaneGroup, simplePoints=TRUE)

# Subject the initially uniform orientations to severe, mostly coaxial, transpression. This may take a minute. The ellipsoids are really clumping around a preferred orientation.
fabric005080 <- fabricFinalOris(fabric000000, 0.5, -8.0)
oriEqualVolumePlot(fabric005080, group=oriLineInPlaneGroup, simplePoints=TRUE)

# Subject the initially uniform orientations to transpression that is severe in both coaxial and simple shear components. This may take a minute. The ellipsoid orientations are spread out along a curve (girdled).
fabric105080 <- fabricFinalOris(fabric000000, 10.5, -8.0)
oriEqualVolumePlot(fabric105080, group=oriLineInPlaneGroup, simplePoints=TRUE)



### AUTOMATION ###

# Doing experiments with isolated values of (gamma, log k) will give us only so much insight. So let's automate the process. We need some way to measure how spread-out orientations are. The Frechet variance is an obvious candidate. Using the four preceding examples, let's get a sense for what its values mean. (I got 0.89, 0.74, 0.04, 0.38.)
oriMeanVariance(fabric000000, group=oriLineInPlaneGroup, numSeeds=5)
oriMeanVariance(fabric005005, group=oriLineInPlaneGroup, numSeeds=5)
oriMeanVariance(fabric005080, group=oriLineInPlaneGroup, numSeeds=5)
oriMeanVariance(fabric105080, group=oriLineInPlaneGroup, numSeeds=5)

# The following code creates a contour plot of variance as a function of gamma and log k. This block of code took 11 hours to execute on my laptop. If you accidentally execute it, then press the Stop button in the upper-right corner of the Console pane. Just move on to the next part of this tutorial...
fabricStartDate <- date()
fabricGammas <- seq(from=0, to=20, by=0.5)
fabricLogKs <- seq(from=-8, to=0, by=0.25)
fabricVariance <- function(gamma, logK) {
  print(c(date(), as.character(gamma), as.character(logK)))
  finals <- fabricFinalOris(fabric000000, gamma, logK)
  meanVar <- oriMeanVariance(finals, group=oriLineInPlaneGroup, numSeeds=5)
  meanVar$variance
}
fabricVariances <- outer(fabricGammas, fabricLogKs, Vectorize(fabricVariance))
fabricEndDate <- date()
c(fabricStartDate, fabricEndDate)
contour(x=fabricGammas, y=fabricLogKs, z=fabricVariances)

# Unfortunately, variance is a very simple measure of dispersion. It can't capture how dispersion differs in differing directions. The following code block creates four contour plots, of the four rotMeanScatter dispersion measures, as functions of gamma and log k. It will probably take hours to execute. If you activate it accidentally, then press the Stop button in the Console pane to cancel it. Just move on to the next part of this tutorial...
fabricStartDate <- date()
fabricGammas <- seq(from=0, to=20, by=0.5)
fabricLogKs <- seq(from=-8, to=0, by=0.25)
fabricScatters <- function(gamma, logK) {
  print(c(date(), as.character(gamma), as.character(logK)))
  finals <- fabricFinalOris(fabric000000, gamma, logK)
  meanVar <- oriMeanVariance(finals, group=oriLineInPlaneGroup, numSeeds=5)
  finals <- oriNearestRepresentatives(finals, meanVar$mean, group=oriLineInPlaneGroup)
  meanScatter <- rotMeanScatter(finals)
  meanScatter$values
}
fabricLambda1s <- matrix(0, length(fabricGammas), length(fabricLogKs))
fabricLambda2s <- matrix(0, length(fabricGammas), length(fabricLogKs))
fabricLambda3s <- matrix(0, length(fabricGammas), length(fabricLogKs))
fabricLambda4s <- matrix(0, length(fabricGammas), length(fabricLogKs))
for (j in 1:length(fabricLogKs))
  for (i in 1:length(fabricGammas)) {
    scatters <- fabricScatters(fabricGammas[[i]], fabricLogKs[[j]])
    fabricLambda1s[[i, j]] <- scatters[[1]]
    fabricLambda2s[[i, j]] <- scatters[[2]]
    fabricLambda3s[[i, j]] <- scatters[[3]]
    fabricLambda4s[[i, j]] <- scatters[[4]]
  }
fabricEndDate <- date()
c(fabricStartDate, fabricEndDate)
contour(x=fabricGammas, y=fabricLogKs, z=fabricLambda1s)
contour(x=fabricGammas, y=fabricLogKs, z=fabricLambda2s)
contour(x=fabricGammas, y=fabricLogKs, z=fabricLambda3s)
contour(x=fabricGammas, y=fabricLogKs, z=fabricLambda4s)

# Our PDF notes include plots of the first, second, and third + fourth scatter measures. Look there.

# When dealing with steady deformations like this, you can trade intensity for time. For example, running the (gamma, log k) = (3, -1) transpression for two units of time is equivalent to running the (6, -2) transpresion for one unit of time. In other words, to observe how fabric develops over time, make a straight-line transect across the plot, starting at the origin in the upper-left corner. As you move along the transect, the deformation moves through time, and the contours tell you how the variance changes over time.

# For example, progressive simple shearing refers to the horizontal (log k = 0) transect. The contours suggest that the ellipsoids become non-uniform, then return to near-uniformity, then repeat. This periodic, pulsating behavior has been noted by various authors.

# For example, focus on the bottom half of the plot. lambda3 + lambda4 is small, and all valleys in lambda1 are compensated for by ridges in lambda2. Hence the ellipsoids are either clustered or girdled; they do not exhibit any other behaviors. So if a geologist observes other behaviors in the field, then she can conclude that log k is small.

# The same calculations could obviously be done for axial ratios other than 1 : 2 : 3. They could be done on sets of ellipsoids of varying axial ratios. They could be done with more complicated deformations. They could also be done on spheroids or deformable ellipsoids.



### WHAT HAVE WE LEARNED? ###

# Having a way to quantify orientation dispersion lets you automate the study of how fabrics develop.



### SEE ALSO ###

# We discuss the deformable case in ellipsoid tutorial 5dispersion.R.


