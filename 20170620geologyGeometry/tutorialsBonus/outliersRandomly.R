


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# This tutorial explains how we generated the synthetic outlier example used in direction tutorial 4notDirections.R and ellipsoid tutorial 4plottingVectors.R. Because the code uses random numbers, the results are different every time. The outlier is always similar to the blue data in one aspect and to the green data in another.

# This tutorial doesn't use any real data. It just generates synthetic data.



### PRELIMINARY WORK ###

# If you have just restarted R, then set the working directory to the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc. Execute the following line of code to load our custom library and its dependencies.
source("library/all.R")



### FOLD ORIENTATIONS ###

# First we pick lines l1, p1, l2, p2 such that (A) l1 and p1 are perpendicular, (B) l2 and p2 are perpendicular, and (C) l1 and p2 are perpendicular. That's the key.
synthL1 <- lineUniform()
synthP1 <- lineUniform()
synthP1 <- rayNormalized(synthP1 - (synthP1 %*% synthL1) * synthL1)
synthP2 <- lineUniform()
synthP2 <- rayNormalized(synthP2 - (synthP2 %*% synthL1) * synthL1)
synthL2 <- lineUniform()
synthL2 <- rayNormalized(synthL2 - (synthL2 %*% synthP2) * synthP2)

# Then we generate a curve of pairs stretching from (l1, p1) to (l2, p2).
synthR1 <- rotProjectedMatrix(rbind(synthP1, synthL1, cross(synthP1, synthL1)))
synthR2 <- rotProjectedMatrix(rbind(synthP2, synthL2, cross(synthP2, synthL2)))
synthR3 <- rotProjectedMatrix(rbind(synthP2, synthL1, cross(synthP2, synthL1)))
synthRs <- rotGeodesicPoints(synthR1, synthR2, 10)
synthRss <- list()
for (r in synthRs)
  synthRss <- c(synthRss, rotIsotropicFisher(r, 8, 5))

# Color roughly the first half blue, the second half green, and the last one red.
synthNumBlue <- floor(length(synthRss) / 2)
synthColors <- c(replicate(synthNumBlue, "blue"), replicate(length(synthRss) - synthNumBlue, "green"), "red")

# Here is the equal-area plot. The red circle should be similar to the blue circles. The red square should be similar to the green squares.
synthRs <- c(synthRss, list(synthR3))
synthPs <- lapply(synthRs, function(r) r[1,])
synthLs <- lapply(synthRs, function(r) r[2,])
lineEqualAreaPlot(points=c(synthPs, synthLs), colors=synthColors,
                  shapes=c(replicate(length(synthRss), "c"), "c", replicate(length(synthRss), "s"), "s"))

# Here is the equal-volume plot. Depending on how the random numbers went, the red point might or might not be obvious as an outlier.
oriEqualVolumePlot(synthRs, group=oriLineInPlaneGroup, colors=synthColors)

# The quickest way to try again is to copy and paste all of this code (lines 26-55) into the RStudio Console, press Return, inspect the results, press up-arrow and Return in the Console to run the code again, inspect the results, and so on.



### ELLIPSOIDS ###

# Make a swath of (9 + 1) * 3 = 30 orientations.
rs <- rotGeodesicPoints(rotUniform(), rotUniform(), 9)
rs <- unlist(lapply(rs, function(r) rotIsotropicFisher(r, 5, 3)), recursive=FALSE, use.names=FALSE)

# Make a swath of 30 octahedral shear strains and a swath of 30 Lode's parameters.
logEss <- seq(from=runif(1, -1, 1), to=runif(1, -1, 1), length.out=10)
ess <- as.numeric(sapply(logEss, function(logEs) exp(rnorm(3, logEs, 0.1))))
tanNu <- runif(1, -5, 5)
tanNus <- seq(from=tanNu, to=-tanNu, length.out=10)
nus <- as.numeric(sapply(tanNus, function(tanNu) atan(rnorm(3, tanNu, 1)) / (pi / 2)))

# Mate the orientations with the shapes to make normalized ellipsoids.
logAs <- lapply(1:30, function(i) ellLogAFromVEsNu(c(1, ess[[i]], nus[[i]])))
ellipsoids <- lapply(1:30, function(i) ellEllipsoidFromRotationLogA(rs[[i]], logAs[[i]], doNormalize=TRUE))

# Make the outlier by mating the first orientation with the last shape, with a small bit of noise.
r <- rotIsotropicFisher(rs[[1]], 10)
logA <- logAs[[30]] + rnorm(3, 0, 0.01)
outlier <- ellEllipsoidFromRotationLogA(r, logA, doNormalize=TRUE)
ellipsoids <- c(ellipsoids, list(outlier))

# The red orientation is similar to the blue orientations, but the red shape is similar to the green shapes.
# Depending on your particular random example, the outlier might or might not be obvious in the vector pair plots.
colors <- c(replicate(15, "blue"), replicate(15, "green"), "red")
ellEqualVolumePlot(lapply(ellipsoids, function(ell) ell$rotation),
                   lapply(ellipsoids, function(ell) ell$a),
                   colors=colors)
ellHsuNadaiPlot(lapply(ellipsoids, function(ell) ell$logA), colors=colors)
ellPairsPlot(lapply(ellipsoids, function(ell) ell$vector), colors=colors)

# Here is the fastest way to try it again: Copy and paste all of those lines of code (lines 32-60) into the R console, and press Return. Inspect the results. To do it again, click in the R console, press the up-arrow key, and press Return.



### WHAT HAVE WE LEARNED? ###

# When your data are high-dimensional, you can't view all of the dimensions at once. You can rig up a 'hidden outlier' example by making it similar to half of the data in some dimensions and the other half of the data in the other dimensions.


