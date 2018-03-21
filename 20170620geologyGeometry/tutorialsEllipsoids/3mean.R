


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# In this short exercise we compute the mean of a set of ellipsoids. It's not difficult. This notion of mean enjoys good properties in some special cases.

# For data we again use anisotropy of magnetic susceptibility (AMS) ellipsoids from rocks from the Troodos ophiolite, Cyprus (Titus et al.). We also use some synthetic data.



### PRELIMINARY WORK ###

# If you have just restarted R, then set the working directory to the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc. Execute the following line of code to load our custom library and its dependencies.
source("library/all.R")

# Load the Cyprus AMS data set again.
cyprusAMSData <- geoEllipsoidDataFromIRMFile("data/cyprus_AMS_groupF.tsv")

# Compute the mean by converting to vectors, taking the mean, and converting back.
cyprusAMSMean <- ellEllipsoidFromVector(arithmeticMean(cyprusAMSData$vector))
cyprusAMSMean

# It's that easy. But we've packaged the calculation into a single convenient function anyway.
cyprusAMSMean <- ellMean(cyprusAMSData$vector)
cyprusAMSMean

# Let's view the result in various plots. The mean is green, like Mean Joe Green.
ellEqualVolumePlot(c(cyprusAMSData$rotation, list(cyprusAMSMean$rotation)),
                   c(cyprusAMSData$a, list(cyprusAMSMean$a)),
                   colors=c(replicate(nrow(cyprusAMSData), "white"), "green"))
ellEqualAreaPlot(c(cyprusAMSData$rotation, list(cyprusAMSMean$rotation)),
                 c(cyprusAMSData$a, list(cyprusAMSMean$a)),
                 colors=c(replicate(nrow(cyprusAMSData), "black"), "green"))
ellHsuNadaiScalarPlot(c(cyprusAMSData$logA, list(cyprusAMSMean$logA)),
                      log(sapply(c(cyprusAMSData$logA, list(cyprusAMSMean$logA)), ellVolume)) / 1000,
                      colors=c(replicate(nrow(cyprusAMSData), "white"), "green"), es=0.05)

# A question to ponder: In the Hsu-Nadai plot, the mean ellipsoid does not plot in the middle of the other ellipsoids. It plots closer to the vertex. What does that mean, geometrically? And why should it be true?



### NICE PROPERTIES IN SPECIAL CASES ###

# If all of your ellipsoids have the same volume (perhaps because they've been volume-normalized), then the mean has that same volume.
cyprusAMSData <- geoEllipsoidDataFromIRMFile("data/cyprus_AMS_groupF.tsv", doNormalize=TRUE)
cyprusAMSMean <- ellMean(cyprusAMSData$vector)
ellHsuNadaiScalarPlot(c(cyprusAMSData$logA, list(cyprusAMSMean$logA)),
                      log(sapply(c(cyprusAMSData$logA, list(cyprusAMSMean$logA)), ellVolume)) / 100,
                      colors=c(replicate(nrow(cyprusAMSData), "white"), "green"), es=0.05)

# If all of your ellipsoids have the same orientation, then the mean has that same orientation.
synthR <- rotUniform()
synthLogAs <- replicate(10, sort(rnorm(3)), simplify=FALSE)
synthEllipsoids <- lapply(synthLogAs, function(logA) ellEllipsoidFromRotationLogA(synthR, logA))
synthEllipsoidsMean <- ellMean(lapply(synthEllipsoids, function(ell) ell$vector))
oriVariance(lapply(synthEllipsoids, function(ell) ell$rotation),
            synthEllipsoidsMean$rotation, group=oriLineInPlaneGroup)

# Also, if all of your ellipsoids have the same orientation, then the mean ellipsoid's semi-axis lengths are the geometric mean of the ellipsoids' semi-axis lengths. Equivalently, the mean ellipsoid's log-semi-axis lengths are the arithmetic mean of the ellipsoids' log-semi-axis lengths.
arithmeticMean(lapply(synthEllipsoids, function(ell) ell$logA))
synthEllipsoidsMean$logA

# Lastly, suppose that your ellipsoids are the finite strain ellipsoids of some steady homogeneous deformations. One can average such deformations by averaging their velocity gradient tensors. If the deformations are coaxial, then the mean of the finite strain ellipsoids is the finite strain ellipsoid of the mean deformation.



### WHAT HAVE WE LEARNED? ###

# Conceptualizing ellipsoids in terms of their ellipsoid vectors makes computations such as the sample mean quite simple.


