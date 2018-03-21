


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# In this tutorial, we learn how to quantify a steady, 'geodesic' tendency in rotation based on another, scalar variable --- in this case, distance from a spreading ridge. Our treatment includes a test of statistical significance.

# The data treated include rotations (orientations with trivial symmetry) inferred from paleomagnetic data from the Troodos ophiolite, Cyprus (Titus et al.).



### PRELIMINARY WORK ###

# If you have just restarted R, then set the working directory to the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc. Execute the following line of code to load our custom library and its dependencies.
source("library/all.R")

# Load some dike data and paleomagnetic data from the Troodos ophiolite, Cyprus. All of these field sites are west of the Solea Graben, which is believed to be a fossil spreading ridge. In this crude map, they are colored from red (west) to magenta (east).
cyprusData <- geoDataFromFile("data/cyprusPmagNW.csv")
names(cyprusData)
nrow(cyprusData)
plot(x=cyprusData$easting, y=cyprusData$northing, col=hues(cyprusData$easting))

# Here are the mean dike poles and mean paleomagnetic directions. The geoDataFromFile function has automatically inferred them from the strikes, dips, declinations, and inclinations. That function has also inferred rotations that we don't want, so delete them.
lineEqualAreaPlot(cyprusData$pole, colors=hues(cyprusData$easting))
rayEqualAreaPlot(cyprusData$direction, colors=hues(cyprusData$easting))
cyprusData$rotation <- NULL



### DEDUCING ROTATIONS ###

# We want to infer the rotation R that has occurred at each station. We want R to rotate the Troodos mean vector to the observed paleomagnetic direction. (By the way, this Troodos mean vector has a Fisher alpha95 confidence angle of 7 degrees.)
cyprusTroodosMeanVector <- geoCartesianFromTrendPlungeDeg(c(274, 36))

# Just knowing that R takes the Troodos mean vector to the observed paleomagnetic direction is not enough. There are infinitely many rotations that do that. So we need additional information or assumptions. In our orientation statistics paper, we impose the additional assumption that the dikes were initially vertical (Allerton and Vine, 1987). But in this tutorial let's do something simpler: Let's assume that the rotation was as small as possible.
cyprusData$rotation <- lapply(cyprusData$direction,
                              function(dir) rotSmallestRotationFromTwoRays(cyprusTroodosMeanVector, dir))
rotEqualAnglePlot(cyprusData$rotation, colors=hues(cyprusData$easting))

# If you view that plot from the correct angle, then a rainbow pattern appears. This leads us to suspect that there is a systematic tendency in the rotations, based on their easting. Here's the analogous equal-area plot of the axes without their angles.
cyprusData$axes <- lapply(cyprusData$rotation, function(r) rotAxisAngleFromMatrix(r)[1:3])
cyprusData$angles <- sapply(cyprusData$rotation, function(r) rotAxisAngleFromMatrix(r)[[4]])
rayEqualAreaPlot(cyprusData$axes, colors=hues(cyprusData$easting))

# In both of those plots there is a strong rainbow pattern, which suggests that the rotations vary systematically based on easting (which is distance to the spreading ridge). Our goal in this tutorial is to quantify that tendency using regression.



### INTERLUDE: SYNTHETIC EXAMPLE ###

# Compared to other rotational data sets, the Cyprus rotations look a little weird. They seem to plot on a particular surface in the equal-angle rotation plot, and their axes seem to plot on a particular curve in the equal-area plot. It's all because of the special way in which these rotations were chosen. They all take the same vector (the Troodos mean vector) to another given vector (the paleomagnetic direction at a station).

# To understand regression in a more typical (but prettier) case, let's work with a synthetic data set. Notice the rainbow pattern.
source("data/synthRotationsRegression.R")
rotEqualAnglePlot(synthRs, colors=hues(synthXs))

# Do geodesic regression. Ignore warnings. Check that error is 0 and minEigenvalue is positive. (If not, then something has gone wrong in the optimization that underlies the regression. See the documentation for rotGeodesicRegression.)
synthGeod <- rotGeodesicRegression(synthXs, synthRs)
synthGeod$error
synthGeod$minEigenvalue

# Visualize the result. You may choose to include 'spurs' showing how the data connect to their predictions.
rotGeodesicRegressionPlot(synthXs, synthRs, synthGeod)
rotGeodesicRegressionPlot(synthXs, synthRs, synthGeod, spurs=synthRs)

# To make a summary in words, convert the infinitesimal rotation synthGeod$m into a scaled rotation axis. Compute the trend, plunge, and rotation angle in degrees. We find that for each unit of synthX, there is 33 degrees of rotation about an axis trending 108 degrees and plunging 42 degrees. From start to finish, there is 245 degrees of rotation in total.
w <- rotVectorFromAntisymmetric(synthGeod$m)
geoTrendPlungeDegFromCartesian(rayNormalized(w))
geoTrendPlungeDegFromCartesian(lower(rayNormalized(w)))
sqrt(dot(w, w)) / degree
(max(synthXs) - min(synthXs)) * sqrt(dot(w, w)) / degree

# Check the quality of the fit in words. The R^2 statistic measures how much of the variation in the data is captured by the regression. If R^2 = 0, then the curve is no better at describing the data than the Frechet mean (a single point). If R^2 = 1, then the curve perfectly describes the data, in that the data lie along the curve. In this synthetic example, we get R^2 = 0.90, which is pretty high.
synthGeod$rSquared

# A separate issue is whether the result is significant. Does the regression capture a meaningful tendency in the data, or could the result be an artifact of random variability ('noise')? To address this question, we do a permutation test, much like the ones done in earlier directional statistics sections. We randomly permute the x-values, do the regression again, and observe how R^2 changes. After 10,000 such permutations (say), we let p be the fraction of permutations that produced greater R^2 than the original unpermuted data. If p < 0.05, then the result of the original regression is significant at the 95% level.

# For the sake of time, we'll do only 100 permutations, even though that is usually not enough to tease out a definitive p-value. When I did this with 1,000 permutations, I got p = 0 / 997 = 0. So we have found a statistically significant tendency in the data. That's not surprising, given this synthetic data set with a clear rainbow pattern in the equal-volume plot.
rSquareds <- rotGeodesicRegressionPermutations(synthXs, synthRs, numPerms=100)
length(rSquareds)
sum(rSquareds > synthGeod$rSquared)



### BACK TO CYPRUS ###

# So let's do the same thing with the Cyprus rotations. Check that error is 0 and minEigenvalue is positive.
cyprusGeod <- rotGeodesicRegression(cyprusData$easting, cyprusData$rotation, numSteps=1000)
cyprusGeod$error
cyprusGeod$minEigenvalue

# Visualize the result. Unsurprisingly, R^2 = 0.47 is not nearly as good as in our synthetic example. Intuitively, there is a lot of variation in the data in the directions perpendicular to the geodesic curve.
rotGeodesicRegressionPlot(cyprusData$easting, cyprusData$rotation, cyprusGeod)
rotGeodesicRegressionPlot(cyprusData$easting, cyprusData$rotation, cyprusGeod, spurs=cyprusData$rotation)
cyprusGeod$rSquared

# The regression reports 0.0054 degrees of rotation per m of easting, or 5.4 degrees of rotation per km of easting, about an axis trending 190 degrees and plunging 28 degrees, for a total of 63 degrees of rotation across the field area.
w <- rotVectorFromAntisymmetric(cyprusGeod$m)
geoTrendPlungeDegFromCartesian(rayNormalized(w))
sqrt(dot(w, w)) / degree
(max(cyprusData$easting) - min(cyprusData$easting)) * sqrt(dot(w, w)) / degree

# Here's the permutation test. When I did this with 1,000 permutations, I got p = 0 / 1000 = 0. So the dependence of rotation on easting is statistically significant.
rSquareds <- rotGeodesicRegressionPermutations(cyprusData$easting, cyprusData$rotation, numPerms=100)
length(rSquareds)
sum(rSquareds > cyprusGeod$rSquared)

# We can use the regression to predict or extrapolate. For example, the Solea graben is near easting 484765. Plug that into the regression curve to get a rotation. Apply that rotation to the Troodos mean vector, to predict what future paleomagnetic measurements at the graben should be.
cyprusPredRot <- cyprusGeod$prediction(484765)
cyprusPredRot
geoTrendPlungeDegFromCartesian(cyprusPredRot %*% cyprusTroodosMeanVector)

# Whether it is geologically reasonable to extrapolate into the graben is another question!



### WHAT HAVE WE LEARNED? ###

# Geodesic regression can pick out simple tendencies in rotational/orientational data depending on another variable.



### SEE ALSO ###

# In tutorialsBonus, rotKernelRegression.R deals with the relative motions of the Farallon and Pacific plates. Their motion is not well-described by geodesic regression, but kernel regression captures it nicely.


