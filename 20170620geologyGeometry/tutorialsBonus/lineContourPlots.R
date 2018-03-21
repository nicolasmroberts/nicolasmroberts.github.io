


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# Directional statistics software used in structural geology (e.g., Stereonet by Allmendinger and Cardozo, Orient by Vollmer) have a feature called Kamb contouring, which depicts the density of lines in an equal-area plot. And our library includes some basic Kamb contouring functionality too. But the underlying mechanism can be used to contour-plot any function of lines --- not just the Kamb density function. This tutorial gives some rudimentary examples, which might empower you to do something more interesting.

# For data we use dike poles from the Troodos ophiolite, Cyprus (Titus et al.), as well as synthetic data sets.



### PRELIMINARY WORK ###

# If you have just restarted R, then set the working directory to the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc. Execute the following line of code to load our custom library and its dependencies.
source("library/all.R")



### KAMB CONTOURING ###

# Load some dike poles from the Troodos ophiolite, Cyprus.
cyprusData <- geoDataFromFile("data/cyprus_inside_corner_field_data.tsv")
nrow(cyprusData)
lineEqualAreaPlot(cyprusData$pole)

# Kamb-contour with the default settings.
lineKambPlot(cyprusData$pole)

# There are various ways to control the smoothness of the contours. Most important is the numNonAdapt parameter, which controls the fineness of the mesh where the density is sampled. The default value is 4. For each increment of numNonAdapt, the contours get a little smoother, but the time and memory requirements go up by a factor of four. So don't set numNonAdapt too high.
lineKambPlot(cyprusData$pole, numNonAdapt=5)

# By default, lineKambPlot displays four contours: 3-sigma, 6-sigma, 9-sigma, and 12-sigma. You can change these numbers c(3, 6, 9, 12) to whatever you want using the multiples parameter.
lineKambPlot(cyprusData$pole, multiples=c(2, 4, 6, 8, 10, 12, 14))

# If you look closely at the bottom of that last plot, you see a tiny defect: a little 'star' of lines sitting along a contour. These stars appear where the mesh is not fine enough to resolve the contours. To get rid of them, try increasing numNonAdapt.
lineKambPlot(cyprusData$pole, multiples=c(2, 4, 6, 8, 10, 12, 14), numNonAdapt=5)



### GIRDLE ###

# The lineEqualAreaFunctionPlot function, which underlies lineKambPlot, lets you contour-plot any function you want. First you write the function. As input it should take a single line u (a unit 3D Cartesian vector). As output it should return a single number.

# For example, the Cyprus dike poles seem to form a girdle. From direction tutorial 1oneDirection.R, here's how to compute the pole of the girdle.
cyprusGirdlePole <- lineMeanScatter(cyprusData$pole)$vectors[,3]
geoTrendPlungeDegFromCartesian(lower(cyprusGirdlePole))

# Which data are within 10 degrees of the girdle? We can plot that region by plotting the lines that are 80 degrees away from the pole.
contourFunc <- function(u) {
  lineDistance(u, cyprusGirdlePole)
}
lineEqualAreaFunctionPlot(contourFunc, levels=c(80 * degree), points=cyprusData$pole)

# This example is not really a great illustration of lineEqualAreaFunctionPlot, because the contour consists of two small circles, which we could compute explicitly using other means.
lineEqualAreaPlot(cyprusData$pole, curves=list(raySmallCircle(cyprusGirdlePole, 80 * degree)))



### KINEMATIC VORTICITY OF TRICLINIC TRANSPRESSIONS ###

# How does the kinematic vorticity of triclinic transpression (Lin et al., 1998) depend on the relative movement direction of the rigid blocks bounding the shear zone? Here we consider a vertical, EW-striking shear plane and let u be the movement direction.
contourFunc <- function(u) {
  defKinematicVorticity(defTriclinicVGT(u))
}

# Again the contours happen to be small circles, because the kinematic vorticity depends only on the angle between u and the shear plane. But notice that the contours are more widely spaced near the shear plane. That's because kinematic vorticity asymptotically approaches 1, as u approaches the shear plane.
lineEqualAreaFunctionPlot(contourFunc, levels=seq(from=0, to=1, by=0.1))
lineEqualAreaFunctionPlot(contourFunc, levels=seq(from=0.9, to=1, by=0.01), numNonAdapt=5)



### DISTANCE FROM A DATA SET ###

# In this admittedly artificial example, we define contourFunc(u) to be the distance from u to the nearest Cyprus dike pole.
contourFunc <- function(u) {
  min(sapply(cyprusData$pole, lineDistance, u))
}

# So here are the points that are 10 degrees away from the Cyprus dike data set. Notice the star defects.
lineEqualAreaFunctionPlot(contourFunc, levels=c(10 * degree), points=cyprusData$pole)

# As you increase numNonAdapt, the stars get smaller, but they persist even at numNonAdapt = 6. They live in the 'cusps' of the contours, which are inherently difficult to resolve.
lineEqualAreaFunctionPlot(contourFunc, levels=c(10 * degree), numNonAdapt=6, points=cyprusData$pole)

# Here's a trippier example, like a brain coral.
lineEqualAreaFunctionPlot(contourFunc, levels=c(5 * degree, 10 * degree, 15 * degree, 20 * degree, 25 * degree),
                          numNonAdapt=5, points=cyprusData$pole)



### VORONOI DIAGRAM ###

# This artificial example displays a Voronoi diagram on the sphere. First we choose some 'seed' lines randomly. Think of them as capital cities of kingdoms. Each kingdom consists of all points that are closer to its capital than to any other capital. We define contourFunc to be nonnegative, and zero wherever two or more kingdoms touch.
contourSeeds <- lineUniform(11)
contourFunc <- function(u) {
  dists <- sort(sapply(contourSeeds, lineDistance, u))
  dists[[2]] - dists[[1]]
}

# Then we approximate the borders by asking where contourFunc is small.
lineEqualAreaFunctionPlot(contourFunc, levels=c(2 * degree), numNonAdapt=6, points=contourSeeds)

# You can make the borders tighter by asking for a smaller contour level. But then resolution-dependent defects often appear. Due to the unstable nature of the function we're contouring, the star defects are difficult to mitigate, even at high numNonAdapts.
lineEqualAreaFunctionPlot(contourFunc, levels=c(2 * degree, 1 * degree), numNonAdapt=6, points=contourSeeds)

# A question to ponder: How would you modify contourFunc, so that the resulting plot showed the triple junctions instead of the boundaries?



### WHAT HAVE WE LEARNED? ###

# This library includes some basic Kamb contouring of lines, and you can use the underlying mechanism to contour-plot any function of lines.


