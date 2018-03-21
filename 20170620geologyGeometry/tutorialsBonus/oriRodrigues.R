


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# This exercise demonstrates the Rodrigues plot of orientational data. Like the equal-volume plot, the Rodrigues plot is a distortion of the angle-axis plot. The Rodrigues plot enjoys some convenient geometric properties. It is apparently popular in crystallography and materials science.

# Data sets used include foliation-lineation pairs from the western Idaho shear zone (Giorgis and Tikoff, 2004), slickenside orientations from the Troodos ophiolite, Cyprus (Titus et al.), quartz orientations from the Moine thrust zone, Scotland (Strine and Wojtal, 2004; Michels et al., 2015), and synthetic data.



### PRELIMINARY WORK ###

# If you have just restarted R, then set the working directory to the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc. Execute the following line of code to load our custom library and its dependencies.
source("library/all.R")



### EXAMPLES ###

# Much like the equal-angle and equal-volume plots, the Rodrigues plot is a distortion of the angle-axis plot. Instead of plotting the rotation with angle a and axis u at the point a * u, we plot it at the point tan(a / 2) * u. Thus the space of rotations plots as a ball of infinite radius, filling out all of 3D space. When all symmetric copies of an orientation are shown, the ones near the origin often appear very compressed together, so you have to zoom in!

# Let's revisit our hidden outlier example from an earlier section. Most of the data are crunched up against the origin, so using the simplePoints option is a good idea.
outlierData <- geoDataFromFile("data/synthFoldsOutlier.csv")
oriRodriguesPlot(outlierData$rotation, group=oriLineInPlaneGroup, simplePoints=TRUE)

# Here are the foliation-lineations from the western Idaho shear zone (Giorgis and Tikoff, 2004).
wiszData <- geoDataFromFile("data/wiszFollins.tsv")
oriRodriguesPlot(wiszData$rotation, group=oriLineInPlaneGroup, simplePoints=TRUE)

# Here are 20 slickenside orientations (2-fold symmetry) from a single field site in Cyprus (Titus et al., unpublished).
slickData <- geoDataFromFile("data/cyprusSlicks2008005.tsv")
slickRots <- lapply(slickData$rotation, geoPoleVorticityFromPoleHanging)
oriRodriguesPlot(slickRots, group=oriRayInPlaneGroup, simplePoints=TRUE)

# Here are 761 quartz orientations (6-fold symmetry) obtained by electron backscatter diffraction (EBSD) of a single grain (Michels et al., 2015) taken from the Moine thrust zone, Scotland (Strine and Wojtal, 2004).
michelsData <- geoDataFromFile("data/moine_one_grainABCxyz.tsv")
oriRodriguesPlot(michelsData$rotation, group=oriTrigonalTrapezohedralGroup, simplePoints=TRUE)



### A NICE GEOMETRIC PROPERTY ###

# In the equal-angle and equal-volume plots, any steady progressive rotation through I plots as a straight line through the origin. However, steady progressive rotations that do not pass through I plot as gentle curves that are not straight. In the Rodrigues plot, ALL steady progressive rotations plot as straight lines.
rodR <- rotUniform()
rodU <- rayUniform()
rodCurve <- lapply(0:360, function(deg) rotMatrixFromAxisAngle(c(rodU, deg * degree)) %*% rodR)
rotEqualAnglePlot(rodCurve)
rotEqualVolumePlot(rodCurve)
rotRodriguesPlot(rodCurve, simplePoints=TRUE)

# This kind of plot has some other nice geometric properties. See Frank (1988). They don't include equal-angle or equal-volume properties.



### WHAT HAVE WE LEARNED? ###

# The Rodrigues plot is another way to view orientational data. When all symmetric copies of an orientation are shown, the ones near the origin appear very compressed.


