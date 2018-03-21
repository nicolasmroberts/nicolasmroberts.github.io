


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# Once we convert orientational data into rotations, we need a suitable way to plot those rotations. Perhaps the most popular system for describing rotations in geology is Euler angles. This exercise shows that Euler angle plots have some advantages over equal-area plots of orientational data. However, Euler angle plots can give an extremely distorted view of those data. So we use better plots in orientation tutorials such as 2plots.R.

# In addition to synthetic data sets, this tutorial uses three natural data sets: slickenside orientations from Cyprus (Titus et al.), quartz crystallographic orientations from the Moine thrust zone, Scotland (Strine and Wojtal, 2004; Michels et al., 2015), and foliation-lineation pairs from the western Idaho shear zone (Giorgis and Tikoff, 2004).



### PRELIMINARY WORK ###

# If you have just restarted R, then set the working directory to the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc. Execute the following line of code to load our custom library and its dependencies.
source("library/all.R")



### EULER ANGLES ###

# There are many different conventions for Euler angles. All of them share the same idea, of describing an overall rotation as a sequence of three rotations, one after the other. The convention we will use in this course is: First, we rotate about the geographic x-axis, through an angle alpha. Second, we rotate about the geographic z-axis through an angle beta. Third, we rotate about the geographic x-axis through an angle gamma.

# Here is an example worked out in detail, with alpha = 24 degrees, beta = 67 degrees, gamma = -12 degrees.
eulerFirst <- rotMatrixAboutX(24 * degree)
eulerSecond <- rotMatrixAboutZ(67 * degree)
eulerThird <- rotMatrixAboutX(-12 * degree)
eulerThird %*% eulerSecond %*% eulerFirst

# Our code offers easier functions for converting between Euler angles and rotation matrices.
eulerRot <- rotMatrixFromXZXAngles(c(-12, 67, 24) * degree)
eulerRot
rotXZXAnglesFromMatrix(eulerRot) / degree

# Converting from a rotation back to Euler angles requires us to choose angles according to certain conventions. So if you convert Euler angles to a matrix and then back to Euler angles, you might not get the angles with which you started.
rotXZXAnglesFromMatrix(rotMatrixFromXZXAngles(c(-23, -10, 41) * degree)) / degree

# Anyway, the idea in this tutorial is: Plot orientational data by converting to Euler angles and plotting those angles as points in 3D. The resulting plot measures (2 pi) x (pi) x (2 pi) or (360 degrees) x (180 degrees) x (360 degrees).



### OUTLIER EXAMPLE REVISITED ###

# Let's revisit our hidden outlier example from an earlier section.
outlierData <- geoDataFromFile("data/synthFoldsOutlier.csv")
lineEqualAreaPlot(c(outlierData$pole, outlierData$direction), colors=as.character(outlierData$color),
                  shapes=c(replicate(length(outlierData$pole), "c"), replicate(length(outlierData$direction), "s")))

# Let's plot the Euler angles for all of the plane-line pairs. Notice that there appear to be four copies of the data set, due to the 4-fold plane-line symmetry. The outlier is now clearly visible.
oriEulerAnglePlot(outlierData$rotation, group=oriLineInPlaneGroup)



### THREE NATURAL DATA SETS ###

# Here are 20 slickenside orientations (2-fold symmetry) from a single field site in Cyprus (Titus et al., unpublished). Notice that the data are few and spread out.
slickData <- geoDataFromFile("data/cyprusSlicks2008005.tsv")
slickRots <- lapply(slickData$rotation, geoPoleVorticityFromPoleHanging)
slickPoles <- lapply(slickRots, function(r) lower(r[1,]))
slickVorts <- lapply(slickRots, function(r) r[2,])
rayEqualAreaPlotTwo(slickPoles, slickVorts)
oriEulerAnglePlot(slickRots, group=oriRayInPlaneGroup)

# Here are 761 quartz orientations (6-fold symmetry) obtained by electron backscatter diffraction (EBSD) of a single grain (Strine and Wojtal, 2004; Michels et al., 2015). Notice that the data are numerous and tightly concentrated.
michelsData <- geoDataFromFile("data/moine_one_grainABCxyz.tsv")
miller100s <- lapply(michelsData$rotation, function(r) r[1,])
miller010s <- lapply(michelsData$rotation, function(r) r[2,])
miller001s <- lapply(michelsData$rotation, function(r) r[3,])
lineEqualAreaPlotThree(miller100s, miller010s, miller001s)
oriEulerAnglePlot(michelsData$rotation, group=oriTrigonalTrapezohedralGroup, simplePoints=TRUE)

# Orientation/rotation statistics has been used heavily in the crystallography community for processing EBSD data. Some of the methods are asymptotic, relying on large sample size or tight concentration. The quartz data above demonstrate that EBSD data often fulfill these requirements. But the slickenside data above demonstrate that structural field data often do not. So we have to be careful in choosing statistical methods that work for small, widely dispersed data sets.

# Here are the 23 foliation-lineations from the western Idaho shear zone (Giorgis and Tikoff, 2004). They look pretty crazy, due to a problem we'll describe next.
wiszData <- geoDataFromFile("data/wiszFollins.tsv")
oriEulerAnglePlot(wiszData$rotation, group=oriLineInPlaneGroup)



### GIMBAL LOCK ###

# Euler angles suffer from a defect called 'gimbal lock'. When the middle angle is any multiple of 180 degrees, the other two angles effectively collapse down to one angle. In these first three examples, notice that only the sum of gamma and alpha matters.
rotMatrixFromXZXAngles(c(53, 0, 8) * degree)
rotMatrixFromXZXAngles(c(41, 0, 20) * degree)
rotMatrixFromXZXAngles(c(24, 0, 37) * degree)

# In these next three examples, notice that only the difference of gamma and alpha matters.
rotMatrixFromXZXAngles(c(53, 180, 7) * degree)
rotMatrixFromXZXAngles(c(61, 180, 15) * degree)
rotMatrixFromXZXAngles(c(76, 180, 30) * degree)

# Geometrically, gimbal lock means that entire lines in Euler angle space represent a single rotation. On the beta == 0 and beta == 180 degrees boundary planes, there is an infinite-to-one correspondence between Euler angles and rotations. Consequently, the Euler angle plot is extremely distorted near these boundary planes. For example, consider the following synthetic data set of tightly concentrated alpha-quartz orientations (6-fold symmetry). Four of the symmetric copies appear tightly concentrated, as is correct. But the other two are 'smeared out' along the gimbal lock lines.
gimbalData <- geoDataFromFile("data/synthGimbalLock.csv")
oriEulerAnglePlot(gimbalData$rotation, group=oriTrigonalTrapezohedralGroup)

# Here are the same data in the equal-volume plot. There are no weird effects.
oriEqualVolumePlot(gimbalData$rotation, group=oriTrigonalTrapezohedralGroup)



### WHAT HAVE WE LEARNED? ###

# Although Euler angles might be familiar, they are not particularly well-behaved, and we should strongly consider using other conventions for describing rotations and orientations.


