


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# Now that we have a system for converting plane-line pairs into rotations (up to 4-fold symmetry), we want to be able to plot those rotations. In this tutorial we introduce two systems for plotting rotational/orientational data: equal-angle and equal-volume. Like the equal-angle and equal-area hemispherical plots, the equal-angle and equal-volume rotation plots depict the shape and density of data accurately, respectively. In fact, these plots can be regarded as hemispherical plots 'one dimension up'.

# Like the preceding tutorial, this one uses foliation-lineation data from the western Idaho shear zone (Giorgis and Tikoff, 2004). It also uses the synthetic data set from direction tutorial 4notDirections.R.



### PRELIMINARY WORK ###

# If you have just restarted R, then set the working directory to the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc. Execute the following line of code to load our custom library and its dependencies.
source("library/all.R")



### ANGLE-AXIS PLOT ###

# Before learning the equal-angle and equal-volume plots, it is helpful to know the angle-axis plot. Every rotation of 3D space can be conceptualized as an axis of rotation and an angle of rotation about that axis. The axis can be described as a unit vector. Rotation happens about that vector counterclockwise, according to the right-hand rule. With these conventions in place, angles between 0 and pi = 180 degrees suffice for describing all rotations.

# It is common to unit the axis vector u and the angle a into a scaled axis vector a * u. In fact, the animations shown in 1linesInPlanes.R included this scaled axis vector exactly. Here is that example again. We're just going to take the four scaled axes, copy them to a new 3D plot, and mark their endpoints.
wiszData <- geoDataFromFile("data/wiszFollins.tsv")
wiszR <- wiszData$rotation[[15]]
wiszUA <- rotAxisAngleFromMatrix(oriLineInPlaneGroup[[1]] %*% wiszR)
wiszV1 <- wiszUA[1:3] * wiszUA[[4]]
wiszUA <- rotAxisAngleFromMatrix(oriLineInPlaneGroup[[2]] %*% wiszR)
wiszV2 <- wiszUA[1:3] * wiszUA[[4]]
wiszUA <- rotAxisAngleFromMatrix(oriLineInPlaneGroup[[3]] %*% wiszR)
wiszV3 <- wiszUA[1:3] * wiszUA[[4]]
wiszUA <- rotAxisAngleFromMatrix(oriLineInPlaneGroup[[4]] %*% wiszR)
wiszV4 <- wiszUA[1:3] * wiszUA[[4]]
plot3D(points=list(wiszV1, wiszV2, wiszV3, wiszV4), radius=pi,
       curves=list(list(c(0, 0, 0), wiszV1), list(c(0, 0, 0), wiszV2), list(c(0, 0, 0), wiszV3), list(c(0, 0, 0), wiszV4)))

# This system of plotting rotations is called the 'angle-axis plot'. It depicts the space of rotations as a ball of radius pi = 180 degrees. Points interior to the ball are in one-to-one correspondence with rotations. For example, the origin of the ball corresponds to the trivial rotation I, and the point [30 degrees  0  0]^T corresponds to rotation about [1 0 0]^T through 30 degrees. Points on the boundary of the ball are antipodally identified, just like the points on the boundary of a hemispherical plot.

# Here's the same plot, without the lines from the origin, and with the boundary sphere emphasized. Remember: This is the plot of a single foliation-lineation pair.
oriAxisAnglePlot(list(wiszR), group=oriLineInPlaneGroup, boundaryAlpha=0.25)

# Here's the entire western Idaho shear zone foliation-lineation data set. Notice that there are four symmetric copies of the data. Two of the copies cross the boundary. In each of those, you have to imagine how the two halves connect --- much like in a hemispherical plot.
oriAxisAnglePlot(wiszData$rotation, group=oriLineInPlaneGroup)

# Here's the hidden outlier example from direction exercise 4noDirections.R. Again there are four copies. In each copy, the outlier is unmistakable.
outlierData <- geoDataFromFile("data/synthFoldsOutlier.csv")
oriAxisAnglePlot(outlierData$rotation, group=oriLineInPlaneGroup)

# The angle-axis plot does not have an equal-angle property, so it does not depict the shape of a data set faithfully. The angle-axis plot also does not have an equal-volume property, so it does not depict density of data faithfully.



### EQUAL-ANGLE AND EQUAL-VOLUME PLOTS ###

# It turns out that we can slightly distort the angle-axis plot to produce new plots that do have equal-angle or equal-volume properties. Here are the western Idaho shear zone foliation-lineation data set again.
oriEqualAnglePlot(wiszData$rotation, group=oriLineInPlaneGroup)
oriEqualVolumePlot(wiszData$rotation, group=oriLineInPlaneGroup)

# Here's the synthetic outlier data set again.
oriEqualAnglePlot(outlierData$rotation, group=oriLineInPlaneGroup)
oriEqualVolumePlot(outlierData$rotation, group=oriLineInPlaneGroup)

# You might want to view the angle-axis, equal-angle, and equal-volume plots side-by-side, to see how they compare. They are subtle distortions of each other, in much the same way that the equal-angle and equal-volume hemispherical plots are subtle distortions of each other.



### NOT EASY? ###

# In the rest of these tutorials, when we need to plot rotational/orientational data, we will usually employ the equal-angle or equal-volume plot. For example, the equal-volume plot will help us determine whether a data set is unimodal (one region of high density, centered on a mean) or bimodal (two regions of high density, so that the mean doesn't make sense). The equal-volume plot depicts distance and shape relationships imperfectly, but much better than an equal-area plot of planes and lines separately, for example.

# With all that said, you may still feel that equal-angle and equal-volume plots of rotational versions of structural data is all too strange. Here is my advice:

# A. Equal-angle and equal-area hemispherical plots were not obvious either, when you first saw them. They required practice. So do these plots. And they require more practice than the equal-area plot, because orientations are inherently more complicated than directions.

# B. Even if you were in love with the equal-angle and equal-volume plots, I would not recommend that you use them exclusively. Always view your data in multiple kinds of plot, to minimize your chances of overlooking important features of your data set.

# C. Even if you've already vowed that you'll never use any of these plots again, the rest of these tutorials still offer a variety of useful statistical methods that have little to do with plotting.



### WHAT HAVE WE LEARNED? ###

# You can plot rotational/orientational data, but the most obvious or easiest way may not be the best way in the long term.



### SEE ALSO ###

# The bonus exercise oriEulerAngles.R describes plotting rotational/orientational data in terms of their Euler angles. That plot is not nearly as well-behaved as the angle-axis plot, let alone the equal-angle and equal-volume plots.

# The bonus exercise oriRodrigues.R describes the Rodrigues plot, which is a different distortion of the angle-axis plot, popular in crystallography and materials science.

# The bonus exercise publicationPlots.R describes how to adjust these plots for publication (white backgrounds, etc.).


