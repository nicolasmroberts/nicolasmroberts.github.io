


# Copyright 2016-2017 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# This tutorial is an introduction to the basics of ellipsoids: How to visualize them, how to describe their orientation, size, and shape, etc. We introduce the crucial idea of degrees of freedom. There are hints that some of the standard ways of conceptualizing ellipsoids are going to lead to mathematical/statistical trouble (in the next exercise).

# The data used here include spinel grains obtained by X-ray computed tomography of rocks from New Caledonia (Titus et al., Chatzaras) and small synthetic data sets.



### PRELIMINARY WORK ###

# If you have just restarted R, then set the working directory to the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc. Execute the following line of code to load our custom library and its dependencies.
source("library/all.R")

# An ellipsoid has three perpendicular axes. Typically, but not always, they are all of different lengths. In this code, basicA is a vector of three semi-axis lengths. The rows of basicR are the semi-axis directions. The plot is colored so that the x-axis is red, the y-axis is green, and the z-axis is blue.
basicA <- c(5, 3, 1.5)
basicR <- rotUniform()
basicR
ellEllipsoidPlot(rots=list(basicR), as=list(basicA))

# Here are two ellipsoids, with the same semi-axis lengths but different orientations. In the plot, which ellipsoid goes with which orientation?
basicA <- c(5, 3, 1.5)
basicR <- rotUniform()
basicR
basicQ <- rotUniform()
basicQ
ellEllipsoidPlot(rots=list(basicR, basicQ), as=list(basicA, basicA))



### DEGREES OF FREEDOM ###

# Remember that choosing a rotation amounts to choosing three numbers: three Euler angles, three components of a scaled axis-angle vector, etc. We say that rotations have three 'degrees of freedom'. Similarly, directions have two degrees of freedom, such as trend and plunge. Ellipsoids have six degrees of freedom: three for the orientation, and three to describe the semi-axis lengths.

# A 'spheroid' is an ellipsoid in which two of the semi-axis lengths are equal. A spheroid is also called an 'ellipsoid of revolution'. Here is a visual example of two different orientations that produce the same spheroid. So do spheroids have six degrees of freedom?
basicA <- c(5, 5, 1)
basicR <- rotUniform()
basicR
basicQ <- rotMatrixAboutZ(20 * degree) %*% basicR
basicQ
ellEllipsoidPlot(rots=list(basicR, basicQ), as=list(basicA, basicA))



### MEASURING SIZE-SHAPE ###

# Thus far we've described size-shape in terms of the semi-axis lengths a1, a2, a3. But there are many other descriptors that geologists use. Here are a few. Most of them are most easily expressed in terms of the natural logarithms of the semi-axis lengths.
# * Volume is sometimes denoted V. V >= 0.
# * Octahedral shear strain is denoted e_s. e_s >= 0. e_s = 0 for spheres.
# * Lode's parameter nu is undefined for spheres. Otherwise -1 <= nu <= 1. nu = -1 (1) for prolate (oblate) spheroids.
# * Jelinek's P is sometimes denoted P' or P_j. P >= 1. P = 1 for spheres.
# * Flinn's K is undefined for a sphere or prolate spheroid. Otherwise K >= 0. K = 0 for non-sphere oblate spheroids.
# Workers in X-ray computed tomography and anisotropy of magnetic susceptibility (AMS) have other measures of ellipsoid shape.
basicA <- c(5, 3, 1.5)
basicLogA <- log(basicA)
basicLogA
ellVolume(basicLogA)
ellOctahedralShearStrain(basicLogA)
ellLodeNu(basicLogA)
ellJelinekP(basicLogA)
ellFlinnK(basicLogA)

# Jelinek's P is redundant with e_s, in that they can be converted into one another by P = exp(sqrt(2) * e_s).
basicA <- c(5, 3, 1.5)
basicLogA <- log(basicA)
basicEs <- ellOctahedralShearStrain(basicLogA)
exp(sqrt(2) * basicEs)
ellJelinekP(basicLogA)

# Except for certain redundancies such as P and e_s, all of the size-shapes measures contain slightly different information about the ellipsoid. In theory, there is no reason to regard the semi-axis lengths or their logs as 'primary' and these other measures as 'secondary'. Almost any three measures of size-shape can be regarded as primary, with the other measures derivable from them. For example, our ellLogsFromVEsNu function computes the semi-axis logs from V, e_s, and nu. Here you can check that it works correctly.
basicA <- c(5, 3, 1.5)
basicLogA <- log(basicA)
basicLogA
basicVEsNu <- c(ellVolume(basicLogA), ellOctahedralShearStrain(basicLogA), ellLodeNu(basicLogA))
basicVEsNu
ellLogAFromVEsNu(basicVEsNu)



### NORMALIZATION ###

# In some geologic problems, ellipsoid volume is not relevant. For example, deformations are commonly assumed to be volume-preserving, so their finite strain ellipsoids all have volume 4 pi / 3, so volume is not a distinguishing feature of these ellipsoids. Or maybe you'd like to compare anisotropy of magnetic susceptibility (AMS) ellipsoids to finite strain ellipsoids, and their units are incomparable, so comparing volumes is meaningless, although comparing shapes and orientations is still meaningful.

# It is common to normalize ellipsoids so that they have the same volume, 4 pi / 3, as the unit sphere. Concretely, we divide each semi-axis length ai by the scalar (a1 * a2 * a3)^(1 / 3). When an ellipsoid is normalized in this way:
# * a1 * a2 * a3 = 1.
# * log(a1) + log(a2) + log(a3) = 0.
# * e_s, nu, P, and K have the same values as for the pre-normalization ellpsoid.
# * nu = 0 only for those ellipsoids that are 'plane-strain', meaning that the intermediate semi-axis length is 1.

# Normalization effectively removes one degree of freedom from size-shape. The normalization above affects only volume, leaving the other two degrees of freedom to describe only shape. E_s and nu are two popular choices for these degrees of freedom. Sometimes geologists speak of E_s as measuring 'strain' and nu as measuring 'shape'. Ignoring the fact that the ellipsoid may not be a strain ellipsoid at all, E_s is just as much about shape as nu is, geometrically speaking. For example, here are two ellipsoids with the same orientation, volume, and nu. They differ only in E_s. Do they have the same shape? I would say not.
basicA <- exp(ellLogAFromVEsNu(c(4 * pi / 3, 1.1, -0.5)))
basicB <- exp(ellLogAFromVEsNu(c(4 * pi / 3, 1.8, -0.5)))
basicR <- diag(c(1, 1, 1))
ellEllipsoidPlot(rots=list(basicR, basicR), as=list(basicA, basicB))

# In other areas of the geosciences, such as AMS or X-ray computed tomography, other normalizations are in use. So you have to be careful when comparing data sets. Also, be aware that some of these normalizations affect multiple size-shape measures --- V, E_s, nu, etc. --- all at once.



### PLOTTING SHAPE ###

# The Hsu-Nadai plot is a way of visualizing ellipsoid shapes. The radial coordinate is e_s. The left boundary of the plot is where nu = -1 (prolate spheroids) and the right boundary is where nu = 1 (oblate spheroids). The left and right boundaries meet where e_s = 0 (spheres). Match up these three verbal descriptions with the corresponding points in the plot.
basicLogAs <- list(c(3, 1, -4), c(2, 2, 1), c(5, 1, 1))
ellHsuNadaiPlot(basicLogAs)

# Here's the same plot, with a third dimension showing log-volume. I've never seen this plot used by geologists, perhaps because publishing 3D images in papers is not easy. But in the age of computers there is no reason not to make such 3D plots, at least for exploring your data.
ellHsuNadaiScalarPlot(basicLogAs, log(sapply(basicLogAs, ellVolume)))

# This R library also has ellFlinnPlot, ellLogFlinnPlot, ellJelinekPlot, etc.



### PLOTTING ORIENTATION ###

# Here are normalized ellipsoids measuring shape preferred orientation of spinel clasts at 31 field sites in New Caledonia. The clasts were measured by X-ray computed tomography and then averaged at each station using ellMean. The file newcalOPXSpinelSPO.R loads these data into a list ncSpinels.
source("data/newcalOPXSpinelSPO.R")
length(ncSpinels)
ncSpinels[[1]]

# As we learned earlier, ellipsoid orientations can be viewed as rotations subject to line-in-plane symmetry. So let's plot the orientations, colored by easting. Frankly they're a mess.
ncColors <- hues(sapply(ncSpinels, function(station) station$en[[1]]))
ncRotations <- lapply(ncSpinels, function(station) station$spinel$rotation)
ncAs <- lapply(ncSpinels, function(station) station$spinel$a)
ellEqualVolumePlot(ncRotations, ncAs, colors=ncColors)

# Here's an equal-area plot of the orientations: circles for short axes, triangles for intermediate axes, and squares for long axes. Of course, you can't tell which axes go with which.
ellEqualAreaPlot(ncRotations, ncAs, colors=ncColors)

# Here are the shapes too. A few of the ellipsoids are close to spheroids, and hence their orientations are dubious in the first place.
ncLogAs <- lapply(ncSpinels, function(station) station$spinel$logA)
ellHsuNadaiPlot(ncLogAs, colors=ncColors)



### WHAT HAVE WE LEARNED? ###

# Ellipsoids have five (if normalized) or six (if not) degrees of freedom, including three for orientation and two or three for size-shape. There are various ways to measure and visualize these degrees of freedom.


