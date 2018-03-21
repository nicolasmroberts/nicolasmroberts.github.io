


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# In earlier tutorials we converted plane-line pairs into rotations. In this tutorial we do the same for other orientational data, including plane-ray pairs and crystallographic orientations. Each orientational data type has a symmetry group, that causes a single orientation to be represented by not just one rotation but a set of them. We also discuss examples of true rotational data.

# Data used in this tutorial include slickenside straie from Cyprus (Titus et al.) and quartz crystallographic orientations (Strine and Wojtal, 2004; Michels et al., 2015).



### PRELIMINARY WORK ###

# If you have just restarted R, then set the working directory to the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc. Execute the following line of code to load our custom library and its dependencies.
source("library/all.R")



### FAULTS WITH SLIP DIRECTION ###

# Suppose that we're looking at a set of faults, where we know the 3D direction of slip along each fault --- for example, from slickenside straie. Here is a single slickenside orientation from Cyprus in strike-dip-rake format, where the rake describes the movement of the hanging wall. The first row of the matrix R is the pole to the fault plane. The second row is the slip direction (hanging wall movement direction).
slickR <- geoRotationFromStrikeDipRakeDeg(c(359, 39, 38))
slickR

# The problem with that format is that it behaves badly for vertical faults, where the roles of foot wall and hanging wall are ill-defined. So instead of the pole and the slip direction, we work with the pole and the vorticity vector. For example, a vertical fault has vorticity [0 0 1]^T if it's sinistral or [0 0 -1]^T if it's dextral.
slickR <- geoPoleVorticityFromPoleHanging(geoRotationFromStrikeDipRakeDeg(c(359, 39, 38)))
slickR

# We can negate the pole vector without changing its physical meaning. But we cannot negate the vorticity vector, because that would change the sense of slip. So we're not talking about a line in a plane. We're talking about a ray (directed line) in a plane. There is a 2-fold symmetry group, called oriRayInPlaneGroup in our software.
oriRayInPlaneGroup[[1]] %*% slickR
oriRayInPlaneGroup[[2]] %*% slickR

# Here are all 20 slickensides from that field site, shown in the equal-volume plot. Each datum appears twice, because of the 2-fold symmetry. So you can see 40 points in total. Unlike some of our earlier data sets, these data are so spread-out that you can't clearly separate the two symmetric copies.
slickData <- geoDataFromFile("data/cyprusSlicks2008005.tsv")
slickRots <- lapply(slickData$rotation, geoPoleVorticityFromPoleHanging)
oriEqualVolumePlot(slickRots, group=oriRayInPlaneGroup)



### BEDDING WITH RIPPLES? ###

# Here's a question to ponder. Suppose that our field area consists of a bunch of deformed sandstone. The bedding tilts in different directions as we move through the field area. The bedding contains ripples. For each bed, you know the younging directin (which way was 'up' when the bed was deposited). What is the symmetry group of this orientational data type? Warning: It depends on whether the ripples are symmetric or asymmetric.



### CRYSTALLOGRAPHIC ORIENTATIONS ###

# Sometimes structural geologists interact with mineralogy, for example when analyzing microstructures in thin section. Michels et al. (2015) analyzed orientations of quartz grains using electron backscatter diffraction (EBSD). The following file contains 761 of these orientations from a single grain from the Moine thrust zone (Strine and Wojtal, 2004). Each orientation is given as three unit Cartesian 3D vectors, indicating the geographic directions of Miller indices 100, 010, 001. Load all of them, but inspect just the first orientation.
michelsData <- geoDataFromFile("data/moine_one_grainABCxyz.tsv")
michelsData[1,]

# The crystallographic point group of alpha-quartz is the trigonal trapezohedral group, which has six elements. So here are the six rotations representing the first orientation.
oriTrigonalTrapezohedralGroup[[1]] %*% michelsData$rotation[[1]]
oriTrigonalTrapezohedralGroup[[2]] %*% michelsData$rotation[[1]]
oriTrigonalTrapezohedralGroup[[3]] %*% michelsData$rotation[[1]]
oriTrigonalTrapezohedralGroup[[4]] %*% michelsData$rotation[[1]]
oriTrigonalTrapezohedralGroup[[5]] %*% michelsData$rotation[[1]]
oriTrigonalTrapezohedralGroup[[6]] %*% michelsData$rotation[[1]]

# Here is the equal-volume plot of all 761 orientations. Each one appears six times due to the 6-fold symmetry. Notice that the data are tightly concentrated, so that the six copies are easily distinguishable. (Because so many points are being displayed, we use the 'simplePoints' option to make each point less computationally demanding.)
oriEqualVolumePlot(michelsData$rotation, group=oriTrigonalTrapezohedralGroup, simplePoints=TRUE)



### AN IMPORTANT SIDE REMARK ###

# Let's look at those slickenside and quartz orientations again. Put these plot windows side-by-side.
oriEqualVolumePlot(slickData$rotation, group=oriRayInPlaneGroup)
oriEqualVolumePlot(michelsData$rotation, group=oriTrigonalTrapezohedralGroup, simplePoints=TRUE)

# There are two crucial things to notice here:

# 1. Structural field data are often few in number --- 10s or 100s. Data derived from automated instruments are often much more numerous --- 100s or 1000s.

# 2. Structural field data are often widely dispersed. EBSD data are often tightly concentrated (partially due to the criterion used to detect grain boundaries).

# Large sample sizes and tight concentration both make statistics easier, because there are many classical methods that are asymptotic in sample size and concentration. Those methods are not appropriate for data sets such as the slickensides above. So we have to choose our methods carefully.



### ACTUAL ROTATIONS ###

# Thus far we have used rotations to describe the orientations of objects in 3D. But sometimes actual rotations show up in geology. They do not have any symmetry group --- or rather, the symmetry group G = {I} is trivial.

# For example, Engebretson et al. (1984) reconstructed the relative motions of tectonic plates adjacent to the Pacific. Their results were stated as a sequence of rotations about various Euler poles. They are analyzed in the bonus tutorial oriKernelRegression.R.



### WHAT HAVE WE LEARNED? ###

# Many kinds of orientational data can be converted into rotations. Each data type has its accompanying symmetry group.



### SEE ALSO ###

# The bonus tutorial oriEulerAngles.R plots these same data sets in terms of Euler angles.

# The bonus tutorial oriRodrigues.R plots these same data sets in the Rodrigues plot.

# The bonus tutorial oriKernelRegression.R analyzes pure rotations without symmetry.


