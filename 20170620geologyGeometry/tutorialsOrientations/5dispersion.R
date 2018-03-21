


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# In elementary statistics, when you report the mean of a data set, you often also report some measure of its dispersion, such as the standard deviation. In this tutorial, we mention several analogous concepts for rotational/orientational data. We use one of them to implement a method called crystallographic vorticity axis analysis (Michels et al., 2015).

# This tutorial uses quartz crystallographic orientations (Strine and Wojtal, 2004; Michels et al., 2015) and some synthetic data sets, randomly generated on the spot.



### PRELIMINARY WORK ###

# If you have just restarted R, then set the working directory to the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc. Execute the following line of code to load our custom library and its dependencies.
source("library/all.R")



### BEST PRACTICE SOMEDAY? ###

# Someday, when statistical practices are better established in structural geology, you will not be able to publish rotational/orientational data without reporting some elementary descriptors such as mean and dispersion. To that end, here are a few ideas for measuring dispersion. I don't want to dwell on them now, because in 2016 the applications below are more interesting and valuable.

# Build a synthetic data set of 100 rotations, by sampling from the matrix Fisher distribution (one of several competing notions of 'normal distribution' on the space of rotations).
practiceM <- rotUniform()
practiceK <- diag(c(15, 7, 3))
practiceRots <- rotFisher(practiceM, practiceK, 100)
rotEqualVolumePlot(practiceRots)

# The Frechet variance is the thing that the Frechet mean minimizes. As a single number, it cannot capture anisotropy --- whether the data are dispersed differently in different directions.
rotMeanVariance(practiceRots)

# The eigensystem of the quaternionic scatter matrix reproduces the projected arithmetic mean, along with four measures of dispersion, the last three of which 'point' toward certain other points in the space of rotations. So these measures of dispersion account for anisotropy. (If you remember the scatter statistics from the line exercise 1oneDirection.R, this is the same thing, but 'one dimension up'.)
rotMeanScatter(practiceRots)

# Given a set of rotations, you can compute the 'maximum likelihood estimate' (MLE) of the matrix Fisher parameters M, K that might have produced that data set. The MLE of M is again the projected arithmetic mean. The MLE of K measures dispersion anisotropically. 
rotFisherMLE(practiceRots)

# When your rotations are tightly concentrated, you can approximate them as points in the tangent space at the mean. Then principal component analysis in that tangent space gives you yet another measure of anisotropic dispersion. Zoom in on this plot to see the principal components.
practicePCA <- rotLeftPrincipalComponentAnalysis(practiceRots, rotFrechetMean(practiceRots), numPoints=5)
practicePCA$magnitudes
practicePCA$directions
rotEqualAnglePlot(practiceRots, curves=practicePCA$curves, simplePoints=TRUE)

# For orientations (rotations with symmetry), only the Frechet variance makes theoretical sense.
oriMeanVariance(practiceRots, group=oriLineInPlaneGroup)

# If your orientations are tightly concentrated, then you could probably regard them as rotations without distorting the truth too much. And then you could use any of the methods above. Whatever you do, be specific when reporting the mean and dispersion of a data set in a paper. Mention the measure that you are using, with a reference, so that readers 30 years later know how to evaluate your claims.

# With all of that out of the way, let's talk about how you might actually use such measures of dispersion to discover something about rocks.



### APPLICATION: CRYSTALLOGRAPHIC VORTICITY AXIS ###

# Reddy and Buchan (2005) proposed that, as a rock deforms, its mineral grains generally rotate about the vorticity axis of the deformation. Therefore, in a thin section of a deformed rock, one should be able to infer the direction of vorticity from how the crystallographic axes are dispersed. Michels et al. (2015) implemented this idea in an objective, quantitative way, using rotation statistics.

# To see how, let's return to this data set of 761 alpha-quartz orientations (6-fold symmetry) from a single grain, obtained by EBSD. In an equal-area plot, the 010 (green) and 001 (blue) axes seem to be dispersed along small circles with a common pole that is nearly vertical. The 100 (red) axes are not as clearly dispersed along a small circle, but that could be because they are close to the pole.
michelsData <- geoDataFromFile("data/moine_one_grainABCxyz.tsv")
michels100s <- lapply(michelsData$rotation, function(r) r[1,])
michels010s <- lapply(michelsData$rotation, function(r) r[2,])
michels001s <- lapply(michelsData$rotation, function(r) r[3,])
lineEqualAreaPlotThree(michels100s, michels010s, michels001s)

# So the question is: How do you fit circles to all three sets of axes simultaneously? To someone familiar with rotation statistics, the obvious answer is to bind the axes together into orientations and fit a curve to those orientations. So, for each alpha-quartz orientation, we put its three axes into the rows of a matrix, which is then subject to trigonal trapezohedral symmetry.
oriEqualAnglePlot(michelsData$rotation, group=oriTrigonalTrapezohedralGroup, simplePoints=TRUE)

# The orientations are tightly concentrated, so it is reasonable to work with one symmetric copy and ignore the others. Compute the mean orientation and choose representative rotations near it. Then the principal components capture the main directions of dispersion. In particular, the first principal component produces a geodesic that best-fits the cloud of orientations. Zoom in on this plot to see the geodesics.
michelsMean <- oriFrechetMean(michelsData$rotation, group=oriTrigonalTrapezohedralGroup)
michelsRots <- oriNearestRepresentatives(michelsData$rotation, michelsMean, group=oriTrigonalTrapezohedralGroup)
michelsPCA <- rotLeftPrincipalComponentAnalysis(michelsRots, michelsMean, numPoints=5)
oriEqualAnglePlot(michelsRots, curves=michelsPCA$curves, group=oriTrigonalTrapezohedralGroup, simplePoints=TRUE)

# Now we dissect the first principal geodesic, along its rows, into three curves of directions. They are small circles about the first principal direction.
michelsFirst <- lapply(michelsPCA$curves[[1]], function(r) r[1,])
michelsSecond <- lapply(michelsPCA$curves[[1]], function(r) r[2,])
michelsThird <- lapply(michelsPCA$curves[[1]], function(r) r[3,])
n <- nrow(michelsData)
geoTrendPlungeDegFromCartesian(michelsPCA$directions[,1])
lineEqualAreaPlot(points=c(michels100s, michels010s, michels001s, list(michelsPCA$directions[,1])),
                  colors=c(replicate(n, "red"), replicate(n, "green"), replicate(n, "blue"), "black"),
                  curves=list(michelsFirst, michelsSecond, michelsThird))

# This grain is one of 74,724 grains in a single sample. Each grain produces a dispersion axis in this manner. See the figure in the notes. So the sample has 74,724 dispersion axes. They end up being widely dispersed, but the mean is at trend-plunge (228, 84) degrees. So Michels et al. (2015) took that direction as the overall vorticity direction of the deformation that deformed this sample. And this direction was similar to estimates of vorticity direction from earlier studies.

# Michels et al. (2015) analyzed two other samples, from the Gem Lake shear zone and the western Idaho shear zone. They found that their technique reproduced earlier estimates of vorticity axis well in all three cases.



### WHAT HAVE WE LEARNED? ###

# There are several ways to measure dispersion. Some capture anisotropy and some do not. Some require tightly concentrated data, while others do not. One of them was just what we needed to implement the crystallographic vorticity axis idea.



### SEE ALSO ###

# The bonus exercise oriRigidFabric.R uses orientation dispersion measures to study the development of fabric in rocks bearing rigid ellipsoidal clasts.


