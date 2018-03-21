


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# In this tutorial we have a simple goal: Compute the average of a set of rotations. Along the way, we learn the concept of the distance between two rotations, and we encounter a situation where we shouldn't average our data.

# Data used in this tutorial include the foliation-lineation pairs of Giorgis and Tikoff (2004) and a couple of synthetic data sets.



### PRELIMINARY WORK ###

# If you have just restarted R, then set the working directory to the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc. Execute the following line of code to load our custom library and its dependencies.
source("library/all.R")



### SOME IDEAS, MOSTLY BAD ###

# When you average angles in 2D, you know to be careful around the 360-degree mark. For example, the average of strikes 003 degrees and 359 degrees is 001 degrees, even though the numerical average of 3 and 359 is 181. Well, these problems only get worse in 3D. Here is a partial list of bad ideas for computing the average of a set of rotations:

# * Express the rotations as three Euler angles each, and average those. In addition to the 360-degree-mark problem just mentioned, this idea suffers from the fact that the rotations don't commute with each other. The results can be terrible.

# * Express the rotations as angles about axes. Average the angles as you would any scalar. Average the axes using rayProjectedMean or lineProjectedMean, as in earlier sections. But divorcing the angles from the axes like this causes lots of bad behavior.

# * Express the rotations as 'infinitesimal rotations', which inhabit a vector space. Average those, and then convert back to finite rotations. This approach suffers from problems near the 180-degree mark. But it forms the basis for computing the Frechet mean below.

# * Express the rotations as matrices, and average those. The average matrix is usually not a rotation at all. But this idea forms the basis for the projected arithmetic mean below.



### DISTANCE BETWEEN TWO ROTATIONS ###

# The 'difference' between two rotations R and Q is a rotation. Namely, if we rotate an object by R, and then rotate it by Q R^T, then the overall rotation is Q R^T R = Q. So Q R^T is the difference between R and Q. The distance between R and Q, as two points in the space of rotations, is simply the angle of that rotation Q R^T. It is also the length of the shortest 'geodesic' curve connecting R to Q. Here is a random example, the distance computed by hand and by a built-in function, and an illustrative plot.
synthR <- rotUniform()
synthQ <- rotUniform()
synthUA <- rotAxisAngleFromMatrix(synthQ %*% t(synthR))
synthUA[[4]] / degree
rotDistance(synthR, synthQ) / degree
synthCurve <- rotGeodesicPoints(synthR, synthQ, 100)
rotEqualAnglePlot(points=list(synthR, synthQ), curves=list(synthCurve))

# Alternatively, we can dissect each of R and Q into three unit vectors along its columns. The rotation Q R^T is the smallest rotation that rotates the first column of R to the first column of Q, and similarly for the second and third. The distance between R and Q is the angle of that smallest rotation.
synthFirsts <- lapply(synthCurve, function(r) r[,1])
synthSeconds <- lapply(synthCurve, function(r) r[,2])
synthThirds <- lapply(synthCurve, function(r) r[,3])
rayEqualAreaPlotThree(list(synthR[,1], synthQ[,1]), list(synthR[,2], synthQ[,2]), list(synthR[,3], synthQ[,3]),
                      curves=list(synthFirsts, synthSeconds, synthThirds))

# A question for you to ponder: What is the relationship between the following vector and that equal-area plot?
geoTrendPlungeDegFromCartesian(lower(synthUA[1:3]))



### FRECHET MEAN ###

# Load the foliation-lineation data from the western Idaho shear zone. Choose representative rotations to be close to the first rotation. Compute the Frechet mean, which is one of two popular notions of mean rotation.
wiszData <- geoDataFromFile("data/wiszFollins.tsv")
wiszRots <- oriNearestRepresentatives(wiszData$rotation, group=oriLineInPlaneGroup)
wiszFrechet <- rotFrechetMean(wiszRots)
wiszFrechet

# Here's your intuition for the Frechet mean. At any point R in the space of rotations, the 'variance' of the data set is a measure of how spread about the data are around R. It is proportional to the sum of the squared distances from R to the data. The Frechet mean is the point R that minimizes this variance. In the following plot, the Frechet mean is at the center of the 'spider', and the variance is proportional to the sum of the squared lengths of the spider legs.
wiszCurves <- lapply(wiszRots, function(r) rotGeodesicPoints(wiszFrechet, r, 10))
rotEqualAnglePlot(points=wiszRots, curves=wiszCurves)



### PROJECTED ARITHMETIC MEAN ###

# The projected arithmetic mean is defined in an entirely different way. First we take the average of the rotations, as 3x3 matrices.
wiszArith <- arithmeticMean(wiszRots)
wiszArith

# Here's a test of how well you remember tutorial 1linesInPlanes.R: What two mathematical properties must you check, to see whether that matrix is a rotation? Write code to check them here, if you like.

# No, wiszArith is not a rotation matrix. But it's nearly a rotation, and we can find the rotation closest to it. That's the projected arithmetic mean.
wiszProjected <- rotProjectedMean(wiszRots)
wiszProjected

# The projected arithmetic mean is often close to, but not identical to, the Frechet mean.
rotDistance(wiszFrechet, wiszProjected) / degree
rotEqualAnglePlot(c(wiszRots, list(wiszFrechet, wiszProjected)),
                  colors=c(replicate(length(wiszRots), "white"), "yellow", "green"))



### MULTIMODALITY ###

# Load a synthetic data set of 105 current ripple orientations (trivial symmetry) from a file. Compare the two means.
rippleData <- geoDataFromFile("data/synthCurrentRipples.csv")
rippleFrechet <- rotFrechetMean(rippleData$rotation)
rippleProjected <- rotProjectedMean(rippleData$rotation)
rotDistance(rippleFrechet, rippleProjected) / degree

# Shockingly, they're 27 degrees apart. Maybe a plot will help us understand what's going on.
rotEqualVolumePlot(c(rippleData$rotation, list(rippleFrechet, rippleProjected)),
                   colors=c(replicate(length(rippleData$rotation), "white"), "yellow", "green"))

# The data set is multimodal, meaning that it contains multiple clumps of high density. So there is no single 'central tendency' in the data. The very concept of 'mean' is meaningless. We should have checked unimodality before even computing the mean.

# There are techniques for detecting and dissecting the clumps. Some of them are described in the bonus tutorial clustering.R.



### WHAT ABOUT SYMMETRY? ###

# When your data are not pure rotations but rather rotations subject to symmetry, you have to compute the mean in a way that pays attention to the symmetry. In this example, the careless answer is 87 degrees away from the right answer.
wiszData <- geoDataFromFile("data/wiszFollins.tsv")
wiszRotMean <- rotFrechetMean(wiszData$rotation)
wiszOriMean <- oriFrechetMean(wiszData$rotation, group=oriLineInPlaneGroup)
oriDistance(wiszRotMean, wiszOriMean, group=oriLineInPlaneGroup) / degree

# The careless answer minimizes the sum of the squared lengths of the legs of this spider (shown four times due to symmetry). But this spider has chosen rotations to represent the orientations poorly. These legs do not depict the true distances between orientations.
wiszRotCurves <- lapply(wiszData$rotation, function(r) rotGeodesicPoints(wiszRotMean, r, 10))
oriEqualAnglePlot(wiszData$rotation, curves=wiszRotCurves, group=oriLineInPlaneGroup)

# The right answer minimizes the sum of the squared lengths of the legs of this other spider.
wiszOriCurves <- lapply(
  oriNearestRepresentatives(wiszData$rotation, wiszOriMean, group=oriLineInPlaneGroup),
  function(r) rotGeodesicPoints(wiszOriMean, r, 10))
oriEqualAnglePlot(wiszData$rotation, curves=wiszOriCurves, group=oriLineInPlaneGroup)

# By the way, you can also use the projected arithmetic mean on orientations. In this example, they differ by half a degree.
wiszOriProjMean <- oriProjectedMean(wiszData$rotation, group=oriLineInPlaneGroup)
wiszOriProjMean
oriDistance(wiszOriMean, wiszOriProjMean, group=oriLineInPlaneGroup) / degree

# In this example, the four symmetric copies of the data are so clearly separated, that one could choose representatives cleverly and then use a pure rotation mean. In fact, we did that on line 62 above. But other data sets are not so amenable.



### SO WHICH ONE DO I USE? ###

# The Frechet mean has some nice mathematical properties, and the projected arithmetic mean has some other nice mathematical properties. There is no objective way to choose one over the other. This multiple-ways-to-do-it kind of situation seems to arise throughout statistics. A prudent approach would be to report which one you are using, and how different it is from the other one.



### WHAT HAVE WE LEARNED? ###

# There are two notions of mean orientation, which behave pretty similarly, if you should be taking the mean at all.


