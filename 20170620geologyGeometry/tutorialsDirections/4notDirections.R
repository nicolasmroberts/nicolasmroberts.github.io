


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# Many geological data types are not simply directions, but pairs (or triples) of directions. As this tutorial demonstrates, we should NOT treat these directions separately, because that would ignore the relationship within each pair.

# This tutorial uses synthetic data sets, both loaded from file and generated on the spot.



### PRELIMINARY WORK ###

# If you have just restarted R, then set the working directory to the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc. Execute the following line of code to load our custom library and its dependencies.
source("library/all.R")



### MEAN PLANE POLE VS. MEAN LINE DIRECTION ###

# With your mouse, select all five of the following lines of code at once. Copy them, paste them into the Console pane, and press Return to execute them. What they do is: Randomly generate a synthetic data set of 10 plane-line pairs, such as foliations-with-lineations or fold-axial-planes-with-hinge-lines. Compute the angle between the mean plane pole and the mean line direction. It 'should' be 90 degrees? In the Console pane, press the up-arrow key and then Return, to execute the entire code block again. Try again and again, to see how the angle varies.
rotations <- rotIsotropicFisher(rotUniform(), 3, 10)
poles <- lapply(rotations, function(r) r[1,])
directions <- lapply(rotations, function(r) r[2,])
lineEqualAreaPlotTwo(poles, directions)
acos(dot(lineProjectedMean(poles), lineProjectedMean(directions))) / degree

# This code block runs 1,000 such experiments, making a histogram of the angles that result. Typically there are many results more than 5 degrees away from perpendicular. Sometimes there are even some results 15 degrees away.
perpendicularMeanTrial <- function(concen, n) {
  rotations <- rotIsotropicFisher(rotUniform(), concen, n)
  poles <- lapply(rotations, function(r) r[1,])
  directions <- lapply(rotations, function(r) r[2,])
  acos(dot(lineProjectedMean(poles), lineProjectedMean(directions))) / degree
}
hist(replicate(1000, perpendicularMeanTrial(3, 10)))

# Computing the sample mean is a really basic operation, that might affect a variety of more sophisticated techniques. So we want our notion of sample mean to be well-behaved, with reliable mathematical properties. Computing the mean directions of planes and lines separately isn't going to be acceptable.



### SYNTHETIC OUTLIER EXAMPLE ###

# Here's a synthetic data set of 28 plane-line pairs, such as foliation-lineations or fold orientations. The plane poles are plotted as circles and the lines are plotted as squares. I claim that there is a dramatic outlier in this data set --- one plane-line pair that is quite different from the others. Can you see it in this plot? (No.)
outlierData <- geoDataFromFile("data/synthFoldsOutlier.csv")
outlierData
lineEqualAreaPlot(c(outlierData$pole, outlierData$direction),
                  shapes=c(replicate(nrow(outlierData), "c"), replicate(nrow(outlierData), "s")))

# The problem with that plot is that it divorces each plane from its associated line. So a lot of information is lost. Here is another plot, that shows each line on its plane. In theory, you should be able to see the outlier now. In practice, it's still difficult to see.
lineEqualAreaPlot(outlierData$direction, curves=lapply(outlierData$pole, rayGreatCircle, 72), shapes="s")

# Now we arbitrarily color some of the plane-line pairs blue, some of them green, and one of them red. Do you see why the red datum is an outlier? (We'll see it even more clearly in the next few sections.)
lineEqualAreaPlot(c(outlierData$pole, outlierData$direction), colors=as.character(outlierData$color),
                  shapes=c(replicate(nrow(outlierData), "c"), replicate(nrow(outlierData), "s")))



### WHAT HAVE WE LEARNED? ###

# If your data consist of multiple inter-dependent directions, then you should not treat them as if they were independent.


