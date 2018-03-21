


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# This tutorial further illustrates the meaning of ellipsoid vectors by plotting special cases. We focus on normalized ellipsoids in 5D plots, rather than unnormalized ellipsoids in 6D plots.

# The data used here are all synthetic and generated on the spot.



### PRELIMINARY WORK ###

# If you have just restarted R, then set the working directory to the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc. Execute the following line of code to load our custom library and its dependencies.
source("library/all.R")



### SPHERE ###

# Here's how you make a normalized sphere. The rotation doesn't matter, so we choose it randomly.
synthEll <- ellEllipsoidFromRotationLogA(rotUniform(), c(0, 0, 0))
synthEll

# Notice that the vector is (0, 0, 0, 0, 0). Among normalized ellipsoids, there is one and only one sphere, and it plots at the origin of 5D space.
synthEll$vector
ellPairsPlot(list(synthEll$vector))



### AMBIENT ORIENTATION ###

# Here's how you make an ellipsoid aligned with the ambient coordinate axes. Notice that the first three coordinates of the vector are zero.
synthEll <- ellEllipsoidFromRotationLogA(diag(c(1, 1, 1)), rnorm(3), doNormalize=TRUE)
synthEll

# Let's make a bunch of them and plot them. Can you foresee what the plot will look like? Everything is in the v4-v5 plane.
synthElls <- replicate(100, ellEllipsoidFromRotationLogA(diag(c(1, 1, 1)), rnorm(3), doNormalize=TRUE), simplify=FALSE)
ellPairsPlot(lapply(synthElls, function(ell) ell$vector))



### OTHER FIXED ORIENTATION ###

# When all of our ellipsoids have the same orientation, but not the ambient one, then we get some other plane in the vector plots.
synthOri <- rotUniform()
synthElls <- replicate(100, ellEllipsoidFromRotationLogA(synthOri, rnorm(3), doNormalize=TRUE), simplify=FALSE)
ellPairsPlot(lapply(synthElls, function(ell) ell$vector))

# To see this plane we really have to look at 3D vector plots, not 2D ones. Try a few combinations of coordinates, not just c(1, 2, 3).
ellVectorPlot(c(1, 2, 3), lapply(synthElls, function(ell) ell$vector))



### SPHEROIDS OF FIXED ORIENTATION ###

# Make 100 increasingly oblate spheroids with the same orientation. (The same thing happens with prolate spheroids.)
synthOri <- rotUniform()
synthAs <- seq(from=1, to=100, by=1)
synthElls <- lapply(synthAs, function(a) ellEllipsoidFromRotationLogA(synthOri, c(a, a, 1), doNormalize=TRUE))

# Can you guess what the plots will look like?
ellPairsPlot(lapply(synthElls, function(ell) ell$vector))
ellVectorPlot(c(1, 2, 3), lapply(synthElls, function(ell) ell$vector))



### SPHEROIDS OF FIXED SHAPE ###

# Make 100 oblate spheroids of the same shape but random orientations.
synthOris <- rotUniform(1000)
synthA <- c(2, 2, 1)
synthElls <- lapply(synthOris, function(r) ellEllipsoidFromRotationLogA(r, synthA, doNormalize=TRUE))

# The plots are more interesting than our recent examples, because ellipsoid vectors depend non-linearly on orientation.
ellPairsPlot(lapply(synthElls, function(ell) ell$vector))

# In the 3D plot, there's a good chance of seeing a https://en.wikipedia.org/wiki/Roman_surface or a https://en.wikipedia.org/wiki/Cross-cap.
ellVectorPlot(c(1, 2, 3), lapply(synthElls, function(ell) ell$vector))

# Here we tweak the plot. After making the plot, we maximize the window. Then we call afterMaximizingWindow to take a screen shot and save it to the file 'figEllVectors.png'. One of the figures in the readme.pdf was made exactly like this. To learn more about saving plots, see the bonus tutorial publicationPlots.R.
ellVectorPlot(c(1, 2, 3), lapply(synthElls, function(ell) ell$vector),
              backgroundColor="white", colors=c("black"))
afterMaximizingWindow(leftName="figEllVectors.png", zoom=0.6)



### WHAT HAVE WE LEARNED? ###

# By examining some special cases, we can improve our intuition for these weird ellipsoid vector plots.


