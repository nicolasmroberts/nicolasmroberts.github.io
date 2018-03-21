


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# This exercise demonstrates various options for tweaking equal-volume and angle-axis plots. You can use these options to adjust color and weighting in preparation for publication. We describe how to capture the plots as raster images. We also discuss tweaking and capturing some 2D plots.

# We don't use any natural data --- just synthetic data sets of rotations and ellipsoids.



### PRELIMINARY WORK ###

# If you have just restarted R, then set the working directory to the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc. Execute the following line of code to load our custom library and its dependencies.
source("library/all.R")

# Here is the synthetic example from orientation tutorial 7geodesicRegression.R.
source("data/synthRotationsRegression.R")
synthGeod <- rotGeodesicRegression(synthXs, synthRs)
rotGeodesicRegressionPlot(synthXs, synthRs, synthGeod)



### COLOR ###

# For publication, you probably want a white background. So you can't have the curve be white too. Let's make the curve black.
rotGeodesicRegressionPlot(synthXs, synthRs, synthGeod,
                          backgroundColor="white", curveColor="black")

# Printers can't handle subtle differences in shade reliably --- or at least journal publishers don't make any promises. So you may want to turn off the boundary sphere by setting its alpha (opacity) to 0.
rotGeodesicRegressionPlot(synthXs, synthRs, synthGeod,
                          backgroundColor="white", curveColor="black",
                          boundaryAlpha=0)

# Unless this is a special occasion, you want grayscale instead of color. So let's color the data by shades instead of hues. You don't want the shades to go all the way from black (0) to white (1). So instead make them go from 0 to 0.75. Also let's make all of the axes black. (You're going to have to keep track of which axis is which in some other way.)
rotGeodesicRegressionPlot(synthXs, synthRs, synthGeod,
                          backgroundColor="white", curveColor="black",
                          boundaryAlpha=0,
                          colors=shades(synthXs, c(0, 0.75)), axesColors=c("black", "black", "black"))

# Now we run into another problem. By default, these plots have a fog effect, that makes distant points appear dimmer than points near the viewpoint. That's a useful depth cue, and you might want to keep it in publication. But it interferes with the shading that we just introduced. So, if you want to turn off the fog effect, do this.
rotGeodesicRegressionPlot(synthXs, synthRs, synthGeod,
                          backgroundColor="white", curveColor="black",
                          boundaryAlpha=0,
                          colors=shades(synthXs, c(0, 0.75)), axesColors=c("black", "black", "black"),
                          fogStyle="none")



### WEIGHT ###

# The regression curve might be too thin. It might not reproduce well on paper. So let's beef it up. (The axes have width 3, and this is not currently customizable.)
rotGeodesicRegressionPlot(synthXs, synthRs, synthGeod,
                          backgroundColor="white", curveColor="black",
                          boundaryAlpha=0,
                          colors=shades(synthXs, c(0, 0.75)), axesColors=c("black", "black", "black"),
                          fogStyle="none",
                          curveWidth=5)

# The points look good to me, but if you want to adjust their weight, there are two ways to go about it. One way is point size. Here's an extreme example.
rotGeodesicRegressionPlot(synthXs, synthRs, synthGeod,
                          backgroundColor="white", curveColor="black",
                          boundaryAlpha=0,
                          colors=shades(synthXs, c(0, 0.75)), axesColors=c("black", "black", "black"),
                          fogStyle="none",
                          curveWidth=5,
                          pointSize=0.05)

# The other way to adjust point weight is to make the points 'simple'. Simple points have the advantage of being faster to render, so we often use them when there are too many data for non-simple points. Their point size is measured in pixels, and zooming does not affect their apparent size. Depending on your particular case, that can be desirable or undesirable.
rotGeodesicRegressionPlot(synthXs, synthRs, synthGeod,
                          backgroundColor="white", curveColor="black",
                          boundaryAlpha=0,
                          colors=shades(synthXs, c(0, 0.75)), axesColors=c("black", "black", "black"),
                          fogStyle="none",
                          curveWidth=5,
                          simplePoints=TRUE, pointSize=10)



### OUTPUT ###

# Once you have your colors and weights configured exactly as you like, maximize your window and take a screen shot. You can subsequently edit the screen shot to crop out junk, add axis labels, etc. To get publication-quality plots, you probably want a large screen.

# If you want to capture your orientation plots from specific viewpoints, so that they are comparable to each other, then the following feature may help. You make your plot, then maximize the window, then call the afterMaximizingWindow function with two file names. (Currently only PNG files are supported maybe.) It chooses a certain viewpoint and zoom level, and captures the window in the first file. Then it chooses another viewpoint, re-zooms, and captures the window in the second file.
rotGeodesicRegressionPlot(synthXs, synthRs, synthGeod,
                          backgroundColor="white", curveColor="black",
                          boundaryAlpha=0,
                          colors=shades(synthXs, c(0, 0.75)), axesColors=c("black", "black", "black"),
                          fogStyle="none",
                          curveWidth=5)
# Maximize the plot window before calling the next line!
afterMaximizingWindow("myExampleLeft.png", "myExampleRight.png")



### ELLIPSOIDS ###

# Following bonus tutorial ellVectors.R, make a synthetic data set of spheroids with equal shapes and uniformly random directions.
synthOris <- rotUniform(1000)
synthA <- c(2, 2, 1)
synthVecs <- lapply(synthOris, function(r) ellEllipsoidFromRotationLogA(r, synthA, doNormalize=TRUE)$vector)

# Here are some vector plots.
ellPairsPlot(synthVecs)
ellVectorPlot(c(1, 2, 5), synthVecs)

# Here is how to capture the 2D plot as a PDF.
pdf("myEllipsoidPairs.pdf")
ellPairsPlot(synthVecs)
dev.off()

# The 3D plot is processed in much the same manner as the rotational example above. We don't need to adjust boundaryAlpha, because there is no boundary sphere. But we do adjust the zoom, so that no points are accidentally cropped out.
ellVectorPlot(c(1, 2, 5), synthVecs,
              colors=c("black"), backgroundColor="white", axesColors=c("black", "black", "black"))
afterMaximizingWindow("myEllipsoidsLeft.png", "myEllipsoidsRight.png", zoom=0.6)



### WHAT HAVE WE LEARNED? ###

# The graphics routines in this library may not have every customization you want (feel free to suggest new ones), but they do give you some options.


