


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# This tutorial analyzes dike directions from many stations dispersed over a large field area. There is no reason to assume that the direction is constant over such a wide area. On the contrary, we use regression to quantify geographic tendencies in the dike directions.

# For data we again use dike directions measured by Titus et al., and taken from earlier published maps, in the Troodos ophiolite, Cyprus.



### PRELIMINARY WORK ###

# If you have just restarted R, then set the working directory to the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc. Execute the following line of code to load our custom library and its dependencies.
source("library/all.R")

# Load our dike data from a file. Inspect them in words and pictures, including their geographic locations.
ourData <- geoDataFromFile("data/cyprus_inside_corner_field_data.tsv")
ourData
lineKambPlot(ourData$pole)
plot(x=(ourData$easting / 1000), y=(ourData$northing / 1000))

# Let's do the same with this other data file, from a published map of the same area. Try mimicking the code above, to plot the data.
mapData <- geoDataFromFile("data/cyprus_inside_corner_map_data.tsv")
# YOUR CODE GOES HERE!

# Assuming that they are comparable to each other (!), combine the two data sets. Plot the poles colored by northing, from red (south) to magenta (north). You see a vague rainbow, which suggests that there is a relationship between northing and dike direction.
allData <- rbind(ourData, mapData)
lineEqualAreaPlot(allData$pole, colors=hues(allData$northing / 1000))

# Also, in a new window, make a 3D plot showing northing versus equal-area. You can resize, zoom, and rotate the 3D view. This is the same information, visualized in a different way. Try to imagine the best-fit curve through the point cloud.
lineEqualAreaScalarPlot(allData$pole, scales(allData$northing / 1000))



### REGRESSION ###

# Perform regression to find a great circle arc that best fits the northing-pole relationship. Several technical problems can arise, so we have to inspect the results to verify that error == 0 and minEigenvalue > 0. The goodness of fit is measured by R^2 = 0.39.
allRegr <- lineRescaledGeodesicRegression(allData$northing / 1000, allData$pole, numPoints=50)
allRegr$error
allRegr$minEigenvalue
allRegr$rSquared

# Inspect the results in words. We find that for each km of northing there are 6 degrees of rotation about an axis with trend-plunge (084, 59) degrees.
geoTrendPlungeDegFromCartesian(lower(allRegr$rotation[,3]))
allRegr$a / degree

# Inspect the results in a picture. Yes, the line seems to travel roughly from the red part of the data set to the blue/magenta part.
lineEqualAreaPlot(allData$pole, colors=hues(allData$northing / 1000), curves=list(allRegr$points))

# Here we make the same plot again, capturing it in a file 'figDirRegr.pdf'. One of the figures in the readme.pdf was made exactly like this. To learn more about saving plots, see the bonus tutorial publicationPlots.R.
pdf("figDirRegr.pdf")
lineEqualAreaPlot(allData$pole, colors=hues(allData$northing / 1000), curves=list(allRegr$points))
dev.off()



### SIGNIFICANCE ###

# Is our regression result significant, or could it be an artifact of the random variability in the data? To address that question, we perform a permutation test, somewhat like the Wellner (1979) test from the previous section. Intuitively, we are wondering whether the relationship between northings and dike poles is meaningful. So we arbitrarily reassign the northings to the poles and perform the regression again. If we get a greater R^2, then the new relationship explains the data better than the original one. And the new relationship is meaningless, so the original relationship was meaningless. If we get a lesser R^2, then that's a little evidence that the original relationship had some meaning.

# Concretely, we compute many (say, 1,000) permutations, let p be the proportion in which R^2 exceeds the original R^2, and call the original result significant if p is small (less than 5%, say). Some of the regressions fail, so be sure to check how many R^2 values we actually got. When I did this with 1,000 permutations, I got p == 0 / 974 == 0. So the result is significant. For the sake of time, here is the test with only 100 permutations. If it's taking too much time and you want to stop it, press the Stop button in the upper right corner of RStudio's Console pane.
allRSquareds <- lineGeodesicRegressionPermutations(allData$northing / 1000, allData$pole, 100)
length(allRSquareds)
sum(allRSquareds > allRegr$rSquared) / length(allRSquareds)



### WHAT HAVE WE LEARNED? ###

# In geology, it is common for data to exhibit spatial or temporal patterns. Regression is one tool for quantifying such patterns. But you need to make sure that you're not over-detecting patterns that aren't really there.


