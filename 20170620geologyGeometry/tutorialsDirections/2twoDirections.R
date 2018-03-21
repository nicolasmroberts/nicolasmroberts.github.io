


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# This tutorial builds on the previous one. We perform a two-sample hypothesis test to determine whether two data sets are significantly different. Actually, we do two versions: Wellner (1979) and bootstrapping the mean. Then we discuss some of the geological implications.

# For data we again use dike directions measured by Titus et al. in the Troodos ophiolite, Cyprus. We also use dikes from earlier published maps of the area.



### PRELIMINARY WORK ###

# If you have just restarted R, then set the working directory to the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc. Execute the following line of code to load our custom library and its dependencies.
source("library/all.R")

# Load our dike data from two data files. Glue the data sets together. Inspect the data in words and pictures. Just for fun, superimpose Kamb contours (3-sigma, 6-sigma, 9-sigma, 12-sigma, by default).
site229 <- geoDataFromFile("data/cyprusDikesSite229.tsv")
site230 <- geoDataFromFile("data/cyprusDikesSite230.tsv")
ourDikes <- rbind(site229, site230)
ourDikes
lineKambPlot(ourDikes$pole)

# Load some more dike data from earlier published maps. Inspect it in words and pictures.
mapDikes <- geoDataFromFile("data/cyprusDikes_ic_three_blocks.tsv")
mapDikes
lineKambPlot(mapDikes$pole)

# Inspect both data sets together. Notice that the map dikes (cyan) seem a little steeper than our dikes (red), but the two data sets overlap each other pretty thoroughly.
lineEqualAreaPlotTwo(ourDikes$pole, mapDikes$pole)

# Compute their means and inspect those. Yes, the map dike mean is a little steeper than our dike mean. Of course, there is no reason to expect the two means to be exactly identical.
ourMean <- lineProjectedMean(ourDikes$pole)
mapMean <- lineProjectedMean(mapDikes$pole)
lineEqualAreaPlotTwo(list(ourMean), list(mapMean))

# We're inspecting these two data sets at the same time, because we're hoping to combine them into a single, larger data set. Are they compatible? Or is there some difference between them that renders them incomparable? To address this issue, we try two approaches, both based on simulation: the permutation test of Wellner (1979), and bootstrapping.



### WELLNER (1979) TEST ###

# Wellner (1979) defined a statistic, denoted T, that measures how different two sets of lines are. Let's compute it for these two data sets. I don't expect the number that comes out to give you any insight right now. But if you used this statistic a lot, on many data sets, then you might start to get some intuition for what it means. And its real purpose is to be used in the permutation test below.
lineWellner(ourDikes$pole, mapDikes$pole)

# We have 31 dikes of our own and 38 dikes taken from a map. We are trying to determine whether the distinction between 'our dikes' and 'map dikes' is a meaningful concept. So we randomly reassign the dikes to these two groups. That is, of the 69 total dikes, we arbitrarily pick 31 to call 'ours', and we call the other 38 'from the map'. Then we recompute T. Intuitively, if the new value of T is greater than the original value of T, then the new distinction between 'ours' and 'from the map' is more meaningful than the original distinction. But the new distinction is meaningless, so the original distinction was probably meaningless too. If the new value of T is less than the original value, then we have a little evidence that the original distinction was meaningful.

# So in practice we compute a large number (say, 10,000) permutations of the data, and count how many of them produce T greater than the original T. If fewer than 5% of them produce greater T, then the original distinction between 'ours' and 'from the map' is statistically significant at the 95% confidence level. 
lineWellnerInference(ourDikes$pole, mapDikes$pole, 10000)

# When I did Wellner's test, I got p = 0.0016 based on 10,000 permutations and p = 0.00193 based on 100,000 permutations. Assuming that your results are similar, we have strong evidence that our dikes and the map dikes come from different populations. We'll discuss what this result means later.



### BOOTSTRAPPING ###

# Let's explore the same question using another approach. For each of the two data sets, we bootstrap the mean, to get an idea of its uncertainty. (See the previous exercise for more explanation of bootstrapping.) If the two clouds of bootstrapped means didn't overlap at all, then we could be very confident that the means were different. But in this example the two clouds do overlap a bit.
ourInf <- lineBootstrapInference(ourDikes$pole, numBoots=10000, numPoints=50)
mapInf <- lineBootstrapInference(mapDikes$pole, numBoots=10000, numPoints=50)
lineEqualAreaPlotTwo(ourInf$bootstraps, mapInf$bootstraps)

# So we need to look at those bootstrapped means more closely. For each data set, we construct a 95% confidence ellipse that contains the middle 95% of those bootstrapped means. Because the two ellipses do not overlap, we can conclude that the two populations are different.
lineEqualAreaPlotTwo(
  list(ourInf$center), list(mapInf$center),
  curves=list(ourInf$points, mapInf$points))

# To illustrate p-values again, let's compute the greatest p-value attained by any pole in the map confidence region, according to the notion of p-value produced by our dikes. And vice-versa. Intuitively, each data set rejects the other.
max(sapply(mapInf$points, ourInf$pvalue))
max(sapply(ourInf$points, mapInf$pvalue))



### WHAT HAVE WE LEARNED? ###

# In this section, we've gained a little more experience working with plotting, description, and inference.

# The two dike populations look pretty similar, but statistics was able to detect a significant difference. Sometimes statistics is sharper than your naked eye.

# What does it mean, geologically, that the dikes from the two studies come from different populations? We were hoping to combine them into one larger data set. Is that still a good idea? I'd say no. Data from the two data sets are not comparable. Perhaps something was different in the methodology of the two studies. For example, maybe they sampled from different locations, and that difference turns out to matter. Or maybe the measurement apparatus (pocket transit compass with human operator) was different in the two studies.


