


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# In this tutorial, we learn about the sample mean and dispersion (the average of a data set, and how to measure the spread of the data about that average). We also study inference about the mean (confidence regions and hypothesis tests for the population mean). We encounter the limitations of asymptotic methods and therefore try a non-asymptotic, bootstrap-based approach instead.

# For data we use dike directions (strikes and dips) measured by Titus et al. in the Troodos ophiolite, Cyprus.



### PRELIMINARY WORK ###

# If you've just restarted R, then you may need to set your working directory. Your current working directory is named next to the word 'Console' at the top of RStudio's Console pane. It should be the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc. If it isn't, then go to the Files pane, navigate to the geologyGeometry directory, click the 'More' menu, and choose 'Set As Working Directory'.

# Then click on the line of code below, and press RStudio's Run button to execute it.
source("library/all.R")

# Similarly, run these lines of code one at a time. What they do is: Load some dike data from a file. (Each dike strike-dip is converted into a unit vector describing the pole to the dike.) Let you inspect the data in the Console pane and the Plots pane. (In the Plots pane, you can click the Zoom button to show the image in a larger window. You can resize that window if it looks poorly shaped.)
site230 <- geoDataFromFile("data/cyprusDikesSite230.tsv")
site230
lineEqualAreaPlot(site230$pole)



### MEAN AND DISPERSION ###

# Compute the sample mean. Inspect it in words and pictures.
site230Mean <- lineProjectedMean(site230$pole)
geoStrikeDipDegFromCartesian(site230Mean)
geoTrendPlungeDegFromCartesian(site230Mean)
lineEqualAreaPlotTwo(site230$pole, list(site230Mean), "black", "gray")

# Compute the scatter of the data about the mean. You get three perpendicular lines. The first line is the mean. The second line defines the direction of greatest girdling from the mean. The third line is the pole to the girdle.
site230Scatter <- lineMeanScatter(site230$pole)
lineEqualAreaPlotTwo(site230$pole,
                     list(site230Scatter$vectors[,1], site230Scatter$vectors[,2], site230Scatter$vectors[,3]),
                     "black", "gray",
                     curves=list(rayGreatCircle(site230Scatter$vectors[,3])))

# The scatter also includes three numbers, which are all non-negative and sum to 1. Roughly, the first number quantifies how concentrated the data are (about the mean). The second number quantifies how girdled the data are (away from the mean). The third number quantifies how spread-out-over-the-sphere the data are (away from the girdle).
site230Scatter$values



### ASYMPTOTIC CONFIDENCE REGION ###

# The basic idea of statistics is that you have some population that you're trying to understand. In this example, it's the dikes at site 230 in the Troodos ophiolite, Cyprus. You measure a small subset of that population. You wish to extrapolate from your data set, to draw inferences about the population of a whole. One of the most basic kinds of inference is: What is the mean of the population? Your data let you guess at where the mean might be --- namely, at the mean of your data set --- but your data also let you evaluate the uncertainty in that guess. A confidence region summarizes that uncertainty.

# In elementary statistics, the normal distribution plays a pivotal role. In directional statistics, there are several analogues of the normal distribution. Examples include the Fisher and Kent distributions for rays and the Watson and Bingham distributions for lines.

# For our Cyprus dikes, we can ask for a 95% confidence region for the Bingham mean. The answer is reported as two angles, in radians, pointing from the mean toward the other two principal directions from lineMeanScatter. Unfortunately, the sofware balks on the angles for two reasons: The sample size is too small (smaller than 25) and the data are too tightly concentrated to use standard Bingham lookup tables (a certain number omega1 is less than 0.02). The same thing happens when you use the 'Axial Bingham Distribution' feature in Allmendinger's and Cardozo's Stereonet software (last I checked).
site230Bing <- lineBinghamInference095(site230$pole)
site230Bing

# This classic method for constructing a directional confidence region is 'asymptotic' in sample size. It works best for large data sets. For small data sets it works too poorly to be taken seriously. So we need another approach.



### BOOTSTRAP CONFIDENCE REGION ###

# Let's try a non-asymptotic approach to inference about the mean: bootstrapping. We resample the data set with replacement, to form a new data set of the same size as the old one, but with some lines repeated and others omitted. Intuitively, this new data set is like the original data set perturbed a bit, and its mean is like the original mean perturbed a bit. We repeat this process many times, to get a bunch of perturbed means, which quantify the uncertainty in the original mean.
site230Inf <- lineBootstrapInference(site230$pole, numBoots=10000, numPoints=50)
lineEqualAreaPlot(site230Inf$bootstraps)

# We take the middle 95% of those bootstrapped means, to form a 95% confidence region for the population mean. Usually it is centered close to, but not exactly at, the original sample mean.
lineEqualAreaPlot(list(site230Mean), curves=list(site230Inf$points))



### BOOTSTRAP HYPOTHESIS TEST ###

# Now suppose, for the sake of argument, that a published paper claimed that the dikes in this region were all horizontal. Then their poles must be vertical, with plunge 90 degrees. We can use the bootstrapping results to test this hypothesis. The following code produces a 'p-value', which is the probability of seeing data like the site 230 data (or data 'more extreme' than them) if the mean dike pole really was vertical. If p < 0.05 (say) then we reject the assumption that the mean dike pole is vertical.
site230Inf$pvalue(geoCartesianFromTrendPlungeDeg(c(0, 90)))

# Confidence regions are intimately related to hypothesis tests. The boundary of the 95% confidence region consists of exactly those hypotheses with p-value p == 0.05. The following code demonstrates this idea by computing the p-value at a bunch of boundary points.
sapply(site230Inf$points, site230Inf$pvalue)

# Similarly, the points outside the region correspond to hypotheses that are rejected (p < 0.05) and points inside the region correspond to hypotheses that are not rejected (p > 0.05). The following code demonstrates this idea by computing the p-value for each bootstrapped mean, and then computing quantiles for those p-values. Notice that 5% of the bootstrapped means produce p < 0.05, and hence 5% of the bootstrapped means lie outside the 95% confidence region, as we intended.
quantile(sapply(site230Inf$bootstraps, site230Inf$pvalue), probs=c(0.00, 0.05, 0.25, 0.50, 0.75, 1.00))



### WHAT HAVE WE LEARNED? ###

# To describe your data set, at a bare minimum you typically want to plot the data, compute the mean, and compute some measure of dispersion.

# To draw inferences from your data about the larger population that it represents is more difficult. Some methods don't work well for small data sets.


