


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# Building on 7geodesicRegression.R in tutorialsOrientations, we learn a second regression technique, called kernel regression. Compared to geodesic regression, this technique is capable of fitting much more sophisticated curves. Unfortunately, the curves can be so complicated that they are difficult to summarize in words.

# The data consist of relative motions of the Farallon and Pacific plates (Engebretson et al., 1984; Prentice, 1987).



### PRELIMINARY WORK ###

# If you have just restarted R, then set the working directory to the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc. Execute the following line of code to load our custom library and its dependencies.
source("library/all.R")

# Load the Engebretson et al. (1984) data as cited by Prentice (1987). These are reconstructed relative plate motions between the Farallon and Pacific plates. Times are in millions of years before present.
source("data/farpacRotations.R")
rotEqualVolumePlot(farpacRs, colors=hues(farpacTs))



### GEODESIC REGRESSION ###

# Do geodesic regression. Check that error is 0 and minEigenvalue is positive. (If not, then something has gone wrong in the optimization. See the documentation for rotGeodesicRegression.)
farpacGeod <- rotGeodesicRegression(farpacTs, farpacRs, numSteps=1000)
farpacGeod$error
farpacGeod$minEigenvalue

# Check the quality of the fit, both in the R^2 statistic and graphically.
farpacGeod$rSquared
rotGeodesicRegressionPlot(farpacTs, farpacRs, farpacGeod)

# A 1,000-permutation test yielded p = 0 / 1000 = 0. So the result seems significant. And R^2 = 0.55 is quite respectable. We have captured quite a tendency in the data.

# On the other hand, the data are clearly following some curve other than the one that we just found. That curve is apparently just too subtle to be described by a geodesic (or else we would have found it). So we shouldn't be satisfied with this result. We need to consider more sophisticated curves.



### KERNEL REGRESSION ###

# Kernel regression is another kind of regression. Instead of producing parameters that describe a smooth curve, it instead produces a finite sequence of points along the curve. What you need to know right now is that the algorithm depends on a number h, called the 'bandwidth', that affects the smoothing of the curve.

# The following commented-out line takes a long time to execute. What it does is select the bandwidth for you automatically, using the algorithm of Davis et al. (2010). I got 0.03976778. So let's go with that.
# farpacH <- rotBandwidthForKernelRegression(farpacTs, farpacRs, dnorm)
farpacH <- 0.03976778

# Do the kernel regression. This might take a minute. Check that errors are 0 and minEigenvalues are positive again.
farpacKerns <- lapply(
  seq(from=min(farpacTs), to=max(farpacTs), by=1),
  rotKernelRegression, farpacTs, farpacRs, farpacH, numSteps=1000)
sapply(farpacKerns, function(regr) regr$error)
sapply(farpacKerns, function(regr) {regr$minEigenvalue > 0})

# Plot the regression curve.
farpacKernCurve <- lapply(farpacKerns, function(regr) regr$r)
rotEqualVolumePlot(farpacRs, curves=list(farpacKernCurve), colors=hues(farpacTs))

# We get a really tight fit, with R^2 = 0.98.
farpacRsquared <- rotRsquaredForKernelRegression(farpacTs, farpacRs, farpacH, numSteps=1000)
farpacRsquared

# Do a permutation test for significance. For the sake of time, we do only 10 permutations, although that is far too few to tell us anything. When I did this with 1,000 permutations, I got p = 0. So the result is significant.
rSquareds <- rotKernelRegressionPermutations(farpacTs, farpacRs, farpacH, numPerms=10)
length(rSquareds)
sum(rSquareds > farpacRsquared$rSquared)

# Here is how you might use such a regression: Predict the relative motion of the Farallon and Pacific plates 100 million years ago. The relative rotation was about axis with trend 350 degrees and plunge 51 degrees, *clockwise* through 84 degrees.
farpac100 <- rotKernelRegression(100, farpacTs, farpacRs, farpacH, numSteps=1000)
farpac100
geoTrendPlungeAngleDegFromMatrix(farpac100$r)



### WHAT HAVE WE LEARNED? ###

# Kernel regression can detect a different kind of dependency than can geodesic regression.


