


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# In this tutorial we perform regression of ellipsoids, based on their locations in a field area. The regression is not very successful, and we don't find a strong pattern, but it's a valid demonstration.

# The data in question here are shape preferred orientation (SPO) ellipsoids of orthopyroxene clasts from New Caledonia (Titus et al.), constructed using the method of Robin (2002).



### PRELIMINARY WORK ###

# If you have just restarted R, then set the working directory to the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc. Execute the following line of code to load our custom library and its dependencies.
source("library/all.R")

# Load a data set of 34 shape preferred orientation (SPO) ellipsoids from New Caledonia. We will focus on the SPOs derived from field measurements of orthopyroxene clasts. Inspect the first datum of the 34 data, just to see what it's like. It's a station name, an easting-northing in km relative to a certain point in the high-strain zone, and an orthopyroxene SPO ellipsoid.
source("data/newcalOPXSpinelSPO.R")
length(ncOrthos)
ncOrthos[[1]]



### EXPLORING POSSIBLE TENDENCIES ###

# Here is a map view of the stations. They lie along a transect that cuts across an oceanic transform fault. The fault is roughly NS-striking. The eastings are relative to the fault, so the fault is at easting 0 in this plot.
plot(t(sapply(ncOrthos, function(station) station$en)))

# Color the stations by easting and plot the normalized ellipsoid vectors.
ncColors <- hues(sapply(ncOrthos, function(station) station$en[[1]]))
ellPairsPlot(lapply(ncOrthos, function(station) station$ortho$vector), colors=ncColors)

# There's too much yellow-green, because most of the stations are clumped near the transform. So let's distort easting in a way that stretches out those stations a bit.
ncColors <- hues(sapply(ncOrthos, function(station) cubeRoot(station$en[[1]])))
ellPairsPlot(lapply(ncOrthos, function(station) station$ortho$vector), colors=ncColors)

# In many of the panels, the reds and purples are near each other and the greens are elsewhere. So the tendency starts in one place, goes to another place in the fault, and returns to the first place. Perhaps distance from the fault, rather than easting itself, is relevant. So we replace easting with an even function of easting.
ncColors <- hues(sapply(ncOrthos, function(station) abs(cubeRoot(station$en[[1]]))))
ellPairsPlot(lapply(ncOrthos, function(station) station$ortho$vector), colors=ncColors)

# We're starting to see some clearer rainbows, especially in the fifth component. So let's look at some other plots.
ellHsuNadaiPlot(lapply(ncOrthos, function(station) station$ortho$logA), es=0.4, colors=ncColors)
ellEqualVolumePlot(lapply(ncOrthos, function(station) station$ortho$rotation),
                   lapply(ncOrthos, function(station) station$ortho$a), colors=ncColors)
ellEqualAreaPlot(lapply(ncOrthos, function(station) station$ortho$rotation),
                 lapply(ncOrthos, function(station) station$ortho$a), colors=ncColors)



### REGRESSION ###

# Here we do regression of ellipsoid vectors V versus easting e. Namely, we fit a curve of the form
#     V  ==  V0 + V1 e + V2 e^2 + V3 e^3 + V4 e^4,
# where V0, V1, V2, V3, V4 are all 5D vectors giving the coefficients of the curve.

# Ordinarily I try to package R's ugliness so that you don't have to see it. But I haven't figured out how to package this kind of thing. So we're going to do it explicitly.
ncVectors <- t(sapply(ncOrthos, function(station) station$ortho$vector))
ncEastings <- sapply(ncOrthos, function(station) station$en[[1]])
ncRegr <- lm(ncVectors ~ 1 + ncEastings + I(ncEastings^2) + I(ncEastings^3) + I(ncEastings^4))

# Here are the best-fit parameters: The first row is V0, the next row is V1, etc.
coefficients(ncRegr)

# Here is a detailed report of the regression results. Fundamentally, fitting a 5D curve amounts to fitting five 1D curves, one to each coordinate independently. So the summary consists of five similar-looking parts. The fourth and fifth parts have small p-values. So we have strong evidence of a correlation between easting and the fourth component of the ellipsoid vector, and we have strong evidence of a correlation between easting and the fifth component too.
summary(ncRegr)

# Plot the components with the smallest p-values. The fitted curve is kind of crazy.
ncPredCurve <- predict(ncRegr, 
                       data.frame(ncEastings=seq(from=min(ncEastings), to=max(ncEastings), length.out=100)))
ncPredCurve <- lapply(1:nrow(ncPredCurve), function(i) ncPredCurve[i,])
ncColors <- hues(sapply(ncOrthos, function(station) station$en[[1]]))
ellVectorPlot(c(4, 5), lapply(ncOrthos, function(station) station$ortho$vector),
              curves=list(ncPredCurve), colors=ncColors)
ellVectorPlot(c(2, 4, 5), lapply(ncOrthos, function(station) station$ortho$vector),
              curves=list(ncPredCurve), colors=ncColors)

# Plot ellipsoid orientation and shape. The fitted curve looks even crazier now.
ncPredEllCurve <- lapply(ncPredCurve, function(v) ellEllipsoidFromVector(v))
ellEqualVolumePlot(lapply(ncOrthos, function(station) station$ortho$rotation),
                   lapply(ncOrthos, function(station) station$ortho$a),
                   list(lapply(ncPredEllCurve, function(ell) ell$rotation)),
                   list(lapply(ncPredEllCurve, function(ell) ell$a)),
                   colors=ncColors)
ellHsuNadaiPlot(lapply(ncOrthos, function(station) station$ortho$logA),
                curves=list(lapply(ncPredEllCurve, function(ell) ell$logA)),
                colors=ncColors, es=0.4)



### USING THE RESULT ###

# Recall that there are no stations between 20 km and 60 km easting.
plot(t(sapply(ncOrthos, function(station) station$en)))

# The regression can be used to predict such missing values.
nc40Vector <- as.numeric(predict(ncRegr, data.frame(ncEastings=40)))
ellEllipsoidFromVector(nc40Vector)

# We'd like to be able to verbally summarize the result of the regression --- something along the lines of, 'for every km of easting, such-and-such happens'. But this is difficult, because we fitted a non-linear curve to the data. If we'd fitted a linear curve, it would have a simple 'for every km of easting...' meaning. But a linear curve would never have fit the data very well.

# We believe that we've found the pattern in the fourth and fifth components, but not in the other components. So let's replace the fourth and fifth components by their predictions from the regression. Intuitively, we're smoothing out the fourth and fifth components, leaving the others untouched.
ncVectors <- lapply(ncOrthos, function(station) station$ortho$vector)
ncPreds <- predict(ncRegr, data.frame(ncEastings=ncEastings))
ncPreds <- lapply(1:nrow(ncPreds), function(i) ncPreds[i,])
ncSmoothed <- lapply(1:length(ncEastings), function(i) ellEllipsoidFromVector(c(ncVectors[[i]][1:3], ncPreds[[i]][4:5])))

# Here are various plots of the smoothed data. Notice that the green points are crammed together in the fourth and fifth components. That's because they have similar eastings, and we have smoothed the data based on easting.
ncColors <- hues(sapply(ncOrthos, function(station) station$en[[1]]))
ellPairsPlot(lapply(ncSmoothed, function(ell) ell$vector), colors=ncColors)
ellHsuNadaiPlot(lapply(ncSmoothed, function(ell) ell$logA), es=0.4, colors=ncColors)
ellEqualVolumePlot(lapply(ncSmoothed, function(ell) ell$rotation),
                   lapply(ncSmoothed, function(ell) ell$a), colors=ncColors)
ellEqualAreaPlot(lapply(ncSmoothed, function(ell) ell$rotation),
                 lapply(ncSmoothed, function(ell) ell$a), colors=ncColors)

# Finally, if the regression had been highly successful --- significant results and low uncertainty in each component, and a good fit to the data --- then it might have steered us toward a fundamental law governing the development of SPO near oceanic transform faults. But this regression does not appear highly successful to me, basically because the data are too messy, probably because there is no simple fundamental law governing them.



### PROBLEMS ###

# When we predict the missing value at easting 40 km, we're extrapolating pretty far from the existing data. The behavior of the regression curve east of about 20 km easting is highly influenced by four very distant data. So you should take this prediction with a grain of salt.

# Actually, the regression results include confidence intervals on the coefficients, which you can use to understand the uncertainty in the fit. But remember that we found statistically significant correlations in only the fourth and fifth components of the ellipsoid vector. So the confidence intervals are going to be pretty big, especially in the first three components. With data this messy, you can't have high confidence and high precision at the same time.

# It would be nice to isolate the fourth and fifth components and re-interpret them in terms of ellipsoid size-shape and orientation. But a normalized ellipsoid's five degrees of freedom are all mixed up in the five components of its ellipsoid vector. Geometrically interpreting two components alone is typically quite difficult.

# Perhaps we should have fit only a degree-2 polynomial, rather than a degree-4 polynomial. The extra coefficients mostly turned out to be insignificantly different from zero. More visual inspection of the data might give us a better sense for what kind of curve to try. It would also help us understand how the variability in the data depend on easting. Regression works best when the variability is fairly constant in easting. But all of these issues are hard to evaluate when the data are five-dimensional.

# Finally, our choice of regression curve would ideally be motivated by some underlying physical theory of how SPO should form near an oceanic transform fault. When we lack a theory guiding the choice of curve, regression can sometimes help us discover it. But this regression certainly did not.



### WHAT HAVE WE LEARNED? ###

# When ellipsoids are rendered into their ellipsoid vectors, you can use standard multivariate regression techniques on them, to (try to) pick out patterns among them.



### SEE ALSO ###

# By changing 'ncOrthos' to 'ncSpinels' and '$ortho' to '$spinel' in the code above, you can do the same exercise for the spinel SPO ellipsoids. The results are similarly inconclusive.

# In tutorialsBonus, ellInferenceTwo.R does a paired two-sample test of whether the New Caledonia orthopyroxene and spinel SPOs agree.


