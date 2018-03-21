


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# Anisotropy of magnetic susceptibility (AMS) is a technique for measuring the magnetic fabric of a rock sample. X-ray computed tomography (XRCT) is a technique for imaging individual clasts inside a sample based on X-ray opacity contrast. Given a rock sample, it is natural to ask whether the AMS ellipsoid is simply the average of the XRCT ellipsoids. Or maybe it is the average of the highly magnetic XRCT ellipsoids? That would be nice to know, because AMS is easier/cheaper/commoner than XRCT. To some extent these questions can be answered using inference on ellipsoids.

# This exercise works through an unpublished data set of AMS and XRCT ellipsoids. The same analysis of other, similar data sets can be found in ellInferenceMore.R in tutorialsBonus.



### PRELIMINARY WORK ###

# If you have just restarted R, then set the working directory to the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc. Execute the following line of code to load our custom library and its dependencies.
source("library/all.R")

# Load one AMS ellipsoid from a file. The AGICO software normalizes the ellipsoids in a strange (to me) way. So we denormalize them and renormalize them. The meanSuscepts argument is needed for the denormalization. You get it from a separate AGICO output file. When there are multiple ellipsoids in the file, there should be one mean susceptibility for each ellipsoid.
ams205 <- geoEllipsoidDataFromAGICOFile("data/HU150205ams.tsv", meanSuscepts=c(8.957E-03), doNormalize=TRUE)
ams205 <- listFromDataFrame(ams205)[[1]]

# Load 4,001 XRCT ellipsoids from the same sample.
xrct205 <- geoEllipsoidDataFromAvizoFile("data/HU150205mafic.tsv", doNormalize=TRUE)
nrow(xrct205)

# Inspect the shapes, with AMS superimposed in green. The AMS ellipsoid is much closer to spherical than the XRCT ellipsoids. Is that consistent with the AMS being the mean of the XRCT?
ellHsuNadaiPlot(c(xrct205$logA, list(ams205$logA)),
                colors=c(replicate(length(xrct205$logA), "black"), "green"))
quantile(sapply(xrct205$logA, ellOctahedralShearStrain))

# Inspect the orientations. They look fairly uniform, although there is a little clumping. You can detect the clumping with Kamb contouring if you want, but it takes a few minutes to get a good image.
ellEqualVolumePlot(xrct205$rotation, xrct205$logA, simplePoints=TRUE)
ellEqualAreaPlot(xrct205$rotation, xrct205$logA)
#lineKambPlot(lapply(xrct205$rotation, function(r) r[1,]), numNonAdapt=5)
#lineKambPlot(lapply(xrct205$rotation, function(r) r[2,]), numNonAdapt=5)
#lineKambPlot(lapply(xrct205$rotation, function(r) r[3,]), numNonAdapt=5)

# Inspect the ellipsoid vectors, with the AMS superimposed in green. It seems plausible that the AMS ellipsoid could be the mean of the XRCT ellipsoids. But keep in mind that regions of high density might not be recognizable in this plot.
ellPairsPlot(c(xrct205$vector, list(ams205$vector)),
             colors=c(replicate(length(xrct205$vector), "black"), "green"))



### HYPOTHESIS TESTS ###

# Now we test whether the AMS ellipsoid is the mean of the XRCT ellipsoids, using four methods. The first is a standard asymptotic multivariate statistics test provided by the R package ICSNP. The next two are apparently sophisticated bootstrap methods provided by the R package FRB. The last one is my own custom bootstrap method, based on percentiles of Mahalanobis distance.
xrctInfT2 <- ellHotellingT2Inference(xrct205$vector, ams205$vector)
xrctInfMM <- ellBootstrapMMInference(xrct205$vector, ams205$vector, numBoots=1000)
xrctInfS <- ellBootstrapSInference(xrct205$vector, ams205$vector, numBoots=1000)
xrctInf <- ellBootstrapInference(xrct205$vector, numBoots=1000)

# Here are the p-values for the null hypothesis that the difference is zero. In plain English, what do they mean about these XRCT and AMS data?
xrctInfT2$p.value
xrctInfMM$p.value
xrctInfS$p.value
xrctInf$pvalue(ams205$vector)



### CONFIDENCE REGIONS ###

# The three bootstrap methods also give information about a confidence region. Here are component-wise 95% confidence intervals from the MM-type bootstrap inference.
xrctInfMM$CI

# Those five confidence intervals define a 5-dimensional box with 2^5 = 32 vertices. Here is a picture.
xrctElls <- lapply(ellCICombinationVectors(xrctInfMM$CI), ellEllipsoidFromVector)
ellPairsPlot(c(lapply(xrctElls, function(ell) ell$vector), list(ams205$vector)),
             color=c(replicate(length(xrctElls), "black"), "green"))

# If that plot doesn't do much for you, here are some more familiar plots.
ellHsuNadaiPlot(c(lapply(xrctElls, function(ell) ell$logA), list(ams205$logA)),
                color=c(replicate(length(xrctElls), "black"), "green"))
ellEqualVolumePlot(c(lapply(xrctElls, function(ell) ell$rotation), list(ams205$rotation)),
                   c(lapply(xrctElls, function(ell) ell$logA), list(ams205$logA)),
                   color=c(replicate(length(xrctElls), "white"), "green"))
ellEqualAreaPlot(c(lapply(xrctElls, function(ell) ell$rotation), list(ams205$rotation)),
                   c(lapply(xrctElls, function(ell) ell$logA), list(ams205$logA)),
                   color=c(replicate(length(xrctElls), "black"), "green"))

# The S-type bootstrap inference can be processed in exactly the same way. Its results are nearly identical. Let's skip it.

# The custom bootstrapping method produces a 5D confidence ellipsoid.
xrctElls <- lapply(ellHighEllipsoidVectors(xrctInf$covarInv, xrctInf$center, xrctInf$q095^2), 
                   ellEllipsoidFromVector)
ellPairsPlot(c(lapply(xrctElls, function(ell) ell$vector), list(ams205$vector)),
             color=c(replicate(length(xrctElls), "black"), "green"))
ellHsuNadaiPlot(c(lapply(xrctElls, function(ell) ell$logA), list(ams205$logA)),
                color=c(replicate(length(xrctElls), "black"), "green"))
ellEqualVolumePlot(c(lapply(xrctElls, function(ell) ell$rotation), list(ams205$rotation)),
                   c(lapply(xrctElls, function(ell) ell$logA), list(ams205$logA)),
                   color=c(replicate(length(xrctElls), "white"), "green"), simplePoints=TRUE)
ellEqualAreaPlot(c(lapply(xrctElls, function(ell) ell$rotation), list(ams205$rotation)),
                 c(lapply(xrctElls, function(ell) ell$logA), list(ams205$logA)),
                 color=c(replicate(length(xrctElls), "black"), "green"))



### WHAT HAVE WE LEARNED? ###

# Even though the plots make AMS look like the mean of XRCT, there are several statistical tests that show it isn't (for this rock sample).



### SEE ALSO ###

# In tutorialsBonus, ellInferenceMore.R carries out the same analysis, on the same rock sample, with XRCT ellipsoids of other compositions (medium and air). We come to the same conclusion.


