


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# This tutorial continues 6inference.R from tutorialsEllipsoids. We apply exactly the same methods to other data sets.

# The data sets are all unpublished. There is one anisotropy of magnetic susceptibility (AMS) ellipsoid. Then there are 'clast' ellipsoids produced by X-ray computed tomography (XRCT) on the same rocks. Some of the clasts have 'medium' composition (between felsic and mafic). Other clasts are actually voids in the rock.



### PRELIMINARY WORK ###

# If you have just restarted R, then set the working directory to the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc. Execute the following line of code to load our custom library and its dependencies.
source("library/all.R")

# Load one AMS ellipsoid from a file. It's the same AMS ellipsoid from ellipsoid tutorial 6inference.R.
ams205 <- geoEllipsoidDataFromAGICOFile("data/HU150205ams.tsv", meanSuscepts=c(8.957E-03), doNormalize=TRUE)
ams205 <- listFromDataFrame(ams205)[[1]]



### HU150205 MEDIUM-COMPOSITION DATA ###

# In brief, here are the same four hypothesis tests for clasts of 'medium' composition in the same sample.
xrct205Medium <- geoEllipsoidDataFromAvizoFile("data/HU150205medium.tsv", doNormalize=TRUE)
nrow(xrct205Medium)
ellHotellingT2Inference(xrct205Medium$vector, ams205$vector)
ellBootstrapMMInference(xrct205Medium$vector, ams205$vector, numBoots=1000)
ellBootstrapSInference(xrct205Medium$vector, ams205$vector, numBoots=1000)
ellBootstrapInference(xrct205Medium$vector, numBoots=1000)$pvalue(ams205$vector)



### HU150205 VOID DATA ###

# Here are the same four hypothesis tests for voids (air pockets) in the same sample.
xrct205Air <- geoEllipsoidDataFromAvizoFile("data/HU150205air.tsv", doNormalize=TRUE)
nrow(xrct205Air)
ellHotellingT2Inference(xrct205Air$vector, ams205$vector)
ellBootstrapMMInference(xrct205Air$vector, ams205$vector, numBoots=1000)
ellBootstrapSInference(xrct205Air$vector, ams205$vector, numBoots=1000)
ellBootstrapInference(xrct205Air$vector, numBoots=1000)$pvalue(ams205$vector)



### WHAT HAVE WE LEARNED? ###

# In all of these examples, we reject the hypothesis that the AMS ellipsoid is the mean of the population from which the clast ellipsoids were drawn. So it appears that we cannot use AMS as a simple, easy proxy for XRCT.


