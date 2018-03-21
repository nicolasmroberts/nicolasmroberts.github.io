


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# This tutorial goes through some examples of how to load a spreadsheet of ellipsoid data into R.



### PRELIMINARY WORK ###

# If you have just restarted R, then set the working directory to the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc. Execute the following line of code to load our custom library and its dependencies.
source("library/all.R")



### TREND-PLUNGE PAIRS AND SEMI-AXIS LENGTHS ###

# Open the file 'data/synthEllipsoids.csv' in your favorite spreadsheet program. Open it in R too, for direct comparison. You can see that geoDataFromFile has computed tons of extra stuff.
importData <- geoDataFromFile("data/synthEllipsoids.csv")
importData

# What happens here is: trend1, plunge1 are taken to describe one semi-axis of the ellipsoid. And a1 is taken as its semi-axis length. And trend2, plunge2, a2 are taken to describe another semi-axis. And then the other semi-axis is perpendicular to the first two, of length a3. (trend3 and plunge3 are ignored, in fact.) Then you can compute everything about the ellipsoid.

# Are these ellipsoids volume-normalized? Yes, because a1 * a2 * a3 == 1 for each one. Or equivalently log a1 + log a2 + log a3 == 0.
importData$a
sapply(importData$a, prod)
sapply(importData$logA, sum)

# But geoDataFromFile doesn't know that they're normalized. It's treating them as having six degrees of freedom.
sapply(importData$vector, length)

# So let's inform geoDataFromFile that it should regard the ellipsoids as normalized.
importData <- geoDataFromFile("data/synthEllipsoids.csv", doNormalize=TRUE)
sapply(importData$a, prod)
sapply(importData$logA, sum)
sapply(importData$vector, length)

# And then you're ready to do statistics. For example, here's the geometric mean ellipsoid.
ellMean(importData$vector)



### AN AMS DATA FILE ###

# Open the file 'data/cyprus_AMS_groupF.tsv' in your favorite spreadsheet application, and open it in R for comparison, using this special loader function.
importData <- geoEllipsoidDataFromIRMFile("data/cyprus_AMS_groupF.tsv")
importData

# This is a file of AMS data from Cyprus, that we processed at the Institute for Rock Magnetism in Minneapolis. The ellipsoid in question is the magnitude ellipsoid, not the susceptibility ellipsoid (Hrouda, 1982). Dmax..in.situ, Imax..in.situ are the trend and plunge of the long semi-axis, and max is its semi-axis length. Similarly, Dmin..in.situ, Imin..in.situ, and min describe the short semi-axis. So geoEllipsoidDataFromIRMFile dissects the data frame accordingly.

# Is it normalized? Yikes, definitely not.
sapply(importData$a, prod)

# Anyway, you're ready to do statistics. For example, here's the geometric mean ellipsoid.
ellMean(importData$vector)



### AN X-RAY COMPUTED TOMOGRAPHY DATA FILE ###

# Open the file 'data/HU150205mafic.tsv' in your favorite spreadsheet application, and open it in R for comparison, using this special loader function. The file is too big to view (4001 data). So let's just look at the first row of the data frame, to see what it's like.
importData <- geoEllipsoidDataFromAvizoFile("data/HU150205mafic.tsv")
nrow(importData)
importData[1,]

# The file contains X-ray computed tomography data from some mafic clasts. The file was produced by an image-processing program made by Avizo. Each row of the file describes an ellipsoid. The eigenvectors are the semi-axis directions, in Cartesian coordinates. The eigenvalues are the squared semi-axis lengths, according to Avizo tech support.

# Are the ellipsoid volume-normalized? No.
hist(sapply(importData$a, prod))
hist(sapply(importData$logA, sum))

# Should they be normalized? It depends on the question you want to study. Anyway, here's the mean.
ellMean(importData$vector)



### ELLIPSOID VECTORS ###

# Open the file 'data/fabricInitials.csv' in your favorite spreadsheet application, and open it in R for comparison. The file is too big to view conveniently (1000 data). So let's just look at the first row of the data frame, to see what it's like.
importData <- geoDataFromFile("data/fabricInitials.csv")
nrow(importData)
importData[1,]

# The file is a synthetic data set, used in one of our ellipsoid exercises. Each row of the file gives the five entries of a normalized log-ellipsoid vector. From that information, geoDataFromFile can reconstruct the entire ellipsoid.

# Here's the normalization and the mean.
sapply(importData$a, prod)
ellMean(importData$vector)



### WHAT HAVE WE LEARNED? ###

# There are various ways to import ellipsoidal data sets into R. Just be careful about how your data are formatted, because there are a number of conventions in use, especially for the semi-axis lengths.



### SEE ALSO ###

# Bonus tutorial dirImportation.R does a similar exercise for directional data. It also goes into more depth about how to 'repair' a data set.

# Bonus tutorial oriImportation.R does a similar exercise for orientational data.


