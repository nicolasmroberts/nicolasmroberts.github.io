


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# This short tutorial demonstrates how to make plots of the 5D or 6D ellipsoid vectors. The plots, being 2D or 3D, cannot show all of the five or six degrees of freedom. So these plots, like the Hsu-Nadai and orientation plots before them, are quite imperfect. You should always view your ellipsoidal data in several different styles of plot, to lessen your chances of missing important features of your data set.

# For data we again use anisotropy of magnetic susceptibility (AMS) ellipsoids from rocks from the Troodos ophiolite, Cyprus (Titus et al.). We also use a synthetic data set.



### PRELIMINARY WORK ###

# If you have just restarted R, then set the working directory to the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc. Execute the following line of code to load our custom library and its dependencies.
source("library/all.R")

# Load that Cyprus AMS data set again, with the ellipsoids volume-normalized.
cyprusAMSData <- geoEllipsoidDataFromIRMFile("data/cyprus_AMS_groupF.tsv", doNormalize=TRUE)

# The preceding exercise argued that ellipsoid vectors are the best way to compute with ellipsoids. So it might behoove us to plot ellipsoids in this format. There's just one problem: The vectors are 5D or 6D, so we can't visualize them.
cyprusAMSData$vector



### PLOTTING ELLIPSOID VECTORS ###

# So we have to resign ourselves to plotting 2D or 3D projections. For starters, the following plot shows the 1st, 3rd, and 4th coordinates, ignoring the 2nd and 5th coordinates. You can change the 'c(1, 3, 4)' to any triple you want, to see different projections of the data. There are 5-choose-3 = 10 projections in all.
ellVectorPlot(c(1, 3, 4), cyprusAMSData$vector)

# The following code makes a 2D plot of the 2nd and 3rd coordinates, ignoring the 1st, 4th, and 5th coordinates. There are 5-choose-2 = 10 of these 2D projections too.
ellVectorPlot(c(2, 3), cyprusAMSData$vector)

# The following plot is like the preceding one, but showing all 10 combinations at once. Actually it shows each combination twice, for a total of 20 plots. (And if you use it on unnormalized ellipsoids you get 30 plots at once.)
ellPairsPlot(cyprusAMSData$vector)

# The vectors are difficult to understand. Each of their five (or six) components is a mixture of orientation and size-shape. The simplest property to note is that the normalized sphere plots at the origin in 5D space. For more properties, check out the bonus exercise ellVectors.R.



### SYNTHETIC EXAMPLE ###

# In an earlier exercise, we examined a synthetic example of plane-line pairs, in which there was an outlier that couldn't be seen in the equal-area plot. Now we do a similar example for ellipsoids. Load a synthetic data set of 31 normalized ellipsoids.
synthOutlierData <- geoDataFromFile("data/synthEllipsoidOutlier.csv")

# Here are the ellipsoid orientations. Does any datum look like an outlier?
ellEqualVolumePlot(synthOutlierData$rotation, synthOutlierData$a)

# Here are the ellipsoid shapes. Does any datum look like an outlier?
ellHsuNadaiPlot(synthOutlierData$logA)

# Here are all 2D projections of the ellipsoid vectors. Does any datum look like an outlier?
ellPairsPlot(synthOutlierData$vector)

# Here are the same three plots, with half of the data arbitrarily colored blue, half colored green, and the outlier colored red. You see that its orientation is similar to the blue ones, while its shape is similar to the green ones. And if you guessed any outlier at all in the vector plots, then you probably guessed this one. But even there the outlier is not completely obvious to me.
synthColors <- c(replicate(15, "blue"), replicate(15, "green"), "red")
ellEqualVolumePlot(synthOutlierData$rotation, synthOutlierData$a, colors=synthColors)
ellHsuNadaiPlot(synthOutlierData$logA, colors=synthColors)
ellPairsPlot(synthOutlierData$vector, colors=synthColors)



### WHAT HAVE WE LEARNED? ###

# Don't rely on any one plotting system. View your data in several kinds of plot, to lessen your chances of missing important features of your data set.



### SEE ALSO ###

# In tutorialsBonus, ellVectors.R plots some important special cases of ellipsoid vectors, to help you get a better sense for what these vectors mean. (But they will never be easy.)

# In tutorialsBonus, ellOutlier.R shows how I made this synthetic example.


