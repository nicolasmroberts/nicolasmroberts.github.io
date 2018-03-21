


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# In a clustering problem, our goal is to group the data into 'clumps' or clusters, such that the points within each cluster are more similar to each other than to the points in other clusters. Because there are many notions of what makes data similar, there are many clustering algorithms, with varying strengths. This tutorial discusses two kinds of clustering: DBSCAN and k-means. The methods are quite general, so we demonstrate them for various data types: geographic locations, directions, orientations, and ellipsoids.

# The data include paleomagnetic and dike data from the Troodos ophiolite, Cyprus (Titus et al.), slickenside orientations from central California (Woodring et al., 1940), and orthopyroxene and spinel shape preferred orientation (SPO) ellipsoids from New Caledonia (Titus et al., Chatzaras).



### PRELIMINARY WORK ###

# If you have just restarted R, then set the working directory to the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc. Execute the following line of code to load our custom library and its dependencies.
source("library/all.R")



### LOCATIONS ###

# Load some data from the Troodos ophiolite, Cyprus. Here are the stations in map view.
cyprusData <- geoDataFromFile("data/cyprusPmagNW.csv")
plot(x=cyprusData$easting, y=cyprusData$northing)

# Before we cluster, we build a matrix of distances between each pair of stations.
cyprusData$location <- thread(c, cyprusData$easting, cyprusData$northing)
cyprusLocDists <- matrixOfDistances(cyprusData$location, euclideanDistance)

# The DBSCAN clustering algorithm depends on two parameters chosen by the user: radius and minPoints. The algorithm declares a data point to be 'dense' if there are at least minPoints data points within the given radius of that point (including the point itself). Two data are declared to be in the same cluster if their distance is less than the given radius and at least one of the two data is dense. One consequence of this criterion is that some isolated points are not incorporated into any cluster. Our software reports those isolated points as a final pseudo-cluster.

# Here is a DBSCAN clustering of the Cyprus station locations. With radius = 2000 and minPoints = 4, we get one cluster of 10 points (shown in red) and 8 isolated points (magenta).
cyprusLocClus <- clusteringDBSCAN(cyprusLocDists, radius=2000, minPoints=4)
sapply(cyprusLocClus, length)
plot(x=cyprusData$easting, y=cyprusData$northing, col=clusteringHues(cyprusLocClus))

# In k-means clustering, the clusters are chosen so that each datum is closer to the average of its cluster than to the average of any other cluster. The user must supply a seed clustering, which the algorithm iteratively improves until it cannot be improved any further. Because differing seeds can end up at differing local optima for the problem, you might want to try multiple seeds. Here we use the DBSCAN cyprusLocClus as our seed.
cyprusLocClus <- clusteringKMeans(cyprusLocDists, cyprusLocClus)$gamma
sapply(cyprusLocClus, length)
plot(x=cyprusData$easting, y=cyprusData$northing, col=clusteringHues(cyprusLocClus))

# By the way, a Voronoi partition can be viewed as another kind of clustering. You pick a set of seeds, each of which produces a cluster. Each datum is placed into the cluster of the nearest seed. Here's an example with four seeds and hence four clusters.
cyprusLocSeeds <- list(c(476000, 3864000), c(476000, 3866000), c(480000, 3864000), c(484000, 3860000))
cyprusLocClus <- clusteringVoronoi(cyprusData$location, cyprusLocSeeds, euclideanDistance)
sapply(cyprusLocClus, length)
plot(x=cyprusData$easting, y=cyprusData$northing, col=clusteringHues(cyprusLocClus))



### DIRECTIONS ###

# In case you skipped the preceding section, load the Cyprus data again. View the paleomagnetic directions and precompute their angular differences.
cyprusData <- geoDataFromFile("data/cyprusPmagNW.csv")
rayEqualAreaPlot(cyprusData$direction)
cyprusPmagDists <- matrixOfDistances(cyprusData$direction, rayDistance)

# This DBSCAN clustering produces two clusters, of five points each, and eight isolated points. Tinker with radius and minPoints, to see how the results change.
cyprusPmagClus <- clusteringDBSCAN(cyprusPmagDists, radius=(10 * degree), minPoints=3)
sapply(cyprusPmagClus, length)
rayEqualAreaPlot(cyprusData$direction, colors=clusteringHues(cyprusPmagClus))

# k-means clustering, seeded from the preceding clustering for convenience.
cyprusPmagClus <- clusteringKMeans(cyprusPmagDists, cyprusPmagClus)$gamma
sapply(cyprusPmagClus, length)
rayEqualAreaPlot(cyprusData$direction, colors=clusteringHues(cyprusPmagClus))

# View the dike poles and precompute their angular differences. Notice that, while paleomagnetic directions are rays, dike poles are lines.
lineEqualAreaPlot(cyprusData$pole)
cyprusDikeDists <- matrixOfDistances(cyprusData$pole, lineDistance)

# This DBSCAN clustering produces one cluster, of size 16, and two isolated points. Tinker with radius and minPoints.
cyprusDikeClus <- clusteringDBSCAN(cyprusDikeDists, radius=(15 * degree), minPoints=4)
sapply(cyprusDikeClus, length)
lineEqualAreaPlot(cyprusData$pole, colors=clusteringHues(cyprusDikeClus))

# k-means clustering, seeded from the preceding clustering for convenience. It doesn't change that clustering at all.
cyprusDikeClus <- clusteringKMeans(cyprusDikeDists, cyprusDikeClus)$gamma
sapply(cyprusDikeClus, length)
lineEqualAreaPlot(cyprusData$pole, colors=clusteringHues(cyprusDikeClus))



### ORIENTATIONS ###

# Load central California slickenside orientations. The rake describes the movement of the hanging wall. As usual, we apply geoPoleVorticityFromPoleHanging to get the orientations into a better-behaved format. Here are the orientations in equal-volume.
slickData <- geoDataFromFile("data/cali_all_kh_faults_latlong.tsv")
slickRots <- lapply(slickData$rotation, geoPoleVorticityFromPoleHanging)
slickDists <- matrixOfDistances(slickRots, oriDistance, verbose=FALSE, oriRayInPlaneGroup)
oriEqualVolumePlot(slickRots, group=oriRayInPlaneGroup)

# This DBSCAN clustering detects two big clusters, three small clusters, and 31 isolated points. I like these results.
slickClus <- clusteringDBSCAN(slickDists, radius=(15 * degree), minPoints=3)
sapply(slickClus, length)
oriEqualVolumePlot(slickRots, group=oriRayInPlaneGroup, colors=clusteringHues(slickClus))

# k-means clustering, seeded from the preceding clustering for convenience. Notice that this clustering splits some big, obvious clumps in the data. k-means clustering often disagrees with my intuition like this.
slickClus <- clusteringKMeans(slickDists, slickClus)$gamma
sapply(slickClus, length)
oriEqualVolumePlot(slickRots, group=oriRayInPlaneGroup, colors=clusteringHues(slickClus))



### ELLIPSOIDS ###

# Load SPO data from New Caledonia. For the sake of this demonstration, we lump the orthopyroxene and spinel SPOs together.
source("data/newcalOPXSpinelSPO.R")
newcalElls <- c(lapply(ncOrthos, function(station) station$ortho), lapply(ncSpinels, function(station) station$spinel))
length(newcalElls)
newcalElls[[1]]

# Inspect the data and precompute the distances.
newcalVecs <- lapply(newcalElls, function(ell) ell$vector)
ellPairsPlot(newcalVecs)
ellVectorPlot(c(1, 2, 3), newcalVecs)
newcalDists <- matrixOfDistances(newcalVecs, euclideanDistance)

# DBSCAN clustering. Four clusters and lots of isolated points. The results are hard to interpret.
newcalClus <- clusteringDBSCAN(newcalDists, radius=0.2, minPoints=3)
sapply(newcalClus, length)
ellPairsPlot(newcalVecs, colors=clusteringHues(newcalClus))
ellVectorPlot(c(1, 2, 3), newcalVecs, colors=clusteringHues(newcalClus))

# k-means clustering, starting from the preceding clustering for convenience. As in the k-means clustering of quartz orientations, the chosen clusters don't agree with my intuition.
newcalClus <- clusteringKMeans(newcalDists, newcalClus)$gamma
sapply(newcalClus, length)
ellPairsPlot(newcalVecs, colors=clusteringHues(newcalClus))
ellVectorPlot(c(1, 2, 3), newcalVecs, colors=clusteringHues(newcalClus))



### WHAT HAVE WE LEARNED? ###

# Clustering can help you dissect data sets for subsequent analysis. There is no one canonical, objective way to cluster. There are multiple approaches, each of which depends on choices by the user.


