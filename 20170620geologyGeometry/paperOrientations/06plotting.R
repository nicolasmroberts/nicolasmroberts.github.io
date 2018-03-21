


# Figure 6: Outliers and gimbal lock.



### SYNTHETIC OUTLIER EXAMPLE ###

# If you have just restarted R, then set the working directory to our code directory (the one that contains subdirectories 'library', etc.).

# Execute the following lines to load our custom library and other libraries.
source("library/all.R")

# Load a synthetic data set of 28 fold orientations. Make two versions of the equal-area plot.
outlierData <- geoDataFromFile("data/synthFoldsOutlier.csv")
lineEqualAreaPlot(c(outlierData$pole, outlierData$direction),
                  shapes=c(replicate(length(outlierData$pole), "c"), replicate(length(outlierData$direction), "s")))
lineEqualAreaPlot(outlierData$direction,
                  curves=lapply(outlierData$pole, rayGreatCircle, 72),
                  shapes=replicate(length(outlierData$direction), "s"))

# An equal-angle orientation plot shows an outlier. 
oriEqualAnglePlot(outlierData$rotation)

# Mahalanobis distance suggests that the last datum is strangely placed relative to the others. That's the outlier.
mu <- oriMeanVariance(outlierData$rotation, group=oriLineInPlaneGroup)$mean
oris <- oriNearestRepresentatives(outlierData$rotation, mu, group=oriLineInPlaneGroup)
covar <- rotLeftCovariance(oris, mu)
covarInv <- solve(covar)
sapply(oris, rotMahalanobisNorm, mu, covarInv)

# The equal-area plot also reveals the outlier, if we use a miraculous choice of coloring.
lineEqualAreaPlot(c(outlierData$pole, outlierData$direction), colors=as.character(outlierData$color),
                  shapes=c(replicate(length(outlierData$pole), "c"), replicate(length(outlierData$direction), "s")))



### OTHER VIEWS OF THE OUTLIER ###

oriEulerAnglePlot(outlierData$rotation, group=oriLineInPlaneGroup,
                  backgroundColor="white", color="black", axesColors=c("black", "black", "black"))

oriEulerAnglePlot(outlierData$rotation, group=oriLineInPlaneGroup, showBoundary=TRUE, curveColor="black",
                  backgroundColor="white", color="black", axesColors=c("black", "black", "black"))

oriAxisAnglePlot(outlierData$rotation, group=oriLineInPlaneGroup,
                 backgroundColor="white", color="black", axesColors=c("black", "black", "black"), boundaryAlpha=0)

oriEqualAnglePlot(outlierData$rotation, group=oriLineInPlaneGroup,
                  backgroundColor="white", color="black", axesColors=c("black", "black", "black"), boundaryAlpha=0)

oriEqualVolumePlot(outlierData$rotation, group=oriLineInPlaneGroup,
                   backgroundColor="white", color="black", axesColors=c("black", "black", "black"), boundaryAlpha=0)

oriRodriguesPlot(outlierData$rotation, group=oriLineInPlaneGroup,
                 backgroundColor="white", color="black", axesColors=c("black", "black", "black"), boundaryAlpha=0)

# This plot was suggested by a reviewer: the pole to the plane through (A) the hinge line and (B) the axial plane.
# Yes, the outlier is revealed. But I don't expect that many geologists would think to plot this information.
lineEqualAreaPlot(lapply(outlierData$rotation, function(r) r[3,]))



### GIMBAL LOCK ###

# Load a synthetic data set of 17 quartz orientations.
gimbalData <- geoDataFromFile("data/synthGimbalLock.csv")

# Making an equal-area plot of this data type with our software takes some work.
gimbalCs <- lapply(gimbalData$rotation, function(r) r[3,])
gimbalCs <- lapply(gimbalCs, lower)
gimbalA1s <- lapply(gimbalData$rotation, function(r) r[1,])
gimbalA2s <- lapply(gimbalData$rotation, function(r) (oriTrigonalTrapezohedralGroup[[2]] %*% r)[1,])
gimbalA3s <- lapply(gimbalData$rotation, function(r) (oriTrigonalTrapezohedralGroup[[3]] %*% r)[1,])
rayEqualAreaPlot(c(gimbalCs, gimbalA1s, gimbalA2s, gimbalA3s),
                 colors=c(replicate(length(gimbalCs), "black"), replicate(length(gimbalA1s), "red"),
                          replicate(length(gimbalA2s), "green"), replicate(length(gimbalA3s), "blue")))

# Each orientation plots as six rotations, because of crystallographic symmetry.
oriAxisAnglePlot(gimbalData$rotation, group=oriTrigonalTrapezohedralGroup)
oriEqualVolumePlot(gimbalData$rotation, group=oriTrigonalTrapezohedralGroup)
oriEqualAnglePlot(gimbalData$rotation, group=oriTrigonalTrapezohedralGroup)
oriRodriguesPlot(gimbalData$rotation, group=oriTrigonalTrapezohedralGroup)

# In the Euler angle plot, four of the symmetric copies are depicted reasonably.
# However, two of the symmetric copies are smeared out along lines of gimbal lock.
oriEulerAnglePlot(gimbalData$rotation, group=oriTrigonalTrapezohedralGroup)


