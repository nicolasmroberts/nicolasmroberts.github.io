


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# This tutorial goes through some examples of how to load a spreadsheet of orientation data into R --- specifically, so that it can be processed with our geoDataFromFile function. You can choose to enter your orientations in several formats.



### PRELIMINARY WORK ###

# If you have just restarted R, then set the working directory to the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc. Execute the following line of code to load our custom library and its dependencies.
source("library/all.R")



### STRIKE-DIP, TREND-PLUNGE ###

# Suppose that our data consist of planes and lines. We could specify the planes as strikes and dips, and the lines as trends and plunges. For an example, open the file 'data/synthFoldsOutlier.csv' in your favorite spreadsheet application. Also load it into R, so that we can compare the two views of these data.
importData <- geoDataFromFile("data/synthFoldsOutlier.csv")
importData

# In the spreadsheet, you see that the first row consists of headers, and every row after that consists of data. (This is a basic requirement of R's read.table function. For more information, see bonus tutorial dirImporting.R.) And you can see exactly these headers in the R data frame, except that 'Deg' has been appended to a few of them.
names(importData)

# Here's just the first row of the data frame. When geoDataFromFile sees strike and dip, it automatically computes a pole vector. When geoDataFromFile sees trend and plunge, it automatically computes a direction vector.
importData[1,]

# Once geoDataFromFile has a pole and a direction, it also computes a rotation matrix from them. Let's compare them, for the first row of the data frame.
importData$pole[[1]]
importData$direction[[1]]
importData$rotation[[1]]

# See how the pole is the first row of the rotation, and the direction is the second row? That's exactly the format used in our orientation tutorials, such as 1linesInPlanes.R. Actually, the rotation matrix's rows are slightly different from the pole and direction, because the rows are forced to be perfectly perpendicular. But you can't tell in this example. Anyway, you can now do all of our usual stuff, such as:
oriFrechetMean(importData$rotation, group=oriLineInPlaneGroup)

# Here's a second example. Open the 'data/wiszFollins.tsv' file in your spreadsheet application, and load it into R for comparison.
importData <- geoDataFromFile("data/wiszFollins.tsv")
importData

# Once again, geoDataFromFile has computed pole, direction, and rotation. Let's inspect them from the first row of the data frame. See how the first two rows of the matrix are slightly different from the pole and direction? That's because they've been forced to perpendicular. The rotation is a 'cleaned up' version of the data.
importData$pole[[1]]
importData$direction[[1]]
importData$rotation[[1]]

# Here's a third example. Let's inspect the first row of this paleomagnetic data set. The first two rows of the rotation are really different from the pole and direction. Why? Well, the pole and direction must have been pretty far away from perpendicular. But why? Because the pole is for a dike, and the direction is for paleomagnetism, and there is no reason why they should be perpendicular.
importData <- geoDataFromFile("data/cyprusPmagNW.csv")
importData
importData$pole[[1]]
importData$direction[[1]]
importData$rotation[[1]]

# The rotation matrix is still a valid mathematical entity. The pole and direction can be regarded as a non-orthogonal rigid frame, and the rotation matrix expresses the orientation of that frame. You might even want to study it, if you thought that dike poles and paleomagnetic directions rotated together, for example in rigid blocks of rock.



### TREND-PLUNGE, TREND-PLUNGE ###

# Alternatively, you can express an orientation as two rays or lines, each in trend-plunge format. If they are called 'trend1', 'plunge1', 'trend2', 'plunge2', then geoDataFromFile will automatically build a rotation matrix for you. Those vectors are the first two rows, 'corrected' so that they are perpendicular.

# For example, open the file 'data/cyprus_AMS_groupF.tsv' in your favorite spreadsheet application, and open it here for comparison.
importData <- geoDataFromFile("data/cyprus_AMS_groupF.tsv")
importData

# This is a data set of AMS ellipsoids. The short axis has trend-plunge (Dmin..in.situ, Imin..in.situ), and the long axis has trend-plunge (Dmax..in.situ, Imax..in.situ). To get the orientation, you have three options.

# Option A. Make a copy of the data file. Change the headers 'Dmin..in.situ', 'Imin..in.situ' to 'trend1', 'plunge1'. Change the headers 'Dmax..in.situ', 'Imax..in.situ' to 'trend2', 'plunge2'. Then load that version of the file with geoDataFromFile. Done.

# Option B. This is much like Option A, but done by hand in R. Maybe you can guess how it works.
importData <- geoDataFromFile("data/cyprus_AMS_groupF.tsv")
importFunc <- function(t1Deg, p1Deg, t2Deg, p2Deg) {
  row1 <- geoCartesianFromTrendPlungeDeg(c(t1Deg, p1Deg))
  row2 <- geoCartesianFromTrendPlungeDeg(c(t2Deg, p2Deg))
  rotProjectedMatrix(rbind(row1, row2, cross(row1, row2)))
}
importData$rotation <- thread(importFunc,
                              importData$Dmin..in.situ, importData$Imin..in.situ,
                              importData$Dmax..in.situ, importData$Imax..in.situ)

# Let's check that it worked correctly, for the first row of the data frame.
geoCartesianFromTrendPlungeDeg(c(importData$Dmin..in.situ[[1]], importData$Imin..in.situ[[1]]))
geoCartesianFromTrendPlungeDeg(c(importData$Dmax..in.situ[[1]], importData$Imax..in.situ[[1]]))
importData$rotation[[1]]

# Option C. Treat the entire data set as an ellipsoid data set, to be loaded with geoEllipsoidDataFromIRMFile. It computes the rotations and a bunch of other information. For example, here is the first orientation. It should exactly match the one we computed in Option B.
importData <- geoEllipsoidDataFromIRMFile("data/cyprus_AMS_groupF.tsv")
importData$rotation[[1]]

# No matter which option you've chosen, you're now ready to go. For example, here is the mean orientation.
oriFrechetMean(importData$rotation, group=oriLineInPlaneGroup)



### STRIKE-DIP-RAKE (AKA PITCH) ###

# There seem to be two conventions for strike-dip-rake.

# First, geoDataFromFile assumes the right-hand rule. The strike of a plane is stated so that the plane dips to the right. Rake is measured down the plane from the strike. In other words, rake is measured clockwise, when the plane is viewed from above. Consequently, strikes vary from 000 degrees to 360 degrees, dips vary from 0 degrees to 90 degrees, and rakes vary from 0 degrees to 360 degrees (for rays, such as hanging wall movement directions) or 0 degrees to 180 degrees (for lines, such as lineations within foliations). 

# Second, some geologists prefer to state strike, dip, dip direction, rake, rake direction, and sense of slip (where applicable). The strike, dip, and dip direction specify the plane. The rake is measured down the plane from one of the two strikes; the rake direction indicates which of the two strikes to measure from. At this point, the rake has specified a line. If we wish to turn the line into a ray --- for example, to distinguish between a normal fault and thrust fault --- then the sense of slip tells how to.

# Let's do an example: 'strike 030, dip 80, dip direction NW, rake 40 NE, thrust' (all in degrees). First we impose the right-hand rule, changing the strike to 210 degrees, so that a dip of 80 degrees goes NW. The NE rake direction indicates that rake is measured from 030 degrees, not 210 degrees. So the line is 40 degrees down the plane (counterclockwise) from the 030 degrees strike. But the thrust says to reverse this direction. So the ray is now 140 degrees clockwise from the 030 degrees strike. So it is 320 degrees clockwise from the 210 degrees strike. So the final strike-dip-rake is (210, 80, 320), in degrees.

# For our first example, open 'data/cali_all_kh_faults_latlong.tsv' in your favorite spreadsheet application, and open it in R for comparison. You can see that geoDataFromFile has appended 'Deg' to some of the headers. And it has added pole, direction, and rotation fields.
importData <- geoDataFromFile("data/cali_all_kh_faults_latlong.tsv")
importData

# These data are faults with slip directions. The rake is describing the movement direction of the hanging wall. As we explain in orientation tutorial 3otherOrientations.R, this is not a good convention for orientations of faults with slip, because the hanging wall is undefined for vertical faults. To clean up the orientations, you apply geoPoleVorticityFromPoleHanging. Then do whatever statistics you want.
importData$cleaned <- lapply(importData$rotation, geoPoleVorticityFromPoleHanging)
oriFrechetMean(importData$cleaned, group=oriRayInPlaneGroup)

# For our second example, open 'data/synthCurrentRipples.csv' in your spreadsheet application and in R, for comparison.
importData <- geoDataFromFile("data/synthCurrentRipples.csv")
importData

# This is a synthetic data set of current ripple orientations. There is no hanging wall involved. So there is no bad behavior at vertical planes. So you do NOT apply geoPoleVorticityFromPoleHanging.

# The pole is the pole to bedding. The younging direction is known, so this pole is a ray, not a line. The rake direction describes the water flow along that plane. Because these are asymmetric ripples, the direction is a ray, not a line. Consequently, current ripples have trivial symmetry. They are effectively just rotations. So here is their mean.
rotFrechetMean(importData$rotation)



### CARTESIAN COORDINATES ###

# geoDataFromFile recognizes one more data format, in which two directions are given as Cartesian vectors of length 1. The first vector is taken from the x1, y1, z1 columns of the data file. The second vector is taken from the x2, y2, z2 columns. They are placed into the first two rows of a rotation matrix, which is then 'corrected' to be special orthogonal.

# For our first example, open 'data/synthGimbalLock.csv' in your favorite spreadsheet application, and in R. Inspect the first row of the data frame.
importData <- geoDataFromFile("data/synthGimbalLock.csv")
c(importData$x1[[1]], importData$y1[[1]], importData$z1[[1]])
c(importData$x2[[2]], importData$y2[[1]], importData$z2[[1]])
importData$rotation[[1]]

# These are quartz orientations. So here is their mean.
oriFrechetMean(importData$rotation, group=oriTrigonalTrapezohedralGroup)

# For our second example, open 'data/moine_one_grainABCxyz.tsv' as a spreadsheet and in R. Inspect the first row of the data frame.
importData <- geoDataFromFile("data/moine_one_grainABCxyz.tsv")
c(importData$x1[[1]], importData$y1[[1]], importData$z1[[1]])
c(importData$x2[[2]], importData$y2[[1]], importData$z2[[1]])
importData$rotation[[1]]

# By coincidence, these data are also quartz orientations. So here is their mean.
oriFrechetMean(importData$rotation, group=oriTrigonalTrapezohedralGroup)




### WHAT HAVE WE LEARNED? ###

# There are four formats, in which you can enter orientation data for importation with geoDataFromFile.



### SEE ALSO ###

# Bonus tutorial dirImportation.R does a similar exercise for directional data. It also goes into more depth about how to 'repair' a data set.

# Bonus tutorial ellImportation.R does a similar exercise for ellipsoidal data.


