


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# This tutorial introduces a few other data types: character strings, matrices, and logicals. We don't go into much depth on any one. We're just trying to give you an idea of what's out there, to de-mystify R a bit.



### PRELIMINARY WORK ###

# If you have just restarted R, then set the working directory to the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc. Execute the following line of code to load our custom library and its dependencies.
source("library/all.R")



### CHARACTER (STRINGS OF TEXT) ###

# Let's load some foliation-lineation pairs from the western Idaho shear zone.
wiszData <- geoDataFromFile("data/wiszFollins.tsv")

# In that line of code, the file name is its own kind of R object, that we haven't discussed yet: a string of characters.
"data/wiszFollins.tsv"

# This data frame contains lots of information. Some of it was in the file. Some of it was computed by geoDataFromFile.
wiszData

# When we analyze geologic data with R, one of the main uses of strings is for specifying colors in plots. This plot shows the field sites in map view, colored green.
plot(x=wiszData$easting, y=wiszData$northing, col="green")

# The hues function from our library helps you generate colors based on a vector of numbers. Check this out.
wiszData$southeasting <- wiszData$easting - wiszData$northing
wiszData$color <- hues(wiszData$southeasting)
wiszData$color
plot(x=wiszData$easting, y=wiszData$northing, col=wiszData$color)

# Coloring a map by southeasting isn't very useful. But coloring an equal-area plot by southeasting might be useful, for helping us detect geographic patterns in dike poles, if there are any.
lineEqualAreaPlot(wiszData$pole, colors=wiszData$color)

# Here's another place in which you might encounter strings: the 'station' column of the WISZ data. This column contains the names of the field sites, as chosen by the geologists who visited them. They are strings, but R goes one step further, interpreting them as 'levels' of a 'factor' --- roughly speaking, different parts of a study.
wiszData$station

# If you don't want R to interpret the station names as levels of a factor, then force them to be mere strings of characters, like this. But most of the time I find myself not caring about the distinction.
wiszData$station <- as.character(wiszData$station)
wiszData$station



### MATRICES ###

# Structural geology uses a lot of matrices, for describing rotations, ellipsoids, deformations, stresses, etc. So our R library uses matrices heavily. For example, here are the orientation matrices of the WISZ foliation-lineation pairs. (See the orientation exercise 1lineInPlane.R for an explanation of what they mean.)
wiszData$rotation

# Here is the 2nd matrix in that list, and its third row, and its first column, and its transpose.
a <- wiszData$rotation[[2]]
a
a[3,]
a[,1]
t(a)

# And here are the determinant and eigensystem of the 12th matrix in that list.
det(wiszData$rotation[[12]])
eigen(wiszData$rotation[[12]])

# And here is the matrix product of the 7th and 8th matrices.
wiszData$rotation[[7]] %*% wiszData$rotation[[8]]



### LOGICAL (TRUE AND FALSE) ###

# The final data type to be discussed here is logical, which has only two possible values: TRUE and FALSE. For example:
9 / 2 == 6
sin(pi / 3) < 0.6

# One of my favorite uses of logicals in R is to dissect vectors, lists, and data frames. For example, which WISZ stations are north of northing 4998000?
wiszData$northing > 4998000

# Let's store that vector of logicals, so that we can use it again and again.
isNorth <- wiszData$northing > 4998000
isNorth

# The '!' operator means 'not'; it switches TRUEs with FALSEs.
!isNorth

# Make separate plots of the dike poles in the northern and southern parts of the field area.
lineEqualAreaPlot(wiszData$pole[isNorth])
lineEqualAreaPlot(wiszData$pole[!isNorth])

# You can even split the wiszData data frame into two parts along its rows, like this.
wiszDataNorth <- wiszData[isNorth,]
wiszDataSouth <- wiszData[!isNorth,]
wiszDataNorth
wiszDataSouth



### BRANCHES (IF-ELSE) ###

# Actually, the biggest use of logicals is in branches (if-else statements). After storing values in memory, branches are the second big stepping stone to programming. In the code block below, one of two paths is taked, depending on whether wiszDataNorth or wiszDataSouth has more rows. The overall effect of this code block is to set wiszDataNew to whichever data frame is larger.
if (nrow(wiszDataNorth) > nrow(wiszDataSouth)) {
  wiszDataNew <- wiszDataNorth
} else {
  wiszDataNew <- wiszDataSouth
}
wiszDataNew

# Here's another example. Structual geologists often like to view directions as points on the lower hemisphere. But sometimes mathematical calculations produce direction vectors on the upper hemisphere. In the following code block, we randomly generate a vector v, and then flip it to the lower hemisphere if needed. Run this code block a few times, to see how its behavior depends on the randomly generated z-coordinate of v.
v <- rnorm(3)
v <- v / sqrt(sum(v * v))
v
if (v[[3]] > 0)
  v <- -v
v



### WHAT HAVE WE LEARNED? ###

# Character strings, matrices, and logicals are three other data types that show up in R programs. Logicals play a key role in branches, which are a key part of programming.


