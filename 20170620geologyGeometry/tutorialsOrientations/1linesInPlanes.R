


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# In this tutorial we convert orientations of plane-line pairs into rotations. This data type has a 4-fold symmetry, in that each orientation corresponds to four rotations. Examples include foliations containing lineations, and fold axial planes containing hinge lines. Examples also include orientations of triaxial ellipsoids, principal stress directions, and double-couple earthquake focal mechanisms.

# The data used in this tutorial include foliation-lineation pairs from the western Idaho shear zone (Giorgis and Tikoff, 2004) and SPO ellipsoid orientations from New Caledonia (Titus et al., Chatzaras).



### PRELIMINARY WORK ###

# If you have just restarted R, then set the working directory to the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc. Execute the following line of code to load our custom library and its dependencies.
source("library/all.R")



### STRIKE-DIP-TREND-PLUNGE TO ROTATION ###

# As an example, let's talk about the plane with strike-dip (348.2 degrees, 65.4 degrees), and the line with trend-plunge (47.6 degrees, 62.4 degrees). They are perpendicular to each other. As in our earlier sections about directional statistics, we express the plane pole and the line direction as unit vectors in 3D Cartesian coordinates.
myPole <- geoCartesianFromStrikeDipDeg(c(348.2, 65.4))
myPole
myDirection <- geoCartesianFromTrendPlungeDeg(c(47.6, 62.4))
myDirection

# Then we compute their cross product. Because the pole and the direction are unit and perpendicular, the cross product is also unit and perpendicular to both of them.
myPoleCrossDir <- cross(myPole, myDirection)
myPoleCrossDir

# Then we put these three vectors into the rows of a 3x3 matrix R.
myR <- rbind(myPole, myDirection, myPoleCrossDir)
myR

# In theory, this matrix R is special orthogonal, meaning that R R^T == I and det(R) == 1.
myR %*% t(myR)
det(myR)

# In practice, the fact that the pole and direction are not exactly perpendicular means that R is not exactly special orthogonal. Let's clean it up a bit, by 'projecting' it onto the nearest special orthogonal matrix.
myR <- rotProjectedMatrix(myR)
myR
myR %*% t(myR)
det(myR)

# Like all matrices, special orthogonal matrices can be viewed as transformations of space. In the special orthogonal case, the transformation is a rotation. So special orthogonal matrices are also called rotation matrices. Viewing the matrix R in this way helps us understand its meaning: R is the rotation that takes the plane pole to the positive x-axis and the line direction to the positive y-axis.

# The 'animation' code below shows a 3D plot of 5 frames of an animation. The animation is of R rotating the given plane-line pair into alignment with the y-z-plane and the y-axis. Superimposed on the plot is a thin white line. This line is the axis of the rotation. It is scaled so that its length indicates the amount of rotation about that axis.
animation <- function(r, size, numSteps) {
  ua <- rotAxisAngleFromMatrix(r)
  scaledAxis <- list(c(0, 0, 0), ua[1:3] * ua[[4]])
  line <- lapply(seq(from=-size, to=size, by=0.1), function(s) c(0, s, 0))
  square <- list(list(c(0, -size, -size), c(0, size, -size), c(0, size, size)),
                 list(c(0, -size, -size), c(0, size, size), c(0, -size, size)))
  points <- list()
  triangles <- list()
  for (s in seq(from=0, to=1, length.out=numSteps)) {
    q <- rotMatrixFromAxisAngle(c(ua[1:3], s * ua[[4]]))
    points <- c(points, lapply(line, function(v) q %*% v))
    triangles <- c(triangles, lapply(square, function(tri) lapply(tri, function(v) q %*% v)))
  }
  plot3D(radius=pi, points=points, curves=list(scaledAxis), triangles=triangles, curveWidth=2)
}
animation(myR, 2, 5)

# In practice, you would never build these rotation matrices by hand. For example, here is a data set of foliation-lineations from the western Idaho shear zone (Giorgis and Tikoff, 2004). The file contains columns for strike, dip, trend, and plunge. The geoDataFromFile function automatically builds the corresponding rotation matrix. The example we've been working above is just the 15th line of the following data file.
wiszData <- geoDataFromFile("data/wiszFolLins.tsv")
wiszData
wiszData[15,]
wiszData$rotation[[15]]

# At this point you might be wondering: What is the point of converting strike-dip-trend-plunge (or strike-dip-rake, or any other popular angles) to a rotation matrix? Well, researchers in numerous fields have been studying how to do statistics on rotations for over four decades. If we want to do statistics on plane-line pairs, it is much easier to convert to rotations and use their work, than to reinvent a bunch of statistics in terms of strike-dip-trend-plunge.



### FOUR ROTATIONS, ACTUALLY ###

# Now some bad news: R is not the only rotation matrix that returns the plane to the y-z-plane and the line to the y-axis. For we can negate the direction vector without changing its geometric meaning. That means that the following matrix R is an equally valid rotational representation of the original plane-line pair.
myR <- rotProjectedMatrix(rbind(myPole, -myDirection, cross(myPole, -myDirection)))
myR
animation(myR, 2, 5)

# But we can also negate the pole vector without changing its physical meaning:
myR <- rotProjectedMatrix(rbind(-myPole, myDirection, cross(-myPole, myDirection)))
myR
animation(myR, 2, 5)

# And finally we can negate both the direction and the pole, to get a fourth equivalent rotation:
myR <- rotProjectedMatrix(rbind(-myPole, -myDirection, cross(-myPole, -myDirection)))
myR
animation(myR, 2, 5)

# The four rotations that represent a given plane-line pair have a certain symmetry to them. They are related by a symmetry group, much like a crystallographic point group in mineralogy. This symmetry group is pre-programmed into our software as oriLineInPlaneGroup. Here are the four rotations again.
oriLineInPlaneGroup[[1]] %*% myR
oriLineInPlaneGroup[[2]] %*% myR
oriLineInPlaneGroup[[3]] %*% myR
oriLineInPlaneGroup[[4]] %*% myR



### ORTHOGONAL LINE TRIPLES ###

# A triaxial ellipsoid has three axes, which are perpendicular to each other and ordered (longest-to-shortest or shortest-to-longest, depending on the convention in use). Similarly, a state of stress has three principal stress directions, which are perpendicular and ordered. These data types can be converted into rotations much as in the previous section: We place the three direction vectors into the rows of a matrix. Each direction can be negated without changing its geometric meaning. So we have exactly the same 4-fold symmetry as in the previous section.

# For example, Vasileios Chatzaras analyzed shape preferred orientation (SPO) of spinel grains in rocks from the Massif du Sud, New Caledonia. Here are the trends and plunges of the three axes of an SPO ellipsoid, and the four representative rotations.
spoTrendPlungeDegs <- c(246, 16, 338, 06, 086, 73)
spoAxis1 <- geoCartesianFromTrendPlungeDeg(spoTrendPlungeDegs[1:2])
spoAxis2 <- geoCartesianFromTrendPlungeDeg(spoTrendPlungeDegs[3:4])
spoAxis3 <- geoCartesianFromTrendPlungeDeg(spoTrendPlungeDegs[5:6])
spoR <- rotProjectedMatrix(rbind(spoAxis1, spoAxis2, spoAxis3))
oriLineInPlaneGroup[[1]] %*% spoR
oriLineInPlaneGroup[[2]] %*% spoR
oriLineInPlaneGroup[[3]] %*% spoR
oriLineInPlaneGroup[[4]] %*% spoR



### EARTHQUAKE FOCAL MECHANISMS? ###

# A double-couple earthquake focal mechanism consists of (A) two perpendicular planes, which divide the ambient space into four quadrants, and (B) a labeling of two opposite quadrants as extensional and the other two quadrants as contractional. How can this data type be represented as a plane-line pair or an orthogonal line triple? (Therefore it has the same 4-fold symmetry as above.)



### WHAT HAVE WE LEARNED? ###

# The orientation of an object in space can be expressed as a rotation --- the rotation that takes the object back to some reference orientation. For many orientational data types there is more than one representative rotation, related by a symmetry group.


