


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# In this tutorial we start thinking about analyzing ellipsoids. Various obstacles arise. For example, how can we do statistics on ellipsoid size-shape? Simple calculations don't work well, fundamentally because size-shapes don't form a vector space. How can we do statistics on orientation? We have lots of tools for that, but those tools don't work well when spheroids (or near-spheroids) are present in the data. We would really like to treat ellipsoids holistically, as unified objects, rather than piecemeal. Fortunately, there is a way to treat ellipsoids holistically as elements in a vector space. This idea forms the basis for all of the methods in subsequent sections.

# In addition to small synthetic data sets, this tutorial uses anisotropy of magnetic susceptibility (AMS) ellipsoids from rocks from the Troodos ophiolite, Cyprus (Titus et al.).



### PRELIMINARY WORK ###

# If you have just restarted R, then set the working directory to the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc. Execute the following line of code to load our custom library and its dependencies.
source("library/all.R")



### COMPUTING WITH SIZE-SHAPES ###

# Here are the semi-axis lengths of two ellipsoids. Let's average them. R gives the right answer, but it's a bad answer. Why?
obstA <- c(5, 1, 2)
obstB <- c(1, 5, 2)
arithmeticMean(list(obstA, obstB))

# The easiest solution is to insist that the semi-axis lengths be recorded in a certain order (ascending or descending). Geologists do this all the time. Frequently they go on to focus on the longest semi-axis or the shortest semi-axis. Suppose that we have a bunch of ellipsoids, and here are their longest semi-axis lengths:
obstA1s <- c(0.4, 0.5, 2.1, 0.3, 0.6, 0.1, 0.2)

# Using a common R function, here is the 95% confidence interval for the mean of the longest axis. Why is this answer bad?
t.test(obstA1s)

# Statisticians run into this kind of problem a lot. If your data are positive numbers, in a context where zero and negative numbers don't even make sense, then you have to use a method that's aware of that fact. Frequently the solution is to (A) transform to the logarithms of the numbers, (B) do statistics there, and (C) transform the results back.
obstLogA1s <- log(obstA1s)
obstLogA1sTest <- t.test(obstLogA1s)
exp(obstLogA1sTest$estimate)
exp(as.numeric(obstLogA1sTest$conf.int))

# Mathematically, the fundamental reason why this tactic works is that the logarithms live in a 'vector space'. A vector space is a context where basic arithmetic operations (addition, subtraction, scaling) make sense. Hence the statistical calculations built from those operations also make sense.

# So to analyze entire size-shapes of ellipsoids, should we just take the logs of their semi-axis lengths, and work with those vectors of logs? Sorry, but no. These log-vectors are not just any 3D vectors. They have to be ordered in descending (or ascending order). Consequently the set of all log-vectors does not form a vector space. Statistics is not going to work well --- at least, not without some weird gymnastics.



### COMPUTING WITH ORIENTATIONS ###

# Load a data set of 27 anisotropy of magnetic susceptibility (AMS) ellipsoids from Cyprus. Their orientations are quite widely dispersed.
cyprusAMSData <- geoEllipsoidDataFromIRMFile("data/cyprus_AMS_groupF.tsv", doNormalize=FALSE)
ellEqualVolumePlot(cyprusAMSData$rotation, cyprusAMSData$a)

# Plot their axes: circles for short axes, triangles for intermediate axes, squares for long axes. Yep, dispersed.
ellEqualAreaPlot(cyprusAMSData$rotation, cyprusAMSData$a)

# To some degree, the axes tend to align in three dense regions of the equal-area plot. Strangely, each dense region contains a mixture of short, intermediate, and long axes. It's almost as if we mixed up the lengths of the axes, without messing up their orientations.

# Do these plots help you understand why the short, intermediate, and long axes might be getting mixed up?
ellHsuNadaiPlot(cyprusAMSData$logA)
ellHsuNadaiPlot(cyprusAMSData$logA, es=0.05)
ellHsuNadaiScalarPlot(cyprusAMSData$logA, log(sapply(cyprusAMSData$logA, ellVolume)) / 1000, es=0.05)

# If you're more of a verbal person than a visual person, here is the same information in another format.
t(simplify2array(cyprusAMSData$a))

# If an ellipsoid is clearly triaxial, meaning that its semi-axis lengths are quite different, then its orientation is well-defined (up to 4-fold line-in-plane symmetry). But if an ellipsoid is spheroidal, then its orientation is ill-defined. A circle's worth of 3D orientations, or equivalently a single 3D direction, describes the spheroid's orientation. Therefore, when a data set contains a mixture of triaxial ellipsoids and spheroids, finding a coherent way to analyze their orientations separately from their size-shapes is difficult. And spheres are even worse than other spheroids in this regard.



### THE EASIEST WAY TO OVERCOME THESE OBSTACLES ###

# In a perfect world, we would have a way of talking about ellipsoids that treated them holistically, as unified objects, rather than piecemeal. That treatment would automatically dodge the ill-definedness of orientation near spheroids. Further, the ellipsoids would live in a vector space, so that we could sensibly average them and do other statistics.

# Remarkably, this world exists. There is a way of packaging an ellipsoid's five or six degrees of freedom (depending on whether it is normalized or not) into a 5D or 6D vector. Under this system, each ellipsoid corresponds to one and only one vector, and each vector corresponds to one and only one ellipsoid. The one-to-one correspondence is quite well-behaved; for example, ellipsoids that are similar correspond to vectors that are close together.

# This system facilitates a simple workflow for doing statistics with ellipsoids:
# A. Convert the ellipsoids into 5D or 6D vectors. In our library, you usually do this by asking the ellipsoid for its $vector.
# B. Do any kind of statistics on those vectors that you want. R has hundreds of things for you to try.
# C. Convert the vectors back into ellipsoids. In our library, this is the function ellEllipsoidFromVector, although sometimes more complicated operations are required.
# D. Inspect the ellipsoids to understand the result.
# We'll see our first example in the next section.



### WHAT HAVE WE LEARNED? ###

# Doing statistics on ellipsoid shapes is difficult. Doing statistics on ellipsoid orientation is difficult unless those ellipsoids are all quite triaxial. So we're going to do statistics on entire ellipsoids in a strange setting.


