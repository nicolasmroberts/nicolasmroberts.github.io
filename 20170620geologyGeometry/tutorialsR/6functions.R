


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# In previous tutorials we have learned three crucial elements of programming: saving values in memory, branches (if-else statements), and loops (mainly through sapply and lapply). In this tutorial we learn a final crucial element: user-defined functions. This is a big topic, and this tutorial just scratches the surface, but we hope that it de-mystifies R a bit, at least.



### PRELIMINARY WORK ###

# If you have just restarted R, then set the working directory to the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc. Execute the following line of code to load our custom library and its dependencies.
source("library/all.R")



### WHAT IS A FUNCTION? ###

# R comes with many functions built in. Some are for math, some are for statistics, some are for managing data, etc. Here are some examples: sine, sample mean and standard deviation, and loading a data frame from a file.
sin(-pi / 4)
myData <- c(3, 1, -7, 2, 4, 6)
mean(myData)
sd(myData)
read.table("data/synthFinancial.csv", sep=",", header=TRUE)

# And here are two functions from our library: geoDataFromFile and lineEqualAreaPlot.
wiszData <- geoDataFromFile("data/wiszFollins.tsv")
lineEqualAreaPlot(wiszData$pole)

# But what is a function? A function is essentially a block of code that (A) takes in some inputs, (B) does some computation on them, and then (C) produces some outputs. In most functions (including sin, sd, read.table, and geoDataFromFile), the outputs are 'returned', so you can store them and use them in subsequent computations. For example, sin returns a number between -1 and 1.
x <- sin(-pi / 4)
x
x + 5 * x^3
sin(pi / 12)^2 - 6

# In some other cases, the output of the function is a 'side effect'. For example, the lineEqualAreaPlot function returns nothing (which is called NULL in R), but along the way it also draws a plot.
value <- lineEqualAreaPlot(wiszData$pole)
value

# Functions are useful because they let you package a bunch of work into a unit that can be used again and again. Once you've figured out how to do a task in R, and you find yourself doing that task over and over again, you will want to write a function to package that work.



### LOWER-HEMISPHERE VECTORS ###

# Recall this code snippet from the R tutorial 4otherData.R? The first line randomly generates a 3D vector. The second line scales the vector to have length 1. The third line asks R what the vector is. The fourth and fifth lines negate the vector, if necessary, to make it lower-hemisphere. The sixth line again asks R what the vector is.
v <- rnorm(3)
v <- v / sqrt(sum(v * v))
v
if (v[[3]] > 0)
  v <- -v
v

# Flipping a vector to the lower hemisphere is something that a geologist might want to do often --- for example, as part of making lower-hemisphere plots. So this chunk of code is a good candidate for a function. Here's my attempt at writing one.
# * The first line says, 'We are defining a function called lowerHemisphere. It takes a single input, which we call v. And whatever that v is, the function does the following...'.
# * The second line says, 'Scale v so that it has length 1, and store that answer in a new memory location called 'scaled'.
# * The third line tests whether scaled has positive z-component. The fourth line negates it if so.
# * Now here's the key. The return value of a function is defined to be whatever the last statement in the function body evaluates to. In lowerHemisphere, the last statement is simply 'scaled'. So this function returns whatever scaled equals at the end of the computation.
lowerHemisphere <- function(v) {
  scaled <- v / sqrt(sum(v * v))
  if (scaled[[3]] > 0)
    scaled <- -scaled
  scaled
}

# Once the function is defined (loaded into memory), you can use it as often as you like. Here's how it fits into the example above.
v <- rnorm(3)
lowerHemisphere(v)

# Is it important to name a vector 'v' before passing it to the function? Maybe because the vector is called 'v' inside the function? No. The two uses of 'v' have nothing to do with each other. In fact, you don't have to name the vector outside the function at all. Here are a few examples.
myVector <- rnorm(3)
lowerHemisphere(myVector)
lowerHemisphere(rnorm(3))
lowerHemisphere(wiszData$pole[[1]])

# Similarly, what are the rules for the vector called 'scaled' inside the function? If I have a vector called 'scaled' outside the function, then will this function accidentally change its value? No.
scaled <- c(1, 2, 3)
lowerHemisphere(scaled)
scaled

# The memory locations named inside a function are 'isolated'. They have nothing to do with memory locations outside the function (at least in simple examples like this). That's good. It means that you can write a function once and use it in a variety of contexts, without worrying about trampling on other memory locations.

# In the preceding tutorial we learned about performing loops using sapply and lapply. You can use these operations on your own user-defined functions. Roughly speaking, there are three cases:
# * If your function returns a scalar, then sapply it.
# * If your function returns NULL (like lineEqualAreaPlot), then you probably don't want to sapply or lapply it, because you'll just get a vector or list of NULLs.
# * If your function returns something other than a scalar or NULL, then lapply it. Here's an example.
vs <- replicate(10, rnorm(3), simplify=FALSE)
vs
lapply(vs, lowerHemisphere)



### BUT WHY DID WE DESIGN IT THAT WAY? ###

# Now that you understand how to define and use lowerHemisphere, we need to discuss a harder question: Why have we designed it that way?

# For starters, the original block of code included a line to randomly generate v. Maybe the function should include that part? In which case the function doesn't need to take any inputs at all?
lowerHemisphereNEW <- function() {
  v <- rnorm(3)
  scaled <- v / sqrt(sum(v * v))
  if (scaled[[3]] > 0)
    scaled <- -scaled
  scaled
}

# Here's how that function works. Run it a few times.
lowerHemisphereNEW()

# Is lowerHemisphereNEW okay? Sure. It is useful? Not very. The point of lowerHemisphere is to help us make lower-hemisphere plots of data, and stuff like that. Usually I get the data from geology. So lowerHemisphere shouldn't begin by making up v randomly. My design is: lowerHemisphere should take one vector as input, and it should return one vector as output.

# Once you've made all of your design choices, you also get to choose how the guts of the function are written. Here's a new version of lowerHemisphereNEW, that works exactly like lowerHemisphere, even though the code is a little different. In my opinion, it is equally as good as lowerHemisphere.
lowerHemisphereNEW <- function(v) {
  scaled <- v / sqrt(sum(v * v))
  if (scaled[[3]] > 0)
    -scaled
  else
    scaled
}

# Here's yet another version. It is mathematically equivalent to lowerHemisphere, but much worse in its code. It seems to be computing the same stuff over and over again. That repetition makes the code unnecessarily hard to understand. (It may or may not also make the code slow to run.)
lowerHemisphereNEW <- function(v) {
  if ((v / sqrt(sum(v * v)))[[3]] > 0)
    -(v / sqrt(sum(v * v)))
  else
    (v / sqrt(sum(v * v)))
}

# Do we have to name this function 'lowerHemisphere'? No. We could name it pretty much anything we want, such as 'jimmy'. But 'jimmy' is a dumb name. Don't use 'jimmy'. Use something descriptive.
jimmy <- function(v) {
  scaled <- v / sqrt(sum(v * v))
  if (scaled[[3]] > 0)
    scaled <- -scaled
  scaled
}
jimmy(wiszData$pole[[1]])
jimmy(wiszData$pole[[2]])

# In fact, our library offers a function much like this one, called 'lower'. If you ask R what lower is, it will show you the code for it.
lower

# There are a few cosmetic differences between lowerHemisphere and lower, but there is also one big, important difference. What?

# The function lowerHemisphere scales its input to be length 1. The function lower doesn't, because it assumes that its input v already has length 1. Which definition is right? Neither and both. It's a design choice. lowerHemisphere is safer, because it works even when v is not length-1 already. lower is faster, because it skips the rescaling operation. In our library, lower is mostly used behind the scenes, not by the user. And usually we want to lower vectors that are known to have length 1 already. And we might have to use lower thousands of times in a single equal-area plot. So we favor speed over safety.

# Compared to a professional computer programmer, we're doing very small tasks here. But you can already see that there is some judgment, creativity, and 'taste' in programming.



### ANGLE BETWEEN VECTORS ###

# Here's a function to compute the angle, in radians, between two vectors. It uses our geologyGeometry library's function for the dot product of two vectors.
angleBetweenVectors <- function(v, w) {
  vScaled <- v / sqrt(dot(v, v))
  wScaled <- w / sqrt(dot(w, w))
  cosAngle <- dot(vScaled, wScaled)
  acos(cosAngle)
}
angleBetweenVectors(c(1, 2), c(3, 6))
angleBetweenVectors(c(-3, 2, 1), c(1, 4, 4))
angleBetweenVectors(c(1, 0, -1, 2), c(3, 1, 1, -2))

# Have you noticed that we seem to be repeating a certain task over and over again: scaling vectors to have length 1? We've done it once in lowerHemisphere and twice in angleBetweenVectors. That task is starting to look like a good candidate for packaging into a function. In fact, our library has that function, under the name 'rayNormalized'.
rayNormalized

# So here is a rewrite of angleBetweenVectors, using rayNormalized.
angleBetweenVectors <- function(v, w) {
  vScaled <- rayNormalized(v)
  wScaled <- rayNormalized(w)
  cosAngle <- dot(vScaled, wScaled)
  acos(cosAngle)
}

# And here is another rewrite. It works exactly the same way, using the same amount of computational time, but the code is much more succinct. I like it. But in general don't try to cram too much work onto a single line, if it makes your code less readable.
angleBetweenVectors <- function(v, w) {
  acos(dot(rayNormalized(v), rayNormalized(w)))
}

# Now can you also understand why the function lower doesn't scale its input automatically? In the cases where we really want its input vector pre-scaled, we can easily do that using rayNormalized. In other words, there is a clear division of labor between the functions lower and rayNormalized. Each one does its specific task well, with no overlap. That's a design choice.
lower(rayNormalized(rnorm(3)))



### TREND AND PLUNGE FROM STRIKE AND DIP ###

# Here's a function that takes in a strike and dip (in degrees), as a 2D vector. And it returns trend and plunge (in degrees), as a 2D vector. The relationship between the two is that the strike and dip describe a plane, and the trend and plunge describe the pole to that plane.
trendPlungeDegFromStrikeDipDeg <- function(strDipDeg) {
  trendDeg <- (strDipDeg[[1]] - 90) %% 360
  plungeDeg <- 90 - strDipDeg[[2]]
  c(trendDeg, plungeDeg)
}
trendPlungeDegFromStrikeDipDeg(c(30, 80))
trendPlungeDegFromStrikeDipDeg(c(242, 11))

# You may be wondering what the '%%' operator in the function means. It is remainder division. It forces our angles into the interval [0, 360). These examples should give you the idea.
215 %% 360
007 %% 360
378 %% 360
-10 %% 360

# You may be wondering why 'Deg' appears so much in the function. Remember that R works in radians by default (because math works better in radians). And you don't want to confuse radians with degrees, or your program will have errors, which might be difficult to track down. So whenever a function or a memory location has anything to do with degrees, I label it with 'Deg'. It's like a skull and crossbones on a bottle of poison.

# Why didn't I design the function to take strike and dip as separate inputs, like this?
trendPlungeDegFromStrikeDipDegNEW <- function(strDeg, dipDeg) {
  trendDeg <- (strDeg - 90) %% 360
  plungeDeg <- 90 - dipDeg
  c(trendDeg, plungeDeg)
}
trendPlungeDegFromStrikeDipDeg(c(30, 80))
trendPlungeDegFromStrikeDipDegNEW(30, 80)

# That would be okay too. Define your functions however is best for the rest of your workflow. Unfortunately, you might not be able to judge that, until you've written a bunch of other code that uses the functions. Experience helps.

# In fact, because the trend and plunge calculations have nothing to do with each other, you could even design two separate functions for them, if that's what's useful to you.
trendDegFromStrikeDeg <- function(strDeg) {
  (strDeg - 90) %% 360
}
plungeDegFromDipDeg <- function(dipDeg) {
  90 - dipDeg
}
trendPlungeDegFromStrikeDipDeg(c(30, 80))
trendPlungeDegFromStrikeDipDegNEW(30, 80)
trendDegFromStrikeDeg(30)
plungeDegFromDipDeg(80)

# Is any version of this function in our library? No, because we don't find this task particularly common or useful. To do this task with our library requires three steps: convert from strike and dip to Cartesian coordinates, make sure it's lower-hemisphere, and then convert Cartesian coordinates to trend and plunge.
geoTrendPlungeDegFromCartesian(lower(geoCartesianFromStrikeDipDeg(c(30, 80))))

# Why is our library designed in this way? When presented with geologic data, we immediately convert to Cartesian coordinates. Then we compute in those coordinates --- possibly for a long time, doing various complicated tasks. Then, once we have the answer, we convert back to a geological format.

# For example, both of the following lines of code show us the first line of the file 'data/wiszFollins.tsv'. The built-in R function read.table just reads the file. Our geoDataFromFile function reads the file and immediately computes Cartesian coordinates. (geoDataFromFile also appends 'Deg' to 'strike', 'dip', 'trend', and 'plunge', because of the poison.)
read.table("data/wiszFollins.tsv", sep="\t", header=TRUE)[1,]
geoDataFromFile("data/wiszFollins.tsv")[1,]



### A LARGER EXAMPLE ###

# There's no limit to how big functions can get. For example, let's inspect lineEqualAreaPlot's code.
lineEqualAreaPlot

# Notice that a function f is defined inside lineEqualAreaPlot. Is that allowed? Yes. Because of the importance of sapply and lapply, defining temporary functions like this is pretty common in R programs.

# lineEqualAreaPlot also uses a function called 'rayCurvesUpperLower', to break curves up at the boundary of the plot. Let's see what that function looks like. (It's pretty ugly, actually. I'm not proud of it.)
rayCurvesUpperLower

# lineEqualAreaPlot ends by calling plotEqualArea. In fact, everything up to that point is really a preamble, to get the data arranged for plotEqualArea, which actually does the drawing.
plotEqualArea

# So you can see how arbitrarily complicated features can be constructed from smaller pieces. To make it all work, you just need to make sure that each piece fulfills its role in the design.



### WHAT HAVE WE LEARNED? ###

# Functions are a way of packaging code for re-use. A lot of judgment goes into designing exactly how the function should work, to make it as useful as possible to the rest of your code.


