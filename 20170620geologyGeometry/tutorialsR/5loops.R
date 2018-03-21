


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# In previous tutorials we have learned two crucial elements of programming: saving values in memory, and branches (if-else statements). In this tutorial we learn a third crucial element: loops. And R handles loops differently from some other programming languages.



### PRELIMINARY WORK ###

# If you have just restarted R, then set the working directory to the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc. Execute the following line of code to load our custom library and its dependencies.
source("library/all.R")



### FOR ###

# Computers excel at doing mindless tasks over and over again at high speed. Therefore most computer programming languages have a notion of 'loop', that lets you repeat a calculation many times. For example, here are the square roots of the first 65 natural numbers.
roots <- c()
for (i in 0:64)
  roots <- c(roots, sqrt(i))
roots

# As an aside, if you wanted to see how the vector of roots grows over time, you could add a print statement.
roots <- c()
for (i in 0:64) {
  print(roots)
  roots <- c(roots, sqrt(i))
}
roots



### SAPPLY ###

# I'm not going to explain that preceding example thoroughly, because that kind of loop is actually pretty rare in R programs. Here is the preferred R way to do that computation.
roots <- sapply(0:64, sqrt)
roots

# Let's break it down. First, what is '0:64'? It's a quick way to make a vector of numbers from 0 to 64. That's not hard.
0:64

# Now what does sapply mean? It means, roughly: 'Apply this function to each element of this vector or list.' In this example, we apply sqrt to each element of the vector 0:64, and get a vector of the square roots.
sapply(0:64, sqrt)

# Here's another example, where we apply sin to each element of a vector.
myAngles <- c(0, pi / 6, pi / 4, pi / 3, pi / 2)
sapply(myAngles, sin)

# Here's another example, where we apply det to each matrix in a list of 23 matrices. (These are rotation matrices, so they all have determinant 1.)
wiszData <- geoDataFromFile("data/wiszFollins.tsv")
wiszData$rotation
sapply(wiszData$rotation, det)



### LAPPLY ###

# Here's a similar example, where we apply t (transposition) to the same list of 23 matrices. The result is a crazy 9x23 matrix.
sapply(wiszData$rotation, t)

# Each column of that matrix is one of the original 23 matrices, transposed, with its 9 elements listed in a particular order. That's not a crime, but I would rather have a list of 23 matrices. I can get what I want by using lapply instead of sapply. Briefly, lapply means: 'Apply this function to each element of this vector or list, and give me the results as a list.'
lapply(wiszData$rotation, t)

# The difference between lapply and sapply is: sapply does lapply to get a list of results, and then simplifies that list down to a vector or matrix. The following direct comparison should make the distinction clear. In this example, I favor sapply, because the elements of the output list are just numbers.
sapply(wiszData$rotation, det)
lapply(wiszData$rotation, det)

# Here's another comparison, where I favor sapply.
sapply(0:64, sqrt)
lapply(0:64, sqrt)

# And here's another comparison, where I favor lapply, because the list elements are not merely numbers but vectors.
sapply(wiszData$direction, geoTrendPlungeDegFromCartesian)
lapply(wiszData$direction, geoTrendPlungeDegFromCartesian)



### ANONYMOUS FUNCTIONS ###

# So what if I want to get the third component of each vector in wiszData$direction? Kind of like this...
wiszData$direction[[1]][[3]]
wiszData$direction[[2]][[3]]
wiszData$direction[[3]][[3]]

# ...and so on. You can do it as follows. You write a little function 'function(v) v[[3]]' to pick off the third component of a vector v. Then you sapply that function to the list of vectors.
sapply(wiszData$direction, function(v) v[[3]])

# If you find this 'function' thing confusing, don't worry too much. The next tutorial really teaches functions in detail. Here we just look at a few simple examples, to get an idea of how functions are used in sapply and lapply.

# For another example, suppose that you want to get the length of each vector in wiszData$direction. The length of a vector v is sqrt(sum(v * v)). For example:
v <- c(0, 4, -3)
v * v
sum(v * v)
sqrt(sum(v * v))

# So you write a little function to compute the length of any vector v, and you sapply that function. (In this example, every vector has length 1, because these are normalized direction vectors.)
sapply(wiszData$direction, function(v) sqrt(sum(v * v)))

# Or maybe you want to compute the trend of each vector in wiszData$direction. You can do this by getting the trend and plunge, and then picking the first of those two numbers. So you write a little function and sapply it.
sapply(wiszData$direction, function(v) geoTrendPlungeDegFromCartesian(v)[[1]])

# A question for you: How would you alter the preceding example, to return the plunge of each vector? To return the strike of the plane perpendicular to that vector? To return the dip of that plane?

# Or maybe you want to compute the trace of each rotation matrix in wiszData$rotation?
sapply(wiszData$rotation, function(r) (r[[1, 1]] + r[[2, 2]] + r[[3, 3]]))

# Or maybe you want to pick off the third row of each rotation matrix in wiszData$rotation? Use lapply, not sapply.
lapply(wiszData$rotation, function(r) r[3,])

# The functions in this section are called 'anonymous', because we haven't bothered to give them names, because we're using them only once and then forgetting about them. In the next tutorial, we will write bigger functions and store them in memory, so that we can use them again and again.



### WHAT HAVE WE LEARNED? ###

# Loops in R are usually accomplished through sapply or lapply (or other, similar functions). sapply returns its results in a vector (or matrix), while lapply returns them in a list.


