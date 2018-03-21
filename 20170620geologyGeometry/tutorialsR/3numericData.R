


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# This tutorial describes three ways of organizing data: vectors, lists, and data frames.



### PRELIMINARY WORK ###

# If you have just restarted R, then set the working directory to the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc. Execute the following line of code to load our custom library and its dependencies.
source("library/all.R")



### VECTORS ###

# In the preceding tutorial, we entered some financial data like this.
revenue <- 13000 + 15000 + 14000
cost <- 16000 + 12000 + 5000

# A real R programmer would never manage data like that. First of all, she would store the data in memory before starting computations on it. Here we store the data in a 'vector', meaning roughly a sequence of numbers constructed using the 'c' operator.
revenues <- c(13000, 15000, 14000)
costs <- c(16000, 12000, 5000)
revenues
costs

# Then she might use the 'sum' function to perform the addition.
revenue <- sum(revenues)
cost <- sum(costs)
revenue
cost
revenue - cost

# Or she might want to know the profit in each month. So she would subtract the costs vector from the revenues vector, to get a vector of profits.
profits <- revenues - costs
profits
sum(profits)

# Or she might plot revenue vs. cost, like this. The plot appears in RStudio's Plots pane.
plot(x=costs, y=revenues)

# Or she might want to inspect just the third revenue.
revenues[[3]]

# Or she might realize that the second revenue has an error. It should actually be 16000. So she changes it, inspects the updated revenues, and recomputes profits.
revenues[[2]] <- 16000
revenues
revenues - costs



### DATA FRAMES ###

# Why is this two-stage process --- load the data, then compute on the data --- a useful way of thinking? Because usually we don't type in data at all. Usually we load the data from a file. And then we do LOTS of computations on it.

# The data file 'synthFinanicial.csv' contains 13 months of revenues and costs. Feel free to open it in a spreadsheet or text editor, to see what it looks like. The following line of code loads that data file into R.
financialData <- read.table("data/synthFinancial.csv", sep=",", header=TRUE)
financialData

# Notice that financialData is an entire table of data --- what R calls a 'data frame'. You can pick off the revenue and cost columns like this. Notice that they are vectors.
financialData$revenue
financialData$cost

# Then you might compute the profit like this.
revenue <- sum(financialData$revenue)
cost <- sum(financialData$cost)
revenue - cost

# Or you might do it all in one line, like this.
sum(financialData$revenue) - sum(financialData$cost)

# Or you might plot revenue vs. cost.
plot(x=financialData$cost, y=financialData$revenue)

# Or you might decide that you want the profits appended to the data frame. Here's how.
financialData$profit <- financialData$revenue - financialData$cost
financialData



### DATA FRAMES, CONTINUED ###

# Maybe you want to plot costs over time. This is a little tricky, because there are two kinds of time: years and months. So we're going to join them together into a single kind of time: decimal years. We slap that onto the data frame as a new column.
financialData$decYear <- financialData$year + (financialData$month - 1) / 12
financialData

# If that code went by too quickly, then run these lines of code one at a time, to see how the pieces fit together.
financialData$month
financialData$month - 1
(financialData$month - 1) / 12
financialData$year + (financialData$month - 1) / 12

# Do those steps make sense? Anyway, now we can plot revenues or costs against decimal years.
plot(x=financialData$decYear, y=financialData$cost)
plot(x=financialData$decYear, y=financialData$revenue)



### LISTS ###

# Here is a more geologic example: 18 dikes from a certain field site in Cyprus.
site230 <- geoDataFromFile("data/cyprusDikesSite230.tsv")
site230

# Here are the strikes and dips (in degrees).
site230$strikeDeg
site230$dipDeg

# As we learn in the directional statistics tutorials, you don't usually compute with strikes and dips. Instead you convert them into pole vectors and compute on those. In our case, the pole vectors have been computed automatically. (The geoDataFromFile command is essentially read.table combined with a bunch of automatic computations like this.)
site230$pole

# The pole column of the site230 data frame is a new kind of beast: a 'list'. It is a list of 18 elements. Each element is not simply a number, but rather a vector. You can access the elements in the list using the same syntax as for accessing the elements in a vector.
site230$pole[[15]]

# And you can access the individual numbers in those vectors as you would access the numbers in any vector.
site230$pole[[15]][[2]]

# When using our R library, or indeed anybody's R code, you need to know exactly how that R code wants to be presented with data. For example, suppose that we want to make an equal-area plot of the dike poles, using the lineEqualAreaPlot command. You might try telling lineEqualAreaPlot the strikes and dips. But this causes an error.
lineEqualAreaPlot(site230$strikeDeg, site230$dipDeg)

# You need to know that lineEqualAreaPlot expects a list. And that list is a list of 3D vectors. (This is all explained in the documentation.) Fortunately, the pole vectors for the dikes are ready to go.
lineEqualAreaPlot(site230$pole)



### WHAT HAVE WE LEARNED? ###

# For our purposes, there are three important structures for storing numerical data: vectors, lists, and data frames. A vector is appropriate for storing a sequence of numbers. A list is appropriate for storing a sequence of more complicated objects, such as vectors. A data frame is appropriate for storing a table of data. Each column in the table is either a vector (if the entries are simple) or a list (if the entries are more complicated).


