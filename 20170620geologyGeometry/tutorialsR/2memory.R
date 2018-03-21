


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# In this tutorial we learn how to store the results of computations in memory. This is a basic and useful feature. It is also a commmon source of programming errors, when the contents of memory aren't what we think they are.



### STORING VALUES IN MEMORY ###

# Modern calculators do much more than isolated simple calculations. You can store the result of a calculation in memory, for later use. This may sound like a small feature, but it's an important stepping stone to programming.

# For example, maybe you prefer to work with angles in degrees, not radians. Here's the cleanest way to proceed.
degree <- 2 * pi / 360
sin(90 * degree)
sin(30 * degree)

# When you use the '<-' operator, R evaluates the stuff on the right side, puts its value into the memory location named on the left side, and DOESN'T show you the answer. But if you want to see what got stored in memory, just ask R.
degree

# For another example, suppose that you run a small business. You want to compute its quarterly profit. (And you've been studying R for only a few minutes.) You might enter three months of revenues and costs like this.
revenue <- 13000 + 15000 + 14000
cost <- 16000 + 12000 + 5000

# To compute your quarterly profit, you might do this.
revenue - cost

# Or, if you think that you might want to recall that profit later, then do this instead.
profit <- revenue - cost
profit



### MEMORY IS TEMPORARY ###

# When you quit RStudio, RStudio tells the underlying R software to stop running too. And when that happens, everything that you've stored in memory is lost. (There is actually an option for saving R's memory, but I recommend that you don't grow dependent on it.) So keep these practices in mind:

# If you're just trying out a little R command here and there, then type it directly into the R interpreter (the RStudio Console).

# If you want to save R code for later use, then type it into the editor. And save that R code to a file, often.



### MEMORY VS. PROGRAM ###

# Run the following lines of code.
x <- 11
x + 4

# Now read, but DO NOT RUN, the following line of code.
x <- 3

# Now try to predict what the following line of code will do. Then run the line of code, to check your intuition.
x + 4

# Does the result make sense? Even though we're looking at a line of code that says 'x <- 3', the value of x in memory is still 11, because we never ran the 'x <- 3' line. There is a disagreement between the text of the program (what we see, as users) and the contents of memory (what R uses, in its calcuations).

# Such disagreements are a REALLY common source of programming errors. You might not believe it yet, because the preceding example was so simple. The errors start cropping up when the code gets more complicated. Try this.
x <- 11
y <- x^2
y + 4

# Now run these lines, and observe what y + 4 is.
x <- 3
y + 4

# Because y was set to x^2 and x was set to 3, you might think that y + 4 would be 13. But the command 'y <- x^2' doesn't establish an ongoing relationship between x and y. It just sets the value of x, based on the value of x at that time. If the value of x changes later, the value of y does not automatically change with it. In other words, just because y equaled x^2 at one point, that doesn't mean that y equals x^2 now.

# Ultimately, you have to keep careful track of which commands have been executed, and in what order. In large programs this becomes impossible. Software developers have a number of strategies for managing such complexity.



### WHAT HAVE WE LEARNED? ###

# Storing values in memory is a double-edged sword: very useful, but also a common source of errors.


