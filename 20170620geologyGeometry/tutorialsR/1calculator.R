


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# This tutorial begins our gentle introduction to R. We see that R is essentially a calculator, and an R program is nothing more than a sequence of commands to give to R.



### GETTING STARTED WITH RSTUDIO ###

# If you're viewing this file in anything but RStudio, then you're making life unnecessarily hard. Launch the RStudio application and open this file using its 'Open File...' menu command.

# RStudio presents a single large window divided into 'panes'. One pane shows plots, another displays online help, another lets you manage your files, etc. In the preferences you can configure how the panes are laid out. Here are the two most important panes:

# In the upper left corner of the RStudio window is the editor pane. It's where you type your R programs, so that you can save them in files for later use.

# In the bottom left corner of the RStudio window is the Console pane. It's where you enter commands into R itself. Usually you see the results here, too. (If your command makes a plot, then it might show up in the Plots pane or a separate window. If your command asks for help, then it might show up in the Help pane.)



### GLORIFIED CALCULATOR ###

# R is essentially a glorified calculator. Type (or copy and paste) the following line of code into RStudio's Console pane, and press the Return key. R evaluates your input and reports an answer.
3 + 2

# Similarly, run each of the following lines of code in R. You can probably guess what they mean. Feel free to try other numbers and operations. Experiment!
(-5)^2
sqrt(16 - 7)

# In RStudio, there is another way to run a line of code. Click anywhere on the line of code that you want to run, and then click the Run button in the top right corner of the editor pane. RStudio sends the line of code to R. RStudio also advances the cursor to the next line of code, so that you can immediately run that if you want. You can also run lines of code from an RStudio menu, which has a keyboard shortcut (e.g., Command-Return on Mac OS X).
log(17)
log(1 + 2^2 + 4 * 2^3)

# As in almost all programming languages and computer systems, angles are in radians by default.
sin(pi / 2)
atan(4 / 3)
atan(-4 / -3)
atan2(-4, -3)

# To view documentation on any R function, enter that function's name into the search bar of RStudio's Help pane. Or use a question mark as in the following line of code. Either way, help appears in the Help pane. You can then click on the Help pane's zoom button to get a larger view of that help.
?atan2

# (This help feature doesn't yet work for our geologyGeometry library. So instead open reference.txt and search within that file for the function name.)



### WORKING DIRECTORY ###

# Your computer's hard drive is a hierarchy of directories (folders), subdirectories, subdirectories of those, etc. At any given time, R has a concept of 'working directory' that it uses for accessing files. Whenever you read a file into R (for example, a data file or a file of R code), that file's location is specified relative to the working directory. Whenever you write a file from R (for example, a file of results), that file's location is specified relative to the working directory.

# In RStudio, your current working directory is always printed next to the word 'Console' in the Console pane. You can change your working directory in the 'Set Working Directory' menu item. Alternatively, you can set it by navigating to a directory in the Files pane, pressing the 'More' button there, and then choosing 'Set As Working Directory'.

# Right now, set your working directory to the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc.



### AN R PROGRAM ###

# A program is nothing more or less than a sequence of statements stored in a file. You can run those statements one at a time, by hand. Or you can ask R to run an entire file's worth of statements --- that is, a program --- all at once.

# For example, this next line of code runs the entire third R tutorial, all at once. If R can't file the file, then your working directory isn't set correctly. What should happen is: Some plots appear in RStudio's Plots pane, and an error appears (because there is an intentional error in the program, on line 129).
source("tutorialsR/3numericData.R")

# As your R programs get longer and longer, they tend to get more and more confusing. You can alleviate this confusion by adding comments. Any line beginning with '#' is a comment. It is for human consumption only. It is ignored by the R interpreter. For example, try running this next line. Nothing happens, unless you remove the '#'.
#125^(1 / 3)

# Comments should not just replicate the code. They shouldn't give a 'play-by-play account' of what the code is doing. Instead, they should describe the programmer's intent: What the goal of the code is, or its role in the larger problem.



### WHAT HAVE WE LEARNED? ###

# At its heart, R is essentially a big calculator. It does mathematical computations. But over the next few tutorials we will see how to build larger and larger structures out of those computations, until they become useful programs.


