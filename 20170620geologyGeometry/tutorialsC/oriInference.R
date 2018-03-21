


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# This tutorial demonstrates Markov chain Monte Carlo (MCMC) simulation for orientational data. It is not self-contained. Rather, it is intended to be studied immediately after tutorial 6inference.R in tutorialOrientations. Like all tutorials in tutorialsC, this tutorial requires compilation of the C part of our R library.

# The data themselves are foliation-lineation orientations from the western Idaho shear zone (Giorgis and Tikoff, 2004).

# Warning: It is not easy for R to stop C code while it is running. Pressing the Stop button in RStudio may not immediately stop the program. Eventually an interface may appear, giving you the option of killing R entirely. So activate a C routine only if you're sure that you want to.



### PRELIMINARY WORK ###

# This tutorial is not self-contained. You are expected to have run tutorial 6inference.R from tutorialsOrientations immediately before this tutorial. So you should already have the data and most of the required libraries loaded. Execute the following line of code to load the C part of our libraries.
source("libraryC/all.R")



### CREDIBLE REGION ###

# Remember that the sample size is n = 23 and the Fisher concentration tensor K-hat has eigenvalues 33, 11, 0.000003. Based on the numerical experiments reported in our orientation statistics paper, we proceed by Markov chain Monte Carlo simulation. The number of MCMC samples collected is 100 * 10,000 = 1,000,000.
wiszMCMC <- oricWrappedTrivariateNormalMCMCInference(wiszData$rotation, group=oriLineInPlaneGroup, numCollection=100)

# Although the MCMC credible region is computed based on all 1,000,000 samples, only 10,000 of those samples are passed back to us for inspection. We're supposed to check that they form a tight, ellipsoidally shaped cloud. Yep.
oriEqualAnglePlot(wiszMCMC$ms, group=oriLineInPlaneGroup, simplePoints=TRUE)
rotEqualVolumePlot(wiszMCMC$ms, simplePoints=TRUE)

# Here are those samples in equal-area. The foliation poles are nearly horizontal and the lineation directions are nearly vertical.
lineEqualAreaPlot(c(lapply(wiszMCMC$ms, function(r) r[1,]), lapply(wiszMCMC$ms, function(r) r[2,])), shapes=c("."))

# Here's the ellipsoidal 95% credible region, containing the middle 95% of the MCMC samples. It is much closer to spherical than the bootstrap confidence region is. It is also larger.
rotEllipsoidPlot(wiszMCMC$ms, wiszMCMC$mBar, wiszMCMC$leftCovarInv, wiszMCMC$q095^2, simplePoints=TRUE, numNonAdapt=5)



### HYPOTHESIS TESTS ###

# This plot shows the predicted foliation-lineations missing the 95% credible region.
wiszPredictions <- oriNearestRepresentatives(wiszPredictions, wiszMCMC$mBar, group=oriLineInPlaneGroup)
rotEllipsoidPlot(wiszPredictions, wiszMCMC$mBar, wiszMCMC$leftCovarInv, wiszMCMC$q095^2, numNonAdapt=4, simplePoints=TRUE)

# In other words, if we did a hypothesis test with any one of these predictions as the hypothesized mean, then that hypothesis would be rejected with a p-value less than 0.05.  More precisely, here is the range of p-values attained. (The particular numbers you see will depend on exactly how your MCMC went.)
range(sapply(wiszPredictions, wiszMCMC$pvalue))

# So we reject the entire proposed class of deformations, as an explanation for these data.



### WHAT HAVE WE LEARNED? ###

# In this problem, MCMC simulation produces results similar to, but not identical to, those of the bootstrapping simulation.


