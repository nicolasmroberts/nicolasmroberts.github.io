


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# This tutorials demonstrates Kamb density contouring of orientational data. The level surfaces of Kamb density are plotted in the equal-volume rotation plot. Like all tutorials in tutorialsC, this tutorial requires compilation of the C part of our R library.

# The data themselves are slickenside orientations gathered from maps of central California (Woodring et al., 1940).

# Warning: It is not easy for R to stop C code while it is running. Pressing the Stop button in RStudio may not immediately stop the program. Eventually an interface may appear, giving you the option of killing R entirely. So activate a C routine only if you're sure that you want to.



### PRELIMINARY WORK ###

# If you have just restarted R, then set the working directory to the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc. Execute the following lines of code to load our custom library and its dependencies, including the C part of the library.
source("library/all.R")
source("libraryC/all.R")

# Load central California slickenside orientations. The rake describes the movement of the hanging wall. We apply geoPoleVorticityFromPoleHanging to get the orientations into a better-behaved format. Here are the orientations in equal-area and equal-volume.
slickData <- geoDataFromFile("data/cali_all_kh_faults_latlong.tsv")
slickRots <- lapply(slickData$rotation, geoPoleVorticityFromPoleHanging)
length(slickRots)
rotEqualAreaTickPlot(slickRots)
oriEqualVolumePlot(slickRots, group=oriRayInPlaneGroup)



### 6-SIGMA KAMB LEVEL SURFACE ###

# By default, the 6-sigma Kamb level surface is plotted. Notice that the surface ends abruptly at the boundary of the plot, but always continues seamlessly on the other side of the plot. Also notice that the surface consists of two symmetric copies, due to the two-fold symmetry inherent in this data type.
oricKambPlot(slickRots, group=oriRayInPlaneGroup)

# In that plot, you might see a couple of holes. They arise when the non-adaptive mesh, that initially samples the density of the data, is not fine enough. To fix the holes, increase numNonAdapt from its default of 3. But don't increase it too crazily, because each increment causes the time and memory requirements of the algorithm to increase by a factor of 8.
oricKambPlot(slickRots, group=oriRayInPlaneGroup, numNonAdapt=4)

# You can also fiddle with numAdapt, which controls the fineness of the level surface after the initial, non-adaptive mesh. Like numNonAdapt, it defaults to 3, and increasing it can drasticaly increase the time and memory requirements. Here we decrease it, to make the plot cruder and faster.
oricKambPlot(slickRots, group=oriRayInPlaneGroup, numNonAdapt=4, numAdapt=0)

# To get an idea for how the algorithm works, try starting at really low-quality and increasing from there.
oricKambPlot(slickRots, group=oriRayInPlaneGroup, numNonAdapt=1, numAdapt=0)
oricKambPlot(slickRots, group=oriRayInPlaneGroup, numNonAdapt=2, numAdapt=0)
oricKambPlot(slickRots, group=oriRayInPlaneGroup, numNonAdapt=3, numAdapt=0)



### OTHER DENSITY LEVELS ###

# Here is the Kamb 3-sigma density level surface. Notice the hole.
oricKambPlot(slickRots, group=oriRayInPlaneGroup, multiple=3)

# To remove the hole, increase numNonAdapt by 2. To avoid taking too much time, decrease numAdapt by 2.
oricKambPlot(slickRots, group=oriRayInPlaneGroup, multiple=3, numNonAdapt=5, numAdapt=1)

# Just for fun, here is the 18-sigma level surface.
oricKambPlot(slickRots, group=oriRayInPlaneGroup, multiple=18, numNonAdapt=4, numAdapt=2)



### FOR PUBLICATION ###

# Here we tweak the plot of the 3-sigma level surface. After making the plot, we maximize the window. Then we call afterMaximizingWindow to take a screen shot and save it to the file 'figKamb.png'. One of the figures in the readme.pdf was made exactly like this. To learn more about saving plots, see the bonus tutorial publicationPlots.R.
oricKambPlot(slickRots, group=oriRayInPlaneGroup, multiple=3, numNonAdapt=5, numAdapt=1,
             boundaryAlpha=0, backgroundColor="white", colors=c("black"))
afterMaximizingWindow(leftName="figKamb.png")



### WHAT HAVE WE LEARNED? ###

# You can plot one Kamb density level surface at a time. If you want high quality, you might have to tweak some parameters.


