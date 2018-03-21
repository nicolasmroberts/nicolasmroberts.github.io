


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# This tutorial goes through a detailed example of how to load a spreadsheet of directional data into R --- specifically, so that it can be processed with our geoDataFromFile function.



### PRELIMINARY WORK ###

# If you have just restarted R, then set the working directory to the 'geologyGeometry' directory, which contains subdirectories 'data', 'library', etc. Execute the following line of code to load our custom library and its dependencies.
source("library/all.R")

# In the 'data' subdirectory, find the file 'synthBroken.csv'. Make a copy of it, named 'synthFixed.csv'. Open synthFixed.csv in your favorite spreadsheet. We're going to tinker with it, leaving synthBroken.csv untouched.



### THREE IMPORTANT CONCEPTS ###

# First, I don't know whether R can load data in Microsoft Excel (.xls, .xlsx) formats. I never store my documents in any proprietary file formats, because I am worried that they will be unreadable 20 years later. I recommend that you store your data in a plain text format such as comma-separated values (.csv) or tab-separated values (.tsv). You can export from Excel to those formats.

# Second, R's main function for loading such data is read.table. It has some quirks, but they're not hard to deal with. Mainly you just need to configure your column headers correctly.

# Third, our library offers custom loading functions, based on read.table, designed to expedite the loading of certain kinds of geologic data. By pre-formatting your data a little bit in a spreadsheet, you can use these custom loaders and save yourself a lot of R trouble.



### COLUMN HEADERS ###

# For starters, let's try reading the file using R's read.table function. We get a data frame, but it is sad. The first couple of rows in the file are headers, but read.table has interpreted them as data.
read.table("data/synthFixed.csv", sep=",")

# read.table expects the first row to be headers, and every row after that to be data. So, in your spreadsheet, merge the first row into the second. Change 'number' to 'site number'. Change 'dec lat' to 'site dec lat'. Change the first 'strike' to 'dike strike'. And so on. Then delete the first row. Save synthFixed.csv. Then try reading it in again, informing read.table that there are proper headers now.
importData <- read.table("data/synthFixed.csv", sep=",", header=TRUE)
importData

# That's really all there is to it. If you're good at R, then you can stop reading this tutorial right now and forge your own path. For example, here are the paleomagnetic vectors, computed from their declinations and inclinations, with their Fisher alpha95 confidence regions, colored by longitude.
importPmags <- lapply(1:nrow(importData),
                      function(i) geoCartesianFromTrendPlungeDeg(c(data$ChRM.dec[[i]], data$ChRM.inc[[i]])))
rayEqualAreaRadiusPlot(points=importPmags,
                       radii=(importData$ChRM.alpha.95 * degree),
                       colors=hues(importData$site.dec.long))

# But things might get easier if you keep reading this tutorial.



### STRIKE AND DIP ###

# In the rest of this tutorial, we use geoDataFromFile instead of read.table. For now, geoDataFromFile produces nearly identical results as read.table, because we haven't tweaked the data file much yet.
geoDataFromFile("data/synthFixed.csv")

# In your spreadsheet, change the headers 'dike strike' and 'dike dip' to just 'strike' and 'dip'. Save the file and load it into R again. geoDataFromFile has renamed 'strike' and 'dip' to 'strikeDeg' and 'dipDeg', and it has computed the pole for us.
geoDataFromFile("data/synthFixed.csv")

# Unfortunately, some of those poles are wrong, because geoDataFromFile expects strikes and dips to follow the right-hand rule, which is being violated at four sites. The (66, 36, n) dike should be (246, 36, n). The (100, 41, ne) dike should be (280, 41, ne). The (91, 24, n) dike should be (271, 24, n). And the (37, 40, ne) dike doesn't make sense. I don't trust that datum.

# In your spreadsheet, change the strikes of the first three dikes just mentioned, so that they obey the right-hand rule. And delete the entire row for the dike that doesn't make any sense. Then reload.
geoDataFromFile("data/synthFixed.csv")

# If there were also rakes (pitches) in the file, then we would need to adjust them for the right-hand rule too. More information about that can be found in the bonus tutorial oriImporting.R.

# Notice that the pole is NA at the fourth site. Why? Because the dip is NA. In your spreadsheet, delete that entire row. Then reload.
importData <- geoDataFromFile("data/synthFixed.csv")
importData

# Possibly we should also delete the row where dike.alpha.95 is NA. It depends on what kind of analysis we're planning. Let's not delete it right now.

# Anyway, the dike data are ready to use. You can plot the dike poles, compute the sample mean, etc.
lineEqualAreaRadiusPlot(points=importData$pole,
                        radii=(importData$dike.alpha.95 * degree),
                        colors=hues(importData$site.dec.long))
geoStrikeDipDegFromCartesian(lineProjectedMean(importData$pole))



### TREND AND PLUNGE ###

# Now let's take a look at the paleomagnetic side of the data set. Notice that we've deleted two perfectly good paleomagnetic data, just because their corresponding dike data were damaged. Depending on what kind of analysis we're planning, maybe we should have kept them, maybe by splitting the paleomagnetic data into a separate file. Let's not worry about it right now.
importData

# In your spreadsheet, change ChRM dec and ChRM inc to trend and plunge. Save and reload here.
importData <- geoDataFromFile("data/synthFixed.csv")
importData

# Now you're ready to work. Here's a plot and a mean trend-plunge in degrees.
rayEqualAreaRadiusPlot(points=importData$direction,
                       radii=(importData$ChRM.alpha.95 * degree),
                       colors=hues(importData$site.dec.long))
geoTrendPlungeDegFromCartesian(rayProjectedMean(importData$direction))



### ROTATION? ###

# Finally, notice that the data frame has a rotation column. It's automatically computed, under the assumption that the pole and the direction are (nearly) perpendicular. In that case, this rotation is the orientation matrix of a plane-line pair or plane-ray pair.
importData

# But in this data set that's not supposed to happen. Each ray is a paleomagnetic direction, and there's no reason to assume that it is perpendicular to its corresponding dike pole. The rotation is still a mathematically valid object. You might even use it, if you wanted to assume that the direction and pole rotated together, as a rigid body. But if you don't want to study that, then feel free to delete the rotation column, like this.
importData$rotation <- NULL
importData



### WHAT HAVE WE LEARNED? ###

# By cleaning up your spreadsheet a little bit, you can help R and our library expedite your data loading process.



### SEE ALSO ###

# Bonus tutorial oriImportation.R does a similar exercise for orientational data.

# Bonus tutorial ellImportation.R does a similar exercise for ellipsoidal data.


