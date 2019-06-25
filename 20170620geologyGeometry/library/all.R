


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# This file is nothing more than a convenient way to load our entire R library and its dependencies (excluding the C part).



# Load external libraries. For example, expm must be loaded before deformations.R, because the latter defines defExp to be expm.
library("rgl")
library("fields")
library("MASS")
library("ICSNP")
#library("FRB") <-- archived and defunct
library("expm")
library("Directional")
library("pracma")

# Load our geologyGeometry library, excluding the C part.
source("library/miscellany.R")
source("library/rays.R")
source("library/lines.R")
source("library/rotations.R")
source("library/orientations.R")
source("library/ellipsoids.R")
source("library/deformations.R")
source("library/geology.R")


