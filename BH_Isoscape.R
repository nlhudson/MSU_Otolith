
#### Isoscape for Big Hole River ####

# install.packages("devtools") # for downlaoding packages from GitHub
# install.packages("remotes") # second option for downlaoding packages from GitHub
library(devtools)
library(remotes)

install_version("SSN", "1.1.15") # installing nessesary dependencies for openStars NEED SSN, rgdal, rgrass7

remotes::install_version("igraph", version = "1.2.6")


install_version()

## Installing openSTARS
devtools::install_github("MiKatt/openSTARS", ref = "dev")

#### THIS DOES NOT WORK ON APPLE SILICON SINCE THE ARCHITECTURE FOR OPENSTARS IN INCOMPATIBLE, AND MUST USE INTEL ARCHITECTURE
####. IF I WANT TO INSTALL OPENSTARS DEPENDENCIES I NEED TO USE ROSETTA 2 TO SWITCH TO X86_64 ARCHITECTURE AND REINSTALL HOMEBREW


