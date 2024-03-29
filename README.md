**This repo is now archived since the same results can be achieved using OpenRefine in a simpler way.**

# EPICC Name Match

A Shiny app to help the NMNH Paleo collection on their digitization efforts. The app tries to match the scientific names from labels in the collection to the taxonomy by EPICC ([Eastern Pacific Invertebrate Communities of the Cenozoic](https://epicc.berkeley.edu/)). 

This app takes the string in the column \"Taxonomy\" of a csv file and matches it with the Taxonomy from EPICC. The process tries to find a match taking into account the variety of ways that a scientific name can appear. 

![epicc](https://user-images.githubusercontent.com/2302171/43724155-4caf6298-9967-11e8-8cf0-ce06fd478b01.png)

The process tries to match the string by looking at possible ways a scientific name is written in a collection label by trying to match:
         
  * Genus
  * Genus (Subgenus)
  * Genus species
  * Genus species Author
  * Genus (Subgenus) species
  * Genus (Subgenus) species Author
  * Synonym

If no matches are found, the system will try to find a match by performing an approximate string matching to a scientific name or a synonym.

## Running in local computer

To test the app locally, without the need of a server, just install R and Shiny. Then, run a command that will download the source files from Github. 

R version 3.3 or better is required. After starting R, copy and paste these commands:

```R
install.packages(c("shiny", "DT", "dplyr", "stringr", "stringdist", "futile.logger"))

library(shiny)
runGitHub("EPICC-name-match", "Smithsonian")
```

Please note that the installation of the required packages may take a few minutes to download and install.

## Installation in Shiny server

This app requires the following R packages:

 * shiny
 * DT
 * dplyr
 * stringr
 * stringdist
 * futile.logger

To install:

```R
install.packages(c("shiny", "DT", "dplyr", "stringr", "stringdist", "futile.logger"))
```
