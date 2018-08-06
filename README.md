# EPICC Name Match

A Shiny app to help the NMNH Paleo collection on their digitization efforts. The app tries to match the scientific names from labels in the collection to the taxonomy by EPICC ([Eastern Pacific Invertebrate Communities of the Cenozoic](https://epicc.berkeley.edu/)). 

This app takes the string in the column \"Taxonomy\" of a csv file and matches it with the Taxonomy from EPICC. The process tries to find a match taking into account the variety of ways that a scientific name can appear. 

The process tries to match the string by looking at possible ways a scientific name is written in a collection label by trying to match:
         
  * Genus
  * Genus (Subgenus)
  * Genus species
  * Genus species Author
  * Genus (Subgenus) species
  * Genus (Subgenus) species Author
  * Synonym

If no matches are found, the system will try to find a match by performing an approximate string matching to a scientific name or a synonym.

This app requires the following R packages:

 * shiny
 * DT
 * dplyr
 * stringr
 * stringdist
 * futile.logger

To install:

```R
install.packages("shiny", "DT", "dplyr", "stringr", "stringdist", "futile.logger")
```