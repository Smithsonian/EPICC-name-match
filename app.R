################################
#Load packages
################################
library(shiny)
library(DT)
library(dplyr)
library(stringr)
library(stringdist)
library(futile.logger)



################################
#Settings
################################

app_name <- "EPICC Name Match"
app_ver <- "0.2.1"
github_link <- "https://github.com/Smithsonian/EPICC-name-match"

options(stringsAsFactors = FALSE)

#how many cpu cores to use for matching
this_cpu_cores <- 2
threshold <- 4
method <- "osa"

logfile <- paste0("logs/", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")


################################
#Functions
################################

find_matching_str <- function(str_to_check, database, method = "osa", no_cores = 2){
  
  #remove these chars
  this_str <- gsub("[?!*]", "", as.character(str_to_check))
  
  if (this_str == "" || this_str == "-" || this_str == "NA" || is.na(this_str)){
    cat("Empty string, returning NAs.")
    results <- as.data.frame(cbind(NA, NA))
    names(results) <- c("match", "score")
    return(results)
  }else{
    
    if (method == "jw"){
      #Jaro-Winkler distance
      str_matches <- stringdist::stringdist(this_str, database[,1], nthread = no_cores, method = method, p = 0.1)
    }else{
      str_matches <- stringdist::stringdist(this_str, database[,1], nthread = no_cores, method = method)
    }
    
    #Add string to scores
    results <- cbind(database, str_matches)
    
    names(results) <- c("match", "score")
    return(results)
  }
}




################################
#UI
################################
ui <- fluidPage(
  
  # App title ----
  titlePanel(app_name),
  p("System to match scientific names to the EPICC Taxonomy."),
  tabsetPanel(type = "tabs",
              tabPanel("Match",  
                       fluidRow(
                         column(width = 8,
                                br(),
                                fluidRow(
                                  column(width = 4, 
                                         fileInput("taxonomy", "Select the EPICC taxonomy file to use",
                                                   multiple = FALSE,
                                                   accept = c("text/csv",
                                                              "text/comma-separated-values,text/plain",
                                                              ".csv")),
                                         uiOutput("taxonomy_info")
                                  ),
                                  column(width = 4, 
                                         fileInput("csv", "Select a csv file to match to the taxonomy",
                                                   multiple = FALSE,
                                                   accept = c("text/csv",
                                                              "text/comma-separated-values,text/plain",
                                                              ".csv")),
                                         uiOutput("csv_info")
                                  ),
                                  column(width = 4, 
                                         uiOutput("downloadData"),
                                         br()
                                  )
                                ),
                                uiOutput("tableheading"),
                                DT::dataTableOutput("table")
                         ),
                         column(width = 4, 
                                uiOutput("displayrow")
                         )
                       )
                       
                       ),
              tabPanel("Help", 
                       br(),
                       fluidRow(
                         column(width = 6, 
                                uiOutput("help1")
                         ),
                         column(width = 6, 
                                uiOutput("help2")
                         )
                      )
              )
          ),
  
  hr(),
  HTML(paste0("<p><a href=\"http://dpo.si.edu\" target = _blank><img src=\"circlelogo.png\"> Digitization Program Office</a> | ", app_name, " ver. ", app_ver, " | <a href=\"", github_link, "\" target = _blank>Source code</a></p>"))
)




################################
#Server
################################
server <- function(input, output) {
  
  #Set logging
  dir.create('logs', showWarnings = FALSE)
  flog.logger("enm", INFO, appender=appender.file(logfile))
  
  output$table <- DT::renderDataTable({
    req(input$taxonomy)
    req(input$csv)
    
    #Progress bar
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    
    progress$set(message = "Processing taxonomy", value = 0.01)
    
    if (!exists("synonym_list") && !exists("taxonomy")){
      taxonomy <- read.csv(input$taxonomy$datapath, header = TRUE, stringsAsFactors = FALSE, encoding = "UTF-8")
      
      taxonomy <- dplyr::select(taxonomy, c(
        "PHYLUM",
        "CLASS",
        "ORDER",
        "FAMILY",
        "GENUS",
        "SUBGENUS",
        "SPECIES",
        "AUTHOR",
        "SYNONYMS"))
      taxonomy$ID <- seq.int(nrow(taxonomy))
      
      #break synonyms
      taxon_list1 <- data.frame(taxonomy$ID, taxonomy$SYNONYMS)
      synonym_list <- data.frame(ID=as.numeric(), SYNONYMS = character(), stringsAsFactors=FALSE)
      
      no_names <- dim(taxon_list1)[1]
      
      progress_fraction <- (0.4/no_names)
      
      for (t in 1:no_names){
        t1 <- strsplit(taxon_list1[t,2], ",")[[1]]
        
        for (tt in 1:length(t1)){
          tsyn <- trimws(t1[tt])
          #cat(tsyn)
          if (length(tsyn)!=0){
            
            synonym_list <- rbind(synonym_list, c(ID = taxon_list1[t,1], SYNONYMS = tsyn))
            
            tsyn_split <- strsplit(tsyn, " ")[[1]]
            
            if (length(tsyn_split) == 3){
              #w subgenus
              synonym_list <- rbind(synonym_list, c(ID = taxon_list1[t,1], SYNONYMS = paste(tsyn_split[1], tsyn_split[3])))
              #w subspecies
              synonym_list <- rbind(synonym_list, c(ID = taxon_list1[t,1], SYNONYMS = paste(tsyn_split[1], tsyn_split[2])))
            }else if (length(tsyn_split) == 4){
              #w subgenus and subspecies
              synonym_list <- rbind(synonym_list, c(ID = taxon_list1[t,1], SYNONYMS = paste(tsyn_split[1], tsyn_split[3])))
            }
            
          }
          
          
        }
        progress$set(message = "Processing taxonomy", detail = paste0(" - ", t, " of ", no_names, " names"), value = round((t * progress_fraction), 3))
      }
      rm(taxon_list1)
      names(synonym_list) <- c("ID", "SYNONYMS")
      synonym_list <<- dplyr::filter(synonym_list, !is.na(SYNONYMS))
      
      taxonomy <<- taxonomy
    }
    
    

    progress$set(message = "Processing CSV file", value = 0.2)
    
    csvfile <- read.csv(input$csv$datapath, header = TRUE, stringsAsFactors = FALSE, encoding = "UTF-8")
    
    csvfile$sciname <- ""
    csvfile$accepted_name <- ""
    csvfile$synonym <- ""
    csvfile$fuzzy_match <- ""
    csvfile$phylum <- ""
    csvfile$class <- ""
    csvfile$order <- ""
    csvfile$family <- ""
    csvfile$genus <- ""
    csvfile$subgenus <- ""
    csvfile$species <- ""
    csvfile$author <- ""
    
    no_rows <- dim(csvfile)[1]
    
    resultsfile <- data.frame(matrix(ncol = dim(csvfile)[2], nrow = 0, data = ""))
    names(resultsfile) <- names(csvfile)

    
    progress$set(message = "Processing CSV file", value = 0.1)
    
    progress_fraction <- (0.6/no_rows)
    
    for (i in 1:no_rows){
      this_row <- csvfile[i,]
      this_taxon <- this_row$Taxonomy
      
      #remove ?
      this_taxon <- gsub(x = this_taxon, pattern = '[?]', replacement = "")
      
      #check if taxon has more than 1 name, pipe-separated
      this_taxon <- strsplit(x = this_taxon, split = "[|]")
      
      no_taxons <- length(this_taxon[[1]])
      
      if (no_taxons > 0){
      
        for (j in 1:no_taxons){
          this_taxon_name <- this_taxon[[1]][j]
          cat(paste(this_taxon_name, "\n"))
          flog.info(this_taxon_name, name = "enm")
          #Remove aff, cf, sp, etc...
          #Remove aff, cf, sp, etc...
          this_taxon_name <- gsub("\"", "", this_taxon_name, fixed = TRUE)
          this_taxon_name <- gsub(" sp.", "", this_taxon_name, fixed = TRUE)
          this_taxon_name <- gsub(" spp.", "", this_taxon_name, fixed = TRUE)
          this_taxon_name <- gsub(" cf. ", " ", this_taxon_name, fixed = TRUE)
          this_taxon_name <- gsub(" aff. ", " ", this_taxon_name, fixed = TRUE)
          this_taxon_name <- gsub(" af. ", " ", this_taxon_name, fixed = TRUE)
          this_taxon_name <- gsub(" near ", " ", this_taxon_name, fixed = TRUE)
          this_taxon_name <- gsub(" n.sp.", "", this_taxon_name, fixed = TRUE)
          this_taxon_name <- gsub(" n.", "", this_taxon_name, fixed = TRUE)
          this_taxon_name <- gsub(" n. sp.", "", this_taxon_name, fixed = TRUE)
          
          
          #Remove &
          if (grepl("[&]", this_taxon_name)){
            this_test <- list(strsplit(this_taxon_name, " ")[[1]])[[1]]
            #remove names on both sides of &
            this_taxon_name <- paste(this_test[-c(which(this_test == "&") - 1, which(this_test == "&"), which(this_test == "&") + 1)], collapse=' ')
          }
          
          #Remove var
          if (grepl("[var.]", this_taxon_name)){
            this_taxon_name <- strsplit(this_taxon_name, "var.")[[1]][1]
          }
          
          #remove whitespace
          this_taxon_name <- trimws(this_taxon_name)
          
          #Save name to try to find in column sciname
          this_row$sciname <- this_taxon_name
          
          if (this_taxon_name == "" || is.na(this_taxon_name)){
            resultsfile <- rbind(resultsfile, this_row)
            next
          }
          
          this_taxon_nowords <- length(strsplit(x = this_taxon_name, split = "[ ]")[[1]])
          
          if (this_taxon_nowords == 1){
            
            #For Genus
            #1 match to: Genus
            taxon_list <- taxonomy$GENUS
            this_match <- which(tolower(taxon_list) == tolower(this_taxon_name))
            if (length(this_match) > 0){
              this_row$accepted_name <- taxon_list[this_match[1]]
              
              this_row$synonym <- ""
              this_row$fuzzy_match <- ""
              
              this_row$phylum <- taxonomy$PHYLUM[this_match[1]]
              this_row$class <- taxonomy$CLASS[this_match[1]]
              this_row$order <- taxonomy$ORDER[this_match[1]]
              this_row$family <- taxonomy$FAMILY[this_match[1]]
              this_row$genus <- taxonomy$GENUS[this_match[1]]
              this_row$subgenus <- taxonomy$SUBGENUS[this_match[1]]
              this_row$species <- taxonomy$SPECIES[this_match[1]]
              this_row$author <- taxonomy$AUTHOR[this_match[1]]

              resultsfile <- rbind(resultsfile, this_row)
              flog.info("Genus match", name = "enm")
              next
            }
            
            
            #Fuzzy match
            #Genus
            taxon_list <- taxonomy$GENUS
            
            fuzzy_match <- find_matching_str(this_taxon_name, database = data.frame(taxon_list), method = method, no_cores = this_cpu_cores)
            results <- data.frame(match = fuzzy_match$match, score = as.numeric(fuzzy_match$score), stringsAsFactors = FALSE)
            #filter
            results_filtered <- dplyr::filter(results, score < threshold)
            #return top match
            results_filtered <- dplyr::top_n(x = results_filtered, wt = -score, n = 1)
            
            if (dim(results_filtered)[1] > 0){
              
              this_match <- which(results_filtered[1,]$match == taxon_list)[1]
              
              matched_row <- taxonomy[taxonomy$ID == this_match,]
              
              this_row$accepted_name <- matched_row$GENUS
              this_row$synonym <- ""
              
              this_row$fuzzy_match <- paste0(results_filtered[1,]$match, " (", results_filtered[1,]$score, ")")
              
              this_row$phylum <- matched_row$PHYLUM
              this_row$class <- matched_row$CLASS
              this_row$order <- matched_row$ORDER
              this_row$family <- matched_row$FAMILY
              this_row$genus <- matched_row$GENUS
              this_row$subgenus <- matched_row$SUBGENUS
              this_row$species <- matched_row$SPECIES
              this_row$author <- matched_row$AUTHOR
              
              resultsfile <- rbind(resultsfile, this_row)
              flog.info("Genus fuzzymatch", name = "enm")
              next
            }
            
            
            
            
          }else if (this_taxon_nowords == 2){
            
            
            #check if Genus subgenus
            #Cardium (Cerastoderma?) sp.
            taxon_list <- paste0(taxonomy$GENUS, " (", taxonomy$SUBGENUS, ")")
            this_match <- which(taxon_list == this_taxon_name)
            if (length(this_match) > 0){
              this_row$accepted_name <- taxon_list[this_match[1]]
              
              this_row$synonym <- ""
              this_row$fuzzy_match <- ""
              
              this_row$phylum <- taxonomy$PHYLUM[this_match[1]]
              this_row$class <- taxonomy$CLASS[this_match[1]]
              this_row$order <- taxonomy$ORDER[this_match[1]]
              this_row$family <- taxonomy$FAMILY[this_match[1]]
              this_row$genus <- taxonomy$GENUS[this_match[1]]
              this_row$subgenus <- taxonomy$SUBGENUS[this_match[1]]
              this_row$species <- taxonomy$SPECIES[this_match[1]]
              this_row$author <- taxonomy$AUTHOR[this_match[1]]
              
              resultsfile <- rbind(resultsfile, this_row)
              flog.info("Genus subgenus", name = "enm")
              next
            }
            
            
            #check if Genus subgenus, remove subgenus
            if (grepl("[(]", this_taxon_name)){
              taxon_partial <- paste(strsplit(this_taxon_name, " ")[[1]][1])
              taxon_list <- paste0(taxonomy$GENUS)
              this_match <- which(taxon_list == taxon_partial)
              if (length(this_match) > 0){
                this_row$accepted_name <- taxon_list[this_match[1]]
                
                this_row$synonym <- ""
                this_row$fuzzy_match <- ""
                
                this_row$phylum <- taxonomy$PHYLUM[this_match[1]]
                this_row$class <- taxonomy$CLASS[this_match[1]]
                this_row$order <- taxonomy$ORDER[this_match[1]]
                this_row$family <- taxonomy$FAMILY[this_match[1]]
                this_row$genus <- taxonomy$GENUS[this_match[1]]
                this_row$subgenus <- taxonomy$SUBGENUS[this_match[1]]
                this_row$species <- taxonomy$SPECIES[this_match[1]]
                this_row$author <- taxonomy$AUTHOR[this_match[1]]
                
                resultsfile <- rbind(resultsfile, this_row)
                flog.info("Genus subgenus, remove subgenus", name = "enm")
                next
              }
            }
            
            
            
            #2 match to: Genus epithet
            taxon_list <- paste0(taxonomy$GENUS, " ", taxonomy$SPECIES)
            this_match <- which(taxon_list == this_taxon_name)
            if (length(this_match) > 0){
              this_row$accepted_name <- taxon_list[this_match[1]]
              
              this_row$phylum <- taxonomy$PHYLUM[this_match[1]]
              this_row$class <- taxonomy$CLASS[this_match[1]]
              this_row$order <- taxonomy$ORDER[this_match[1]]
              this_row$family <- taxonomy$FAMILY[this_match[1]]
              this_row$genus <- taxonomy$GENUS[this_match[1]]
              this_row$subgenus <- taxonomy$SUBGENUS[this_match[1]]
              this_row$species <- taxonomy$SPECIES[this_match[1]]
              this_row$author <- taxonomy$AUTHOR[this_match[1]]

              resultsfile <- rbind(resultsfile, this_row)
              flog.info("Genus epithet", name = "enm")
              next
            }
            
            
            #match to: synonym
            this_match <- which(this_taxon_name == synonym_list[,2])
            if (length(this_match) > 0){
              
              this_id <- synonym_list[this_match[1],1]
              
              matched_row <- taxonomy[taxonomy$ID == this_id,]
              
              this_row$accepted_name <- paste0(matched_row$GENUS, " ", matched_row$SPECIES)
              this_row$synonym <- this_taxon_name

              this_row$fuzzy_match <- ""

              this_row$phylum <- matched_row$PHYLUM
              this_row$class <- matched_row$CLASS
              this_row$order <- matched_row$ORDER
              this_row$family <- matched_row$FAMILY
              this_row$genus <- matched_row$GENUS
              this_row$subgenus <- matched_row$SUBGENUS
              this_row$species <- matched_row$SPECIES
              this_row$author <- matched_row$AUTHOR

              resultsfile <- rbind(resultsfile, this_row)
              flog.info("Genus epithet by synonym", name = "enm")
              next
            }
            
            
            
            
            #Fuzzy match
            #Genus epithet
            taxon_list <- paste0(taxonomy$GENUS, " ", taxonomy$SPECIES)
            
            fuzzy_match <- find_matching_str(this_taxon_name, database = data.frame(taxon_list), method = method, no_cores = this_cpu_cores)
            results <- data.frame(match = fuzzy_match$match, score = as.numeric(fuzzy_match$score), stringsAsFactors = FALSE)
            #filter
            results_filtered <- dplyr::filter(results, score < threshold)
            #return top match
            results_filtered <- dplyr::top_n(x = results_filtered, wt = -score, n = 1)
            
            if (dim(results_filtered)[1] > 0){
              
              this_match <- which(results_filtered[1,]$match == taxon_list)[1]
              
              matched_row <- taxonomy[taxonomy$ID == this_match,]
              
              this_row$accepted_name <- paste0(matched_row$GENUS, " ", matched_row$SPECIES)
              this_row$synonym <- ""
              
              this_row$fuzzy_match <- paste0(results_filtered[1,]$match, " (", results_filtered[1,]$score, ")")
              
              this_row$phylum <- matched_row$PHYLUM
              this_row$class <- matched_row$CLASS
              this_row$order <- matched_row$ORDER
              this_row$family <- matched_row$FAMILY
              this_row$genus <- matched_row$GENUS
              this_row$subgenus <- matched_row$SUBGENUS
              this_row$species <- matched_row$SPECIES
              this_row$author <- matched_row$AUTHOR
              
              resultsfile <- rbind(resultsfile, this_row)
              flog.info("Genus epithet fuzzymatch", name = "enm")
              next
            }
            
            
            
            
            
          }else if (this_taxon_nowords == 3){
          
            #3 match to: Genus subgenus epithet
            taxon_list <- paste0(taxonomy$GENUS, " ", taxonomy$SUBGENUS, " ", taxonomy$SPECIES)
            #taxon_list <- gsub(" [(][)] ", " ", taxon_list)
            
            this_match <- which(taxon_list == this_taxon_name)
            if (length(this_match) > 0){
              this_row$accepted_name <- taxon_list[this_match[1]]
              
              this_row$synonym <- ""
              this_row$fuzzy_match <- ""
              
              this_row$phylum <- taxonomy$PHYLUM[this_match[1]]
              this_row$class <- taxonomy$CLASS[this_match[1]]
              this_row$order <- taxonomy$ORDER[this_match[1]]
              this_row$family <- taxonomy$FAMILY[this_match[1]]
              this_row$genus <- taxonomy$GENUS[this_match[1]]
              this_row$subgenus <- taxonomy$SUBGENUS[this_match[1]]
              this_row$species <- taxonomy$SPECIES[this_match[1]]
              this_row$author <- taxonomy$AUTHOR[this_match[1]]

              resultsfile <- rbind(resultsfile, this_row)
              flog.info("Genus subgenus epithet", name = "enm")
              next
            }
            
            #5 match to: Genus epithet Author
            taxon_list <- paste0(taxonomy$GENUS, " ", taxonomy$SPECIES, " ", taxonomy$AUTHOR)
            taxon_list <- gsub("[  ]", " ", taxon_list)
            
            this_match <- which(taxon_list == this_taxon_name)
            if (length(this_match) > 0){
              this_row$accepted_name <- taxon_list[this_match[1]]
              
              this_row$synonym <- ""
              this_row$fuzzy_match <- ""
              
              this_row$phylum <- taxonomy$PHYLUM[this_match[1]]
              this_row$class <- taxonomy$CLASS[this_match[1]]
              this_row$order <- taxonomy$ORDER[this_match[1]]
              this_row$family <- taxonomy$FAMILY[this_match[1]]
              this_row$genus <- taxonomy$GENUS[this_match[1]]
              this_row$subgenus <- taxonomy$SUBGENUS[this_match[1]]
              this_row$species <- taxonomy$SPECIES[this_match[1]]
              this_row$author <- taxonomy$AUTHOR[this_match[1]]

              resultsfile <- rbind(resultsfile, this_row)
              flog.info("Genus epithet Author", name = "enm")
              next
            }
            
            
            #match to: Genus epithet Author but remove author
            taxon_list <- paste0(taxonomy$GENUS, " ", taxonomy$SPECIES)
            taxon_list <- gsub("[  ]", " ", taxon_list)
            
            taxon_partial <- paste(strsplit(this_taxon_name, " ")[[1]][1], strsplit(this_taxon_name, " ")[[1]][2])
            
            this_match <- which(taxon_list == taxon_partial)
            if (length(this_match) > 0){
              this_row$accepted_name <- taxon_list[this_match[1]]
              
              this_row$synonym <- ""
              this_row$fuzzy_match <- ""
              
              this_row$phylum <- taxonomy$PHYLUM[this_match[1]]
              this_row$class <- taxonomy$CLASS[this_match[1]]
              this_row$order <- taxonomy$ORDER[this_match[1]]
              this_row$family <- taxonomy$FAMILY[this_match[1]]
              this_row$genus <- taxonomy$GENUS[this_match[1]]
              this_row$subgenus <- taxonomy$SUBGENUS[this_match[1]]
              this_row$species <- taxonomy$SPECIES[this_match[1]]
              this_row$author <- taxonomy$AUTHOR[this_match[1]]

              resultsfile <- rbind(resultsfile, this_row)
              flog.info("Genus epithet Author but remove author", name = "enm")
              next
            }
            
            
            
            #match to: Genus subgenus epithet bur remove subgenus
            taxon_list <- paste0(taxonomy$GENUS, " ", taxonomy$SPECIES)
            taxon_list <- gsub("[  ]", " ", taxon_list)
            
            taxon_partial <- paste(strsplit(this_taxon_name, " ")[[1]][1], strsplit(this_taxon_name, " ")[[1]][3])
            
            this_match <- which(taxon_list == taxon_partial)
            if (length(this_match) > 0){
              this_row$accepted_name <- taxon_list[this_match[1]]
              
              this_row$synonym <- ""
              this_row$fuzzy_match <- ""
              
              this_row$phylum <- taxonomy$PHYLUM[this_match[1]]
              this_row$class <- taxonomy$CLASS[this_match[1]]
              this_row$order <- taxonomy$ORDER[this_match[1]]
              this_row$family <- taxonomy$FAMILY[this_match[1]]
              this_row$genus <- taxonomy$GENUS[this_match[1]]
              this_row$subgenus <- taxonomy$SUBGENUS[this_match[1]]
              this_row$species <- taxonomy$SPECIES[this_match[1]]
              this_row$author <- taxonomy$AUTHOR[this_match[1]]

              resultsfile <- rbind(resultsfile, this_row)
              flog.info("Genus subgenus epithet bur remove subgenus", name = "enm")
              next
            }
            
            
            #match to: Genus epithet subspecies but remove subspecies
            taxon_list <- paste0(taxonomy$GENUS, " ", taxonomy$SPECIES)
            taxon_list <- gsub("[  ]", " ", taxon_list)
            
            taxon_partial <- paste(strsplit(this_taxon_name, " ")[[1]][1], strsplit(this_taxon_name, " ")[[1]][2])
            
            this_match <- which(taxon_list == taxon_partial)
            if (length(this_match) > 0){
              this_row$accepted_name <- taxon_list[this_match[1]]
              
              this_row$synonym <- ""
              this_row$fuzzy_match <- ""
              
              this_row$phylum <- taxonomy$PHYLUM[this_match[1]]
              this_row$class <- taxonomy$CLASS[this_match[1]]
              this_row$order <- taxonomy$ORDER[this_match[1]]
              this_row$family <- taxonomy$FAMILY[this_match[1]]
              this_row$genus <- taxonomy$GENUS[this_match[1]]
              this_row$subgenus <- taxonomy$SUBGENUS[this_match[1]]
              this_row$species <- taxonomy$SPECIES[this_match[1]]
              this_row$author <- taxonomy$AUTHOR[this_match[1]]

              resultsfile <- rbind(resultsfile, this_row)
              flog.info("Genus epithet subspecies but remove subspecies", name = "enm")
              next
            }
            
            
            
            #match to: synonym without author/subspecies
            taxon_partial <- paste(strsplit(this_taxon_name, " ")[[1]][1], strsplit(this_taxon_name, " ")[[1]][2])
            this_match <- which(taxon_partial == synonym_list[,2])
            if (length(this_match) > 0){
              
              this_id <- synonym_list[this_match[1],1]
              
              matched_row <- taxonomy[taxonomy$ID == this_id,]
              
              this_row$accepted_name <- paste0(matched_row$GENUS, " ", matched_row$SPECIES)
              this_row$synonym <- taxon_partial

              this_row$fuzzy_match <- ""

              this_row$phylum <- matched_row$PHYLUM
              this_row$class <- matched_row$CLASS
              this_row$order <- matched_row$ORDER
              this_row$family <- matched_row$FAMILY
              this_row$genus <- matched_row$GENUS
              this_row$subgenus <- matched_row$SUBGENUS
              this_row$species <- matched_row$SPECIES
              this_row$author <- matched_row$AUTHOR

              resultsfile <- rbind(resultsfile, this_row)
              flog.info("synonym without author/subspecies", name = "enm")
              next
            }
            
            
            #match to: Genus subgenus epithet to synonym
            taxon_partial <- paste(strsplit(this_taxon_name, " ")[[1]][1], strsplit(this_taxon_name, " ")[[1]][3])
            this_match <- which(taxon_partial == synonym_list[,2])
            if (length(this_match) > 0){
              
              this_id <- synonym_list[this_match[1],1]
              
              matched_row <- taxonomy[taxonomy$ID == this_id,]
              
              this_row$accepted_name <- paste0(matched_row$GENUS, " ", matched_row$SPECIES)
              this_row$synonym <- taxon_partial
              
              this_row$fuzzy_match <- ""
              
              this_row$phylum <- matched_row$PHYLUM
              this_row$class <- matched_row$CLASS
              this_row$order <- matched_row$ORDER
              this_row$family <- matched_row$FAMILY
              this_row$genus <- matched_row$GENUS
              this_row$subgenus <- matched_row$SUBGENUS
              this_row$species <- matched_row$SPECIES
              this_row$author <- matched_row$AUTHOR

              resultsfile <- rbind(resultsfile, this_row)
              flog.info("Genus subgenus epithet to synonym", name = "enm")
              next
            }
            
            
            
            
            }else{
            
            
            #4 match to: Genus subgenus epithet Author
            taxon_list <- paste0(taxonomy$GENUS, " (", taxonomy$SUBGENUS, ") ", taxonomy$SPECIES, " ", taxonomy$AUTHOR)
            taxon_list <- gsub(" [(][)] ", " ", taxon_list)
            taxon_list <- gsub("[  ]", " ", taxon_list)
            
            this_match <- which(taxon_list == this_taxon_name)
            if (length(this_match) > 0){
              this_row$accepted_name <- taxon_list[this_match[1]]

              this_row$synonym <- ""
              this_row$fuzzy_match <- ""              

              this_row$phylum <- taxonomy$PHYLUM[this_match[1]]
              this_row$class <- taxonomy$CLASS[this_match[1]]
              this_row$order <- taxonomy$ORDER[this_match[1]]
              this_row$family <- taxonomy$FAMILY[this_match[1]]
              this_row$genus <- taxonomy$GENUS[this_match[1]]
              this_row$subgenus <- taxonomy$SUBGENUS[this_match[1]]
              this_row$species <- taxonomy$SPECIES[this_match[1]]
              this_row$author <- taxonomy$AUTHOR[this_match[1]]

              resultsfile <- rbind(resultsfile, this_row)
              flog.info("Genus subgenus epithet Author", name = "enm")
              next
            }
            
            
            
            #match to: Genus epithet subspecies Author but remove subspecies
            taxon_list <- paste0(taxonomy$GENUS, " ", taxonomy$SPECIES, " ", taxonomy$AUTHOR)
            taxon_list <- gsub("[  ]", " ", taxon_list)
            
            taxon_partial <- paste(strsplit(this_taxon_name, " ")[[1]][1], strsplit(this_taxon_name, " ")[[1]][2], strsplit(this_taxon_name, " ")[[1]][4])
            
            this_match <- which(taxon_list == taxon_partial)
            if (length(this_match) > 0){
              this_row$accepted_name <- taxon_list[this_match[1]]
              
              this_row$synonym <- ""
              this_row$fuzzy_match <- ""
              
              this_row$phylum <- taxonomy$PHYLUM[this_match[1]]
              this_row$class <- taxonomy$CLASS[this_match[1]]
              this_row$order <- taxonomy$ORDER[this_match[1]]
              this_row$family <- taxonomy$FAMILY[this_match[1]]
              this_row$genus <- taxonomy$GENUS[this_match[1]]
              this_row$subgenus <- taxonomy$SUBGENUS[this_match[1]]
              this_row$species <- taxonomy$SPECIES[this_match[1]]
              this_row$author <- taxonomy$AUTHOR[this_match[1]]

              resultsfile <- rbind(resultsfile, this_row)
              flog.info("Genus epithet subspecies Author but remove subspecies", name = "enm")
              next
            }
            
            
            
            
            
            #match to: Genus epithet subspecies Author but remove subspecies and author
            taxon_list <- paste0(taxonomy$GENUS, " ", taxonomy$SPECIES)
            taxon_list <- gsub("[  ]", " ", taxon_list)
            
            taxon_partial <- paste(strsplit(this_taxon_name, " ")[[1]][1], strsplit(this_taxon_name, " ")[[1]][2])
            
            this_match <- which(taxon_list == taxon_partial)
            if (length(this_match) > 0){
              this_row$accepted_name <- taxon_list[this_match[1]]
              this_row$synonym <- ""
              this_row$fuzzy_match <- ""
              this_row$phylum <- taxonomy$PHYLUM[this_match[1]]
              this_row$class <- taxonomy$CLASS[this_match[1]]
              this_row$order <- taxonomy$ORDER[this_match[1]]
              this_row$family <- taxonomy$FAMILY[this_match[1]]
              this_row$genus <- taxonomy$GENUS[this_match[1]]
              this_row$subgenus <- taxonomy$SUBGENUS[this_match[1]]
              this_row$species <- taxonomy$SPECIES[this_match[1]]
              this_row$author <- taxonomy$AUTHOR[this_match[1]]

              resultsfile <- rbind(resultsfile, this_row)
              flog.info("Genus epithet subspecies Author but remove subspecies and author", name = "enm")
              next
            }
            
            
            
            
            #match to: synonym
            taxon_partial <- paste(strsplit(this_taxon_name, " ")[[1]][1], strsplit(this_taxon_name, " ")[[1]][3])
            this_match <- which(taxon_partial == synonym_list[,2])
            if (length(this_match) > 0){
              
              this_id <- synonym_list[this_match[1],1]
              
              matched_row <- taxonomy[taxonomy$ID == this_id,]
              
              this_row$accepted_name <- paste0(matched_row$GENUS, " ", matched_row$SPECIES)
              this_row$synonym <- taxon_partial
              this_row$fuzzy_match <- ""
              this_row$phylum <- matched_row$PHYLUM
              this_row$class <- matched_row$CLASS
              this_row$order <- matched_row$ORDER
              this_row$family <- matched_row$FAMILY
              this_row$genus <- matched_row$GENUS
              this_row$subgenus <- matched_row$SUBGENUS
              this_row$species <- matched_row$SPECIES
              this_row$author <- matched_row$AUTHOR

              resultsfile <- rbind(resultsfile, this_row)
              flog.info("synonym", name = "enm")
              next
            }
            
            
            
            
            ######
            #Fuzzy match
            ######
            #match Genus [Subgenus] species [author/subspecies/etc] to synonym
            #is second word in paren?
            if (grepl("[(]", strsplit(this_taxon_name, " ")[[1]][2])){
              taxon_partial <- paste(strsplit(this_taxon_name, " ")[[1]][1], strsplit(this_taxon_name, " ")[[1]][3])
            }else{
              taxon_partial <- paste(strsplit(this_taxon_name, " ")[[1]][1], strsplit(this_taxon_name, " ")[[1]][2])
            }
            
            fuzzy_match <- find_matching_str(taxon_partial, database = data.frame(synonym_list[,2]), method = method, no_cores = this_cpu_cores)
            results <- data.frame(match = fuzzy_match$match, score = as.numeric(fuzzy_match$score), stringsAsFactors = FALSE)
            #filter
            results_filtered <- dplyr::filter(results, score < threshold)
            #return top match
            results_filtered <- dplyr::top_n(x = results_filtered, wt = -score, n = 1)
            
            if (dim(results_filtered)[1] == 1){
              
              this_match <- which(results_filtered$match == synonym_list[,2])
              
              this_id <- synonym_list[this_match[1],1]
              
              matched_row <- taxonomy[taxonomy$ID == this_id,]
              
              this_row$accepted_name <- paste0(matched_row$GENUS, " ", matched_row$SPECIES)
              this_row$synonym <- results_filtered$match
              
              this_row$fuzzy_match <- paste0(results_filtered$match, " (", results_filtered$score, ")")
              
              this_row$phylum <- matched_row$PHYLUM
              this_row$class <- matched_row$CLASS
              this_row$order <- matched_row$ORDER
              this_row$family <- matched_row$FAMILY
              this_row$genus <- matched_row$GENUS
              this_row$subgenus <- matched_row$SUBGENUS
              this_row$species <- matched_row$SPECIES
              this_row$author <- matched_row$AUTHOR
              
              resultsfile <- rbind(resultsfile, this_row)
              flog.info("Genus [Subgenus] species [author/subspecies/etc] to synonym fuzzy_match", name = "enm")
              next
            }
            
            
          }
          
          #nothing found, return the same
          this_row$accepted_name <- ""
          this_row$synonym <- ""
          this_row$fuzzy_match <- ""
          this_row$phylum <- ""
          this_row$class <- ""
          this_row$order <- ""
          this_row$family <- ""
          this_row$genus <- ""
          this_row$subgenus <- ""
          this_row$species <- ""
          this_row$author <- ""
          resultsfile <- rbind(resultsfile, this_row)
          flog.info("No match found", name = "enm")
          next
          
        }
          
      }else{
        
        this_row$accepted_name <- ""
        this_row$synonym <- ""
        this_row$fuzzy_match <- ""
        this_row$phylum <- ""
        this_row$class <- ""
        this_row$order <- ""
        this_row$family <- ""
        this_row$genus <- ""
        this_row$subgenus <- ""
        this_row$species <- ""
        this_row$author <- ""
        resultsfile <- rbind(resultsfile, this_row)
        flog.info("No match found", name = "enm")
        next
        
      }
      
      progress$set(message = "Analyzing ", detail = paste0(i, " of ", no_rows, " rows"), value = 0.4 + round((i * progress_fraction), 4))
      
    }
    
    progress$close()
    
    resultsfile <<- resultsfile
    
    resultssummary <- data.frame(irn = resultsfile$irn, Taxonomy = resultsfile$Taxonomy, sciname = resultsfile$sciname, accepted_name = resultsfile$accepted_name, synonym = resultsfile$synonym, fuzzy_match = resultsfile$fuzzy_match)

    DT::datatable(resultssummary, escape = FALSE, options = list(searching = TRUE, ordering = TRUE, pageLength = 15), rownames = FALSE, selection = 'single')
    
    })
  
  
  
  
  
    
  output$tableheading <- renderUI({
    req(input$taxonomy)
    req(input$csv)
    
    tagList(
      h3("Results"),
      p("Click on a row to see the details")
    )
  })
  
  
  
  
  
  
  
  #display row info
  output$displayrow <- renderUI({
    req(input$table_rows_selected)
    
    this_row <- resultsfile[input$table_rows_selected, ]  
    
    list_data <- ""
    
    for (i in 1:dim(this_row)[2]){
      list_data <- paste(list_data, paste0("<dt>", names(this_row[i]), "</dt><dd>", this_row[i], "</dd>"))
    }
    
    tagList(
        HTML("<div class=\"panel panel-success\"> <div class=\"panel-heading\"> <h3 class=\"panel-title\">Row selected</h3> </div> <div class=\"panel-body\"><dl class=\"dl-horizontal\">"),
        HTML(list_data),
        HTML("</dl></div></div>")
      )
    })
  
  
  
  
  
  
  # Downloadable csv of rows with no issues
  output$downloadData2 <- downloadHandler(
    filename = function() {
      paste0("results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      write.csv(resultsfile, file, row.names = FALSE)
    }
  )
  
  
  
  
  output$downloadData <- renderUI({
    req(input$taxonomy)
    req(input$csv)
    tagList(
      HTML("<label>Download csv file with results</label>"),
      br(),
      downloadButton("downloadData2", "Download CSV", class = "btn-primary")  
    )
    
  })
  

  
  
  
  
  
  output$csv_info <- renderUI({
    if (is.null(input$csv)){
      HTML("<div class=\"panel panel-primary\"> <div class=\"panel-heading\"> <h3 class=\"panel-title\">The CSV file should be encoded using \"UTF-8\" and have at least 1 column</h3> </div> <div class=\"panel-body\"> <ul>
          <li>Taxonomy</li>
           </ul></div></div>")
    }
    })
  
  
  
  
  output$taxonomy_info <- renderUI({
    if (is.null(input$taxonomy)){
      HTML("<div class=\"panel panel-primary\"> <div class=\"panel-heading\"> <h3 class=\"panel-title\">The EPICC Taxonomy file should be encoded using \"UTF-8\" and have at least these 9 columns</h3> </div> <div class=\"panel-body\"> <ul>
           <li>PHYLUM</li>
           <li>CLASS</li>
           <li>ORDER</li>
           <li>FAMILY</li>
           <li>GENUS</li>
           <li>SUBGENUS</li>
           <li>SPECIES</li>
           <li>AUTHOR</li>
           <li>SYNONYMS</li>
           </ul></div></div>")
    }
    })
  
  
  
  
  
  output$help1 <- renderUI({
    HTML("<div class=\"panel panel-primary\"> <div class=\"panel-heading\"> <h3 class=\"panel-title\">Process</h3></div><div class=\"panel-body\">
         <p>This system will take the string in the column \"Taxonomy\" and match it with the Taxonomy from EPICC. The process tries to find a match taking into account the variety of ways that a scientific name can appear. 
         <p>The match is done by first removing words and characters like:
         <ul>
         <li>?</li>
         <li>af.</li>
         <li>cf.</li>
         <li>spp.</li>
         <li>n. sp.</li>
         </ul>
         <p>Then, the system tries to match the string by looking at possible ways to match:
         <ul>
         <li>Genus</li>              
         <li>Genus (Subgenus)</li>
         <li>Genus species</li>
         <li>Genus species Author</li>
         <li>Genus (Subgenus) species</li>
         <li>Genus (Subgenus) species Author</li>
         <li>Synonym</li>
         </ul>
        <p>To match another CSV you don't need to reload the taxonomy file. Just browse for a new CSV, the taxonomy is already in memory.</p>
         </div></div>")
  })
  
  output$help2 <- renderUI({
    HTML("
         <div class=\"panel panel-primary\"> <div class=\"panel-heading\"> <h3 class=\"panel-title\">Results file</h3></div><div class=\"panel-body\">
        <p>The results file will return the same columns as the input file, plus these columns:
          <dl class=\"dl-horizontal\">
              <dt>sciname</dt><dd>The string used to find a match</dd>
              <dt>accepted_name</dt><dd>Concatenated (Genus Species) from the taxonomy file</d>
              <dt>synonym</dt><dd>The synonym matched to the sciname</dd>
              <dt>fuzzy_match</dt><dd>Name matched to using a fuzzy matching algorithm (score in parenthesis) *</d>
              <dt>phylum</dt><dd>The Phylum of the accepted_name from the taxonomy file</dd>
              <dt>class</dt><dd>The Class of the accepted_name from the taxonomy file</dd>
              <dt>order</dt><dd>The Order of the accepted_name from the taxonomy file</dd>
              <dt>family</dt><dd>The Family of the accepted_name from the taxonomy file</dd>
              <dt>genus</dt><dd>The Genus of the accepted_name from the taxonomy file</dd>
              <dt>subgenus</dt><dd>The Subgenus of the accepted_name from the taxonomy file</dd>
              <dt>species</dt><dd>The Species of the accepted_name from the taxonomy file</dd>
              <dt>author</dt><dd>The \"AUTHOR\" of the accepted_name from the taxonomy file</dd>
         </dl>
        <p>If the Taxonomy column had more than one name, separated by pipes (|), each name is returned in a separate row.</p>
        <p>* Fuzzy matching uses the \"osa\" method in the function stringdist::stringdist(), for details see van der Loo (2014).</p>
        <pre>van der Loo, Mark (2014). The stringdist Package for Approximate String Matching</a>. The R Journal 6: 111-122.</pre>
         </div></div>")
  })
  
}

# Create Shiny app ----
shinyApp(ui = ui, server = server)
