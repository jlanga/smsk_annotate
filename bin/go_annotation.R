#!/usr/bin/env Rscript
#title: go_annotation.R
#author: Jorge Langa
#description: Set of functions to do the Gene Ontology annotation (tables and 
#   plots) given the table resulting from retrieving a dataframe from Uniprot 
#   using the bioservices python package. The result should be 3 plots (one per
#   category), 3 tables (one per category) and one table with the summary of
#   the three categories.

library(readr)

arguments <- commandArgs(trailingOnly = TRUE)
if(length(arguments) != 1){
    message("ERROR! Pass me only the bioservices file!")
    quit(save = "no", status = 1, runLast = FALSE)
}

# Required functions

read_bioservices_go <- function(filename){
    
    #Load libraries if necessary
    require(dplyr)
    require(tidyr)
    require(readr)
    
    #Column names in which we are interested
    go_columns_old <- c(
        "qseqid", "Entry",
        "Gene ontology (biological process)",
        "Gene ontology (cellular component)",
        "Gene ontology (molecular function)"
    )
    
    go_table <- filename %>%
        read_tsv() %>%
        select(one_of(go_columns_old)) %>% # Take only those columns
        rename( # Rename the ugly ones
            bp = `Gene ontology (biological process)`,
            cc = `Gene ontology (cellular component)`,
            mf = `Gene ontology (molecular function)`
        ) %>%
        transform(bp = strsplit(bp, "; ")) %>% # Split multiple GOs into multiple lines
        unnest(bp) %>% 
        transform(cc = strsplit(cc, "; ")) %>%
        unnest(cc) %>%
        transform(mf = strsplit(mf, "; ")) %>%
        unnest(mf) %>%
        gather(category, temp, bp:mf) %>% # Make the table long
        filter(!is.na(temp)) %>% # Remove NAs
        mutate( # change bp, cc, mf into the full word via factors
            category = factor(
                category,
                labels = c("Biological Process", "Cellular Component", "Molecular Function")
            )
        ) %>%
        group_by(qseqid, Entry, category, temp) %>% # Remove dups
        filter(row_number() == 1) %>%
        ungroup() %>%
        extract( # Use a regex to separate the GO term and the id
            temp,
            into = c("description", "go_id"),
            regex = "^([:print:]+) \\[(GO:[0-9]+)\\]$"
        )
    
    return(go_table)

}

count_gos <- function(go_table){
    go_table %>%
        group_by(description) %>%
        summarise(counts = n()) %>%
        arrange(desc(counts)) %>%
        mutate(percent = counts / sum(counts) * 100)
}

count_gos_by_category <- function(go_table){
    go_table %>%
        group_by(category) %>%
        summarise(counts = n()) %>%
        arrange(desc(counts)) %>%
        mutate(percent = counts / sum(counts) * 100)
}

plot_top <- function(go_table, category_to_filter, top = 25, outfile){
    require(dplyr)
    require(ggplot2)
    
    go_top <- go_table %>%
        filter(category == category_to_filter) %>% # Take only the category we want
        select(description) %>% # Retain only the description
        group_by(description) %>% # Summarise by counts
        summarise(counts = n()) %>%
        arrange(desc(counts)) %>%
        slice(1:top) # take the top
    
    go_top %>% # Make a barplot of description-counts
        ggplot(
            aes(x = description, y = counts)
        ) +
        geom_bar(stat = "identity") +
        ggtitle(category_to_filter) + # Add title
        theme( # Rotate 90º the x labels
            axis.text.x = element_text(
                angle = 90,
                hjust = 1,
                vjust = 1
            )
        )
    
    #Save to file
    ggsave(file = outfile, width = 297, height = 210, units = "mm")
}


# Main part of the script
# Read the raw table and tidy it
filename <- arguments[1]

go_table <- read_bioservices_go(filename)    

# Store it
go_table %>%
    write_tsv("go_annotation.tsv")


# Plot BP
go_table %>%
    plot_top(
        category = "Biological Process",
        top = 25,
        outfile = "bp.pdf"
    )

# Plot CC
go_table %>%
    plot_top(
        category = "Cellular Component",
        top = 25,
        outfile = "cc.pdf"
    )

# Plot MF
go_table %>%
    plot_top(
        category = "Molecular Function",
        top = 25,
        outfile = "mf.pdf"
    )

# Counts for BP
go_table %>%
    filter(category == "Biological Process") %>%
    count_gos() %>%
    write_tsv("bp_counts.tsv")

# Counts for CC
go_table %>%
    filter(category == "Cellular Component") %>%
    count_gos() %>%
    write_tsv("cc_counts.tsv")

# Counts for MF
go_table %>%
    filter(category == "Molecular Function") %>%
    count_gos() %>%
    write_tsv("mf_counts.tsv")

# Counts for categories
go_table %>%
    select(category) %>%
    group_by(category) %>%
    summarise(counts = n()) %>%
    mutate(percent = counts / sum(counts) * 100) %>%
    write_tsv("category_counts.tsv")
