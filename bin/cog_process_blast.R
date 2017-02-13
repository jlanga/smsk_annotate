#!/usr/bin/env Rscript

library("tidyverse")
library("optparse")
library("stringr")

option_list = list(
    make_option(
        opt_str = c("-c", "--full_cog_table"),
        type = "character",
        default = NULL, 
        help = "TSV with the all the cog description cogs (cog_big_table.tsv)",
        metavar = "character"
    ),
    make_option(
        opt_str = c("-b", "--blast_table"),
        type = "character",
        default = NULL,
        help = "A TSV with the blastp agaist the COG fasta (blastp.tsv)",
        metavar = "character"
    ),
    make_option(
        opt_str = c("-o", "--output"),
        type = "character",
        default = NULL,
        help = "A TSV with the cog annotation of the sequences",
        metavar = "character"
    )
)   

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# opt$full_cog_table = "results/cog/cog_full_table.tsv"
# opt$blast_table = "results/cog/blastp.tsv"
# opt$output = "results/cog/cog_results.tsv"

# Read blast tables
read_blast <- function(filename){
    filename %>%
        read_tsv(
            file = .,
            col_names = c(
                "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                "qstart", "qend", "sstart", "send", "evalue", "bitscore")
        )
}

    
    
# Remove duplicated hits (pairs qseqid-sseqid) using dplyr's criterion
remove_duplicated_hits <- function(blast_table){
    blast_table %>% 
        group_by(qseqid, sseqid) %>%
        filter(row_number() == 1)
}

#From the blast results, take the identifier in which we are interested
get_ids <- function(raw_ids){
    # Convert the gi|whatever1|ref|whatever2| to whatever1
    ids <- rep(0, length(raw_ids)) # Initialize
    for(i in 1:length(raw_ids)){
        # Get an ID
        full_id <- raw_ids[i]
        # Split by |, take the first element in list, and then th
        id <- (str_split(full_id, "\\|") %>% unlist())[2]
        ids[i] <- id
    }
    return(as.integer(ids))
}


# Compute the cogs found
analyze_cog <- function(blast_data_frame, cog_full){
    blast_data_frame %>%
        mutate(protein_id = get_ids(sseqid)) %>%
        inner_join(x = cog_full, y = ., by = "protein_id")
}

cog_full <- read_tsv(opt$full_cog_table)

opt$blast_table %>%
    read_blast() %>%
    remove_duplicated_hits() %>%
    analyze_cog(cog_full) %>%
    write_tsv(
        path = opt$output, col_names= TRUE
    )    