#!/usr/bin/env Rscript
library("tidyverse")
library("stringr")
library("optparse")

option_list = list(
    make_option(
        opt_str = c("-c", "--cog_table"),
        type = "character",
        default = NULL, 
        help = "csv with the list of cogs (cog2003-2014.csv)",
        metavar = "character"
    ),
    make_option(
        opt_str = c("-n", "--cog_names"),
        type = "character",
        default = NULL,
        help = "A TSV cog names (cognames2003-2014.tab)",
        metavar = "character"
    ),
    make_option(
        opt_str = c("-f", "--functions"),
        type= "character",
        default = NULL, 
        help= "A TSV with the cog functions (fun2003-2014.tab)",
        metavar="character"
    ),
    make_option(
        opt_str = c("-g", "--genomes"),
        type = "character",
        default = NULL,
        help= "TSV file with the genomes in COG (genomes2003-2014.tab)",
        metavar= "character"
    ),
    make_option(
        opt_str = c("-p", "--proteins"),
        type = "character",
        default = NULL,
        help= "TSV file with the RefSeq - Protein id mapping (proteins2003-2014.tab)",
        metavar= "character"
    ),
    make_option(
        opt_str = c("-o", "--output"),
        type = "character",
        default = NULL,
        help= "TSV resulting of merging the tables",
        metavar= "character"
    )
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


cog <- read_csv(
        file = opt$cog_table,
        col_names = c(
            "domain_id", "genome_name", "protein_id",
            "protein_length", "domain_start", "domain_end",
            "cog_id", "membership_class", "unknown"
        )
    ) %>%
    select(-unknown)



# COG names
# Table schema:
# <COG-id> <functional-class> <COG-annotation>
cog_names <- read_tsv(file = opt$cog_names
    , col_names = FALSE, comment = "#") %>%
    select(
        cog_id = X1,
        functional_class = X2,
        cog_annotation = X3
    ) %>%
    mutate( # Convert functional classes with more than one letter into Multi
        functional_class = ifelse(
            str_length(functional_class) >=2,
            "Multi",
            functional_class
        ),
        cog_annotation = ifelse(
            is.na(cog_annotation),
            "No annotation",
            cog_annotation
        )
    )



# Functional classification table:
# <class-id> <description>
# renamed class-id to funcional_class to match the cog_names table
functions <- read_tsv(
        file = opt$functions,
        col_names = c("functional_class", "description"),
        comment = "#"
    ) %>%
    rbind(
        c("Multi", "Multiple functions")
    )


# Table with the genome source
# <genome-code> <ncbi-taxid> <ncbi-ftp-name>
# Renamed ncbi-ftp-name to genome_name to match the cog table
genomes <- read_tsv(
        file = opt$genomes,
        col_names = c("genome_code", "ncbi_taxid", "genome_name"),
        comment = "#"
    )

# Table with the protein info
# Table schema:
# <protein-id> <refseq-acc>
proteins <- read_tsv(
        file = opt$proteins,
        col_names = c("protein_id", "refseq_acc"),
        comment = "#"
    )




# Merge functions and cog_names
cog_names_function <- full_join(
    x = cog_names,
    y = functions,
    by ="functional_class"
)

# Merge all tables with cog the main one
cog_full <- cog %>%
    left_join(x = ., y = genomes,            by = "genome_name") %>%
    left_join(x = ., y = cog_names_function, by = "cog_id") %>%
    left_join(x = ., y = proteins,           by = "protein_id")


write_tsv(
    x = cog_full,
    path = opt$output,
    col_names = TRUE
)

