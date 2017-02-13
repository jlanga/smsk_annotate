#!/usr/bin/env Rscript

library("optparse")
library("tidyverse")
library("stringr")


option_list = list(
    make_option(
        opt_str = c("-i", "--input"),
        type = "character",
        default = NULL, 
        help = "TSV with the fully annotated results",
        metavar = "character"
    ),
    make_option(
        opt_str = c("-o", "--output"),
        type = "character",
        default = NULL,
        help = "A PDF file where to store the results",
        metavar = "character"
    )
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# opt$input = "results/cog/cog_results.tsv"
# opt$out = "results/cog/functions.pdf"



cog_full <- read_tsv(opt$input)

functions <- cog_full %>%
    select(functional_class, description) %>%
    unique() %>%
    arrange(functional_class)



# Get a dataframe with the data we want: the counts of each category; properly
# formatted for ggplot
compute_counts <- function(cog_analysis){
    cog_analysis %>%
        group_by(functional_class) %>%
        summarise(number = n()) %>%
        full_join(y = functions, by = "functional_class") %>%
        na.exclude() %>%
        mutate(
            legend_description = str_c(
                functional_class, ": ", description,
                sep = "",
                collapse = NULL)
        )
}


# Do the plot
plot_cog_counts <- function(cog_counts, file_out = "plot.pdf", title = "COG Classification"){
    require(ggplot2)
    q <- ggplot(
        data = cog_counts,
        aes(
            x = functional_class,
            y = number,
            fill = functional_class
        )
    ) +
        guides(
            fill = guide_legend(ncol=1)
        ) +
        geom_bar(stat = "identity") +
        labs(
            title = title,
            x = "Functional Classification",
            y = "Number of Isoforms"
        ) +
        scale_fill_manual(
            values = rep("black", nrow(cog_counts)),
            name = "Functional classification",
            breaks = cog_counts$functional_class,
            labels = cog_counts$legend_description
        )
    q
    
    ggsave(filename = file_out, width = 29.7, height = 21.0, units = "cm")
}


read_tsv(file= opt$input) %>%
    compute_counts() %>%
    plot_cog_counts(file_out = opt$out)
    


