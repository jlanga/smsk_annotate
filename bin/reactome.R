require(readr)



annotate_reactome <- function(bioservices_filename, reactome_filename){
    require(readr)
    require(dplyr)
    
    u2r_colnames <- c(
        "Entry", "reactome_id", "url",
        "pathway_name", "evidence_code", "species"
    )
    u2r <- reactome_filename %>%
        read_tsv(col_names = u2r_colnames)
    
    up_acc <- bioservices_filename %>%
        read_tsv() %>%
        select(qseqid, Entry) %>%
        unique()

    annotation <- inner_join(up_acc, u2r, by = "Entry")

    return(annotation)    
}

plot_reactome_top <- function(annotation, top=25, filename_out){
    
    require(dplyr)
    require(ggplot2)
    top_reactome <- annotation %>%
        select(pathway_name) %>%
        group_by(pathway_name) %>%
        summarise(counts = n()) %>%
        arrange(desc(counts)) %>%
        slice(1:top) 
    
    top_reactome %>%
        ggplot(
            aes(x = pathway_name, y = counts)
        ) +
        geom_bar(stat = "identity") +
        ggtitle("Top Reactome Pathways") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 1))
    
    ggsave(file= filename_out, width = 297, height = 210, units = "mm")
}

compute_counts_reactome <- function(annotation){
    require(dplyr)
    
    annotation %>%
        select(pathway_name) %>%
        group_by(pathway_name) %>%
        summarise(counts = n()) %>%
        arrange(desc(counts)) %>%
        mutate(percent = counts/sum(counts) *100)
}


arguments <- commandArgs(trailingOnly = TRUE)
if(length(arguments) != 2){
    message("ERROR! Pass me only the bioservices file and the path to the Uniprot2Reactome_All_Levels.txt file")
    quit(save = "no", status = 1, runLast = FALSE)
}


bioservices_filename <- arguments[1]
reactome_filename <- arguments[2]

annotation <- annotate_reactome(
        bioservices_filename,
        reactome_filename
    )

write_tsv(annotation, "reactome_results.tsv")

plot_reactome_top(annotation, 25, "reactome_top.pdf")

compute_counts_reactome(annotation) %>%
    write_tsv("reactome_counts.tsv")


