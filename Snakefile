shell.prefix("set -euo pipefail;")
configfile: "config.yaml"

snakefiles = "bin/snakefiles/"

include: snakefiles + "folders"
include: snakefiles + "raw"
include: snakefiles + "download"
include: snakefiles + "db"

include: snakefiles + "qc"
include: snakefiles + "filtering"
include: snakefiles + "tissue"

include: snakefiles + "transdecoder"
include: snakefiles + "swissprot"
include: snakefiles + "nrsubset"
include: snakefiles + "uniref90"
include: snakefiles + "pfam"
include: snakefiles + "busco"
include: snakefiles + "panther"
include: snakefiles + "cog"
include: snakefiles + "bioservices"
include: snakefiles + "go"
include: snakefiles + "reactome"

include: snakefiles + "ncrnas"
include: snakefiles + "trnas"




#include: snakefiles + "mirnas"
# include: snakefiles + "qdd"
# include: snakefiles + "interproscan"

rule all:
    input:
        expand(
            transdecoder + "transdecoder.{extension}",
            extension = "bed cds gff3 pep".split()
        ),
        swissprot + "blastp.tsv",
        #uniref90 + "blastp.tsv",
        nr + "blastp.tsv",
        pfam + "hmmscan.tsv",
        busco + "busco_figure.png", 
        ncrnas + "cmsearch.tsv",
        trnas + "trnascanse.tsv",
        #mirnas + "candidates.fa",
        bioservices + "query_swissprot.tsv",
        go +  "go_annotation.tsv",
        reactome + "reactome_top.pdf",
        panther + "pantherScore.tsv",
        #ipro         + "ipro.tsv",
        cog + "functions.pdf",
        tissue + "normalised_tpms.tsv"
        
       
        
rule help: 
    shell:
        """
        echo raw:
            < snakefiles/raw \
            fgrep -h "##" | sed 's/rule /\t/' | sed -e 's/\\$$//' | sed -e 's/##//' | sed 's/\: /:\t/'
        
        echo NR:
            < snakefiles/nr \
            fgrep -h "##" | sed 's/rule /\t/' | sed -e 's/\\$$//' | sed -e 's/##//' | sed 's/\: /:\t/'
        
        echo Uniref90:
            < snakefiles/uniref90 \
            fgrep -h "##" | sed 's/rule /\t/' | sed -e 's/\\$$//' | sed -e 's/##//' | sed 's/\: /:\t/'

        echo Busco \(lineage = arthropoda, eukaryota, metazoa, bacteria, fungi, vertebrata\):
            < snakefiles/busco \
            fgrep -h "##" | sed 's/rule /\t/' | sed -e 's/\\$$//' | sed -e 's/##//' | sed 's/\: /:\t/'
        
        echo COG:
            < snakefiles/cog \
            fgrep -h "##" | sed 's/rule /\t/' | sed -e 's/\\$$//' | sed -e 's/##//' | sed 's/\: /:\t/'
        
        echo TransDecoder:
            < snakefiles/transdecoder \
            fgrep -h "##" | sed 's/rule /\t/' | sed -e 's/\\$$//' | sed -e 's/##//' | sed 's/\: /:\t/'
        
        echo InterProScan:
            < snakefiles/interproscan \
            fgrep -h "##" | sed 's/rule /\t/' | sed -e 's/\\$$//' | sed -e 's/##//' | sed 's/\: /:\t/'
        
        echo miRNAs:
            < snakefiles/mirnas \
            fgrep -h "##" | sed 's/rule /\t/' | sed -e 's/\\$$//' | sed -e 's/##//' | sed 's/\: /:\t/'
        
        echo tRNAs:
            < snakefiles/trnas \
            fgrep -h "##" | sed 's/rule /\t/' | sed -e 's/\\$$//' | sed -e 's/##//' | sed 's/\: /:\t/'
        
        echo ncRNAs:
            < snakefiles/ncrnas \
            fgrep -h "##" | sed 's/rule /\t/' | sed -e 's/\\$$//' | sed -e 's/##//' | sed 's/\: /:\t/'
        """
