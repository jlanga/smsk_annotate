#!/usr/bin/env python3

# This script is supposed to take as input filtered blast (by evalue or whatever)
# And return a file with the full annotation from UniProt.
# The blast database used should be one from UniProt (Swissprot, UnirefXXX).

import sys
import pandas as pd
from bioservices import *



def read_blast_table(filename):
    """
    Function to open the Blast results (outfmt = 6).
    """
    blast = pd.read_table(filename, header=None)
    blast.columns = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    return blast



if __name__ == "__main__":
    
    if len(sys.argv) != 3:
        sys.exit("ERROR! Usage: annotate_with_bioservices_uniprot.py blast.tsv results.tsv")
    
    blast_in = sys.argv[1]
    results  = sys.argv[2]
    
    # Read the table
    raw_blast_table = read_blast_table(blast_in)
    
    # Drop all columns that are not qseqid and sseqid, and drop duplicated rows
    tidy_blast_table = raw_blast_table[["qseqid", "sseqid"]].drop_duplicates()

    # Add the Entry column: for each element in sseqid (Uniref90_XXXXX) select XXXXX
    tidy_blast_table["Entry"] = [x.split("_")[1] for x in tidy_blast_table["sseqid"]]
    
    # Create the connector to Uniprot
    u = UniProt(cache= True)
    u.debugLevel = "INFO"
    u.timeout = 10   # some queries are long and requires much more time; default is 1000 seconds
    
    # Query the UniProt database
    annotation = u.get_df(tidy_blast_table["sseqid"].tolist(), nChunk = 100)
    
    # Merge the blast table and the annotation table
    full_table = pd.merge(tidy_blast_table, annotation, on= "Entry", how = "left")
    
    # Save the results
    full_table.to_csv(results, sep= "\t", index= False)
    
