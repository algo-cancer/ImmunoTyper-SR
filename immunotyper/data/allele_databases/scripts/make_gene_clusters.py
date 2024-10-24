#!/usr/bin/env python

import argparse
import os
import logging
from collections import defaultdict
from itertools import combinations
import pytest
from Bio import SeqIO
from io import StringIO

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def extract_allele_id(record_id, delimiter, key):
    """
    Extracts the allele ID from a record ID using the specified delimiter and key.
    If the key is out of range, it defaults to using the first field.
    """
    fields = record_id.split(delimiter)
    try:
        return fields[key]
    except IndexError:
        logging.error(f"Key {key} is out of range for record ID: {record_id}. Using default key 0.")
        return fields[0]

def get_database(fasta_file, delimiter, key):
    logging.info(f"Reading FASTA file: {fasta_file}")
    genes = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        allele_id = extract_allele_id(record.id, delimiter, key)
        gene_id = allele_id.split('*')[0]
        genes.setdefault(gene_id, []).append(record)
    logging.info(f"Found {len(genes)} unique alleles")
    return genes

def cluster_genes(genes, delimiter, key):
    logging.info("Clustering genes")
    same_cluster = {}
    
    # Iterate through all pairs of gene groups
    for g1, g2 in combinations(genes.values(), 2):
        # Flag to track if we found identical sequences between these groups
        found_identical = False
        
        # Compare each sequence from first gene group
        for allele1 in g1:
            # Compare with each sequence from second gene group
            for allele2 in g2:
                # Count differences between sequences
                print
                differences = 0
                for c1, c2 in zip(allele1.seq, allele2.seq):
                    if c1 != c2:
                        differences += 1
                
                # If sequences are identical (no differences)
                if differences == 0:
                    found_identical = True
                    break
            if found_identical:
                break
        
        # If we found identical sequences, merge the clusters
        if found_identical:
            g1_id = extract_allele_id(g1[0].id, delimiter, key).split('*')[0]
            g2_id = extract_allele_id(g2[0].id, delimiter, key).split('*')[0]
            g1_c = same_cluster.get(g1_id, set())
            g2_c = same_cluster.get(g2_id, set())
            new_cluster = g1_c.union(g2_c).union({g1_id, g2_id})
            for a in new_cluster:
                same_cluster[a] = frozenset(new_cluster)
                
    logging.info(f"Generated {len(set(same_cluster.values()))} clusters")
    return set(same_cluster.values())
def main(fasta_file, output_file, delimiter, key):
    logging.info("Starting gene clustering process")
    genes = get_database(fasta_file, delimiter, key)
    clusters = cluster_genes(genes, delimiter, key)

    logging.info(f"Writing clusters to output file: {output_file}")
    with open(output_file, "w") as f:
        for cl in clusters:
            f.write("\t".join(cl) + "\n")
    logging.info("Gene clustering process completed")

# Unit tests
def test_cluster_genes_identical_sequences():
    fasta_data = """>gene1*allele1
ATGC
>gene2*allele2
ATGC
"""
    genes = get_database(StringIO(fasta_data), "|", 1)
    clusters = cluster_genes(genes, "|", 1)
    assert len(clusters) == 1
    assert {"gene1", "gene2"} in clusters

def test_cluster_genes_different_sequences():
    fasta_data = """>gene1*allele1
ATGC
>gene2*allele2
ATGA
"""
    genes = get_database(StringIO(fasta_data), "|", 1)
    clusters = cluster_genes(genes, "|", 1)
    assert len(clusters) == 0

def test_cluster_genes_complex_case():
    fasta_data = """>gene1*allele1
ATGC
>gene1*allele2
ATGA
>gene1*allele3
ATGG
>gene2*allele1
ATGC
>gene2*allele2
ATGA
>gene3*allele1
CGTA
>gene3*allele2
CGTC
>gene3*allele3
CGTG
>gene3*allele4
CGTT
"""
    genes = get_database(StringIO(fasta_data), "|", 1)
    clusters = cluster_genes(genes, "|", 1)
    
    # Verify there is exactly one cluster
    assert len(clusters) == 1
    
    # Verify the cluster contains gene1 and gene2 only
    cluster = clusters.pop()
    assert cluster == {"gene1", "gene2"}
    
    # Additional assertions to verify the gene groupings
    assert len(genes["gene1"]) == 3  # gene1 has 3 sequences
    assert len(genes["gene2"]) == 2  # gene2 has 2 sequences
    assert len(genes["gene3"]) == 4  # gene3 has 4 sequences

def test_cluster_genes_three_way_cluster():
    fasta_data = """>gene1*allele1
ATGC
>gene1*allele2
ATGA
>gene2*allele1
ATGC
>gene2*allele2
CGTA
>gene3*allele1
CGTA
>gene3*allele2
ATGA
"""
    genes = get_database(StringIO(fasta_data), "|", 1)
    clusters = cluster_genes(genes, "|", 1)
    
    # Verify there is exactly one cluster
    assert len(clusters) == 1
    
    # Verify the cluster contains all three genes
    cluster = clusters.pop()
    assert cluster == {"gene1", "gene2", "gene3"}
    
    # Additional assertions to verify the gene groupings
    assert len(genes["gene1"]) == 2  # gene1 has 2 sequences: ATGC, ATGA
    assert len(genes["gene2"]) == 2  # gene2 has 2 sequences: ATGC, CGTA
    assert len(genes["gene3"]) == 2  # gene3 has 2 sequences: CGTA, ATGA

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cluster genes in IGHV allele database if they share allele sequences")
    parser.add_argument("fasta_file", help="Input FASTA file")
    parser.add_argument("-o", "--output", help="Output TSV file", default=None)
    parser.add_argument("-d", "--delimiter", help="Delimiter between fields in the FASTA ID", default="|")
    parser.add_argument("-k", "--key", type=int, help="Key integer specifying the delimited column containing the allele ID", default=1)
    args = parser.parse_args()

    if args.output is None:
        output_file = f"{os.path.splitext(args.fasta_file)[0]}-gene_clusters.tsv"
    else:
        output_file = args.output

    main(args.fasta_file, output_file, args.delimiter, args.key)
