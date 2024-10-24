#!/bin/bash

# Function to display help message
show_help() {
    echo "Usage: $0 --gene-type <GENE_TYPE> [--novel-orphon-fasta <FASTA_FILE>]"
    echo
    echo "Required arguments:"
    echo "  --gene-type           Specify the gene type (e.g., IGHV, IGKV, IGLV)"
    echo
    echo "Optional arguments:"
    echo "  --novel-orphon-fasta  Path to a FASTA file containing novel orphon sequences"
    exit 1
}

# Initialize variables
GENE_TYPE=""
NOVEL_ORPHON_FASTA=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --novel-orphon-fasta)
            NOVEL_ORPHON_FASTA="$2"
            shift 2
            ;;
        --gene-type)
            GENE_TYPE="$2"
            shift 2
            ;;
        -h|--help)
            show_help
            ;;
        *)
            shift
            ;;
    esac
done

# Check if gene type is provided
if [ -z "$GENE_TYPE" ]; then
    echo "Error: --gene-type is required"
    show_help
fi

# Get IMGT database sequences
wget -qO - "https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+allP" | perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' | grep -A1 "$GENE_TYPE" | grep 'Homo sapiens' -A1 | sed 's/>\([^ ]*\) />\1_/g'

# If novel orphon fasta is provided, append its contents
if [ ! -z "$NOVEL_ORPHON_FASTA" ]; then
    if [ -f "$NOVEL_ORPHON_FASTA" ]; then
        cat "$NOVEL_ORPHON_FASTA"
    else
        echo "Error: Novel orphon FASTA file not found: $NOVEL_ORPHON_FASTA" >&2
        exit 1
    fi
fi
