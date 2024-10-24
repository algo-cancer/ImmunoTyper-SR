#!/bin/bash

gene_types=("IGHV" "IGKV" "IGLV" "TRBV" "TRDV" "TRGV" "TRAV")

cd "immunotyper/data/allele_databases/"

for gene_type in "${gene_types[@]}"; do
    echo "Processing gene type: ${gene_type}"
    
    # Create directory for the gene type
    mkdir -p "${gene_type}"
    
    # Run nextflow command
    if [ "${gene_type}" == "IGHV" ]; then
        nextflow run scripts/update-allele-db.nf -resume \
            --gene_type "${gene_type}" \
            --outdir "./${gene_type}/new" \
            --novel_orphon_fasta "${gene_type}/IGHV-orphon-novel.fasta"
    else
        nextflow run scripts/update-allele-db.nf -resume \
            --gene_type "${gene_type}" \
            --outdir "./${gene_type}/new"
    fi
    
    # Check if nextflow command was successful
    if [ $? -eq 0 ]; then
        echo "Nextflow command completed successfully for ${gene_type}"
        
        # Compare the fasta headers only
        diff_output=$(diff \
            <(grep '^>' "${gene_type}/new/${gene_type}-IMGT-allele-db.fa" | cut -f2 -d'|' | sort) \
            <(grep '^>' "./${gene_type}/${gene_type}-IMGT-allele-db.fa" | cut -f2 -d'|' | sort))
        
        # Count sequences in both files
        new_count=$(grep -c '^>' "${gene_type}/new/${gene_type}-IMGT-allele-db.fa")
        old_count=$(grep -c '^>' "./${gene_type}/${gene_type}-IMGT-allele-db.fa")
        
        echo "Sequence counts for ${gene_type}:"
        echo "  New database: ${new_count} sequences"
        echo "  Old database: ${old_count} sequences"
        
        # Check if there are differences
        if [ -n "${diff_output}" ]; then
            echo "Differences found for ${gene_type}. Saving to ${gene_type}/${gene_type}-new-db-diffs.txt"
            echo "${diff_output}" > "${gene_type}/${gene_type}-new-db-diffs.txt"
            
            # Delete the work directory first
            echo "Deleting work directory..."
            rm -rf "${gene_type}/new/work"
            
            # Move remaining files to parent directory
            mv "${gene_type}/new"/* ./${gene_type}
            rmdir "${gene_type}/new"
        else
            echo "No differences found for ${gene_type}"
            # Remove the subdirectory and its contents since no differences were found
            rm -rf "${gene_type}/new"
        fi
    else
        echo "Nextflow command failed for ${gene_type}. Skipping comparison."
    fi
    
    echo "Finished processing ${gene_type}"
    echo "---"
done

# Check if any nextflow command failed
if [ $? -ne 0 ]; then
    echo "One or more nextflow commands failed. Exiting with error."
    exit 1
fi

echo "All gene types processed successfully."
