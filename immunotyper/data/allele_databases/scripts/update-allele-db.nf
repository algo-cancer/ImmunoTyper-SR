nextflow.enable.dsl=2

params.gene_type = 'IGHV'
params.novel_orphon_fasta = "$projectDir/NO_FILE"  // Default to NO_FILE
params.outdir = '../'  // Default to parent directory when running from scripts/
params.script_dir = "${workflow.projectDir}"
params.get_imgt_script = "${params.script_dir}/get_imgt_allele_database.sh"
params.remove_duplicates_script = "${params.script_dir}/remove_duplicate_sequences.py"
params.make_clusters_script = "${params.script_dir}/make_gene_clusters.py"
params.create_alignment_script = "${params.script_dir}/create_allele_db_alignment.py"


// Process to get IMGT allele database
process GET_IMGT_DB {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    val gene_type
    path novel_orphon_fasta

    output:
    path "${gene_type}-IMGT-allele-db.fa", emit: imgt_db

    script:
    def novel_orphon_arg = novel_orphon_fasta.name != 'NO_FILE' ? "--novel_orphon_fasta ${novel_orphon_fasta}" : ''
    """
    ${params.get_imgt_script} --gene-type ${gene_type} ${novel_orphon_arg} > ${gene_type}-IMGT-allele-db.fa
    """
}

// Process to remove duplicate sequences
process REMOVE_DUPLICATES {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path imgt_db

    output:
    path "${imgt_db.baseName}-no_duplicates.fa", emit: no_duplicates
    path "${imgt_db.baseName}-identical_sequences.txt"

    script:
    """
    python ${params.remove_duplicates_script} --input_fasta ${imgt_db} \
        --output_fasta ${imgt_db.baseName}-no_duplicates.fa \
        2> ${imgt_db.baseName}-identical_sequences.txt
    """
}

// Process to add flanking N nucleotides
process ADD_FLANKING_N {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path no_duplicates

    output:
    path "${no_duplicates.baseName}+Ns.fa", emit: with_ns

    script:
    """
    #!/usr/bin/env python3
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    
    # Add 100 N's to both ends of each sequence
    with open("${no_duplicates}", "r") as input_handle:
        with open("${no_duplicates.baseName}+Ns.fa", "w") as output_handle:
            for record in SeqIO.parse(input_handle, "fasta"):
                new_seq = "n" * 100 + str(record.seq) + "n" * 100
                new_record = SeqRecord(
                    Seq(new_seq),
                    id=record.id,
                    description=record.description
                )
                SeqIO.write(new_record, output_handle, "fasta")
    """
}

// Process to create gene clusters
process MAKE_GENE_CLUSTERS {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path imgt_db

    output:
    path "${params.gene_type}*gene_clusters.tsv"

    script:
    """
    python ${params.make_clusters_script} ${imgt_db}
    """
}

// Process to create BWA index
process BWA_INDEX {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path fasta

    output:
    path "${fasta}.*"

    script:
    """
    bwa index ${fasta}
    """
}

// Process to create Bowtie2 index
process BOWTIE_INDEX {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path fasta

    output:
    path "${fasta}.*"

    script:
    """
    bowtie2-build ${fasta} ${fasta}
    """
}

// Process to create allele database alignment
process CREATE_ALLELE_DB_ALIGNMENT {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path imgt_db

    output:
    path "${imgt_db.baseName}-aligned.fasta", emit: aligned_db
    path "${imgt_db.baseName}-consensus.fasta", emit: consensus_db

    script:
    """
    python ${params.create_alignment_script} --input_fasta ${imgt_db}
    """
}

workflow {
    // Create channel for optional input file
    novel_orphon_file = file(params.novel_orphon_fasta, checkIfExists: true)

    // Run processes
    imgt_db = GET_IMGT_DB(params.gene_type, novel_orphon_file)
    aligned_db = CREATE_ALLELE_DB_ALIGNMENT(imgt_db.imgt_db)
    no_duplicates = REMOVE_DUPLICATES(imgt_db.imgt_db)
    with_ns = ADD_FLANKING_N(no_duplicates.no_duplicates)
    
    // Run parallel processes
    MAKE_GENE_CLUSTERS(imgt_db.imgt_db)
    BWA_INDEX(with_ns.with_ns)
    BOWTIE_INDEX(with_ns.with_ns)
}
