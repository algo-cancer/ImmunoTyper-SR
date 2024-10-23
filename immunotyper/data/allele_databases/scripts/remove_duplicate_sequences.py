from __future__ import print_function
import sys, os
import argparse
from Bio import SeqIO
from src.common import fasta_from_seq, create_temp_file
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)
def oprint(*args, **kwargs):
    print(*args, file=sys.stdout, **kwargs)


def get_priority_score(description):
    """Returns priority score for a sequence description. Higher score = higher priority"""
    score = 0
    
    if '*01' in description:
        score += 4
    if '|F|' in description:
        score += 2
    if 'D*' not in description:
        score += 1
    return score

@pytest.fixture
def sample_records():
    """Fixture providing a set of test sequence records"""
    return [
        SeqRecord(Seq("ATGC"), id="seq1", description="HLA*01|F|"),
        SeqRecord(Seq("ATGC"), id="seq2", description="HLA*02"),
        SeqRecord(Seq("GCTA"), id="seq3", description="D*01"),
        SeqRecord(Seq("GCTA"), id="seq4", description="D*01|F|"),
    ]

@pytest.fixture
def temp_fasta(tmp_path, sample_records):
    """Fixture creating a temporary FASTA file with test sequences"""
    fasta_path = tmp_path / "test.fasta"
    with open(fasta_path, 'w') as f:
        for record in sample_records:
            f.write(f">{record.description}\n{str(record.seq)}\n")
    return str(fasta_path)

def test_get_priority_score():
    """Test priority score calculation for different sequence descriptions"""
    assert get_priority_score("HLA*01|F|") == 7  # *01 (4) + F (2) + not D* (1)
    assert get_priority_score("HLA*02") == 1     # not D* only
    assert get_priority_score("D*01|F|") == 6    # *01 (4) + F (2)
    assert get_priority_score("D*02") == 0       # no priority markers

def test_main_with_duplicates(sample_records, temp_fasta):
    """Test main function with duplicate sequences"""
    result = main(temp_fasta)
    
    # Should have 2 sequences (ATGC and GCTA)
    assert len(result) == 2
    
    # Check that highest priority sequences were kept
    sequences = {str(record.seq): record.description for record in result}
    assert sequences["ATGC"] == "HLA*01|F|"  # Higher priority due to *01, F, and not D*
    assert sequences["GCTA"] == "D*01|F|"    # Higher priority due to F flag

def test_main_empty_file(tmp_path):
    """Test handling of empty input file"""
    empty_fasta = tmp_path / "empty.fasta"
    empty_fasta.write_text("")
    
    result = main(str(empty_fasta))
    assert len(result) == 0

def test_main_single_sequences(tmp_path):
    """Test with sequences that have no duplicates"""
    unique_records = [
        SeqRecord(Seq("ATGC"), id="seq1", description="HLA*01"),
        SeqRecord(Seq("GCTA"), id="seq2", description="HLA*02"),
        SeqRecord(Seq("TAGA"), id="seq3", description="HLA*03"),
    ]
    
    fasta_path = tmp_path / "unique.fasta"
    with open(fasta_path, 'w') as f:
        for record in unique_records:
            f.write(f">{record.description}\n{str(record.seq)}\n")
    
    result = main(str(fasta_path))
    assert len(result) == 3
    assert sorted([str(r.seq) for r in result]) == sorted([str(r.seq) for r in unique_records])

def test_main_invalid_file():
    """Test handling of non-existent file"""
    with pytest.raises(FileNotFoundError):
        main("nonexistent_file.fasta")


def main(input_path):
    """Takes the path of a fasta as input, removes duplicate sequences based on priority rules"""
    
    input_seqs = list(SeqIO.parse(input_path, 'fasta'))
    sequence_mapping = dict()
    
    for record in input_seqs:
        seq = str(record.seq)
        if seq not in sequence_mapping:
            sequence_mapping[seq] = record
        else:
            existing = sequence_mapping[seq]
            existing_priority = get_priority_score(existing.description)
            current_priority = get_priority_score(record.description)
            
            eprint(f"Found duplicate sequence:")
            eprint(f"  Existing: {existing.description}")
            eprint(f"  Current:  {record.description}")
            
            if current_priority > existing_priority:
                eprint(f"  Keeping: {record.description} (higher priority)")
                sequence_mapping[seq] = record
            else:
                eprint(f"  Keeping: {existing.description} (first occurrence or equal/higher priority)")
    
    return list(sequence_mapping.values())

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_fasta', help="Path to fasta used as input to clustal", required=True)
    parser.add_argument('--output_fasta', help="Path to fasta used as input to clustal", required=False, default=None)
    args = parser.parse_args()

    dedupl_seqs = main(args.input_fasta)

    if args.output_fasta:
        eprint('Saving output to {}'.format(args.output_fasta))
        with open(args.output_fasta, 'w') as f:
            f.write(fasta_from_seq(*zip(*[(x.description, x.seq) for x in dedupl_seqs])))
    else:
        oprint('\n'.join([f">{x.description}\n{str(x.seq)}" for x in dedupl_seqs]))
