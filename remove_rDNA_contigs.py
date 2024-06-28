import argparse
import subprocess
import os
from Bio import SeqIO
from intervaltree import Interval, IntervalTree
from dataclasses import dataclass, field

blast_header = "qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp qcovs qcovus"

# Singularity container paths (must be updated regularly)
barrnap_sif = "/usr/local/biotools/b/barrnap:0.9--hdfd78af_4"
seqkit_sif = "/usr/local/biotools/s/seqkit:2.8.2--h9ee0642_0"


@dataclass
class BarrnapResult:
    query: str
    length: int
    hits: IntervalTree = field(default_factory=IntervalTree)

    # def initialize_tree(self):
        # self.hits = intervaltree.IntervalTree()

def parse_arg():
    parser = argparse.ArgumentParser(description='Run BLAST using Biopython.')
    parser.add_argument('--query', '-q', required=True, help='Query sequence file')
    parser.add_argument('--evalue', type=float, default=1e-10, help='E-value threshold')
    parser.add_argument('--prefix', default="rDNA", help='Output file prefix')
    parser.add_argument('--rdna_cov_cutoff', type=int, default=45, help='rDNA coverage threshold (default=50)')
    parser.add_argument('--num_threads', '-n', type=int, default=1, help='Number of threads')
    parser.add_argument('--use_singularity', action='store_true', help='Use Singularity container')
    args = parser.parse_args()

    return args


def run_barrnap(query, prefix, num_threads, evalue, use_singularity):
    # Construct the barrnap command
    barrnap_cmd = [
        'barrnap',
        '--threads', str(num_threads),
        '--kingdom', 'euk',
        '--reject', '0.1',
        '--evalue', str(evalue),
        '--quiet',
        query,
        ">", f"{prefix}_barrnap.gff"
    ]

    if use_singularity:
        barrnap_cmd = ['singularity', 'exec', barrnap_sif] + barrnap_cmd

    # Print the barrnap command (for debugging purposes)
    print(f"Running barrnap with command: {' '.join(barrnap_cmd)}")

    # Run the barrnap command
    subprocess.run(' '.join(barrnap_cmd), shell=True)

def parse_barrnap_gff(query, prefix, rdna_cov_cutoff):
    def get_result_dict(query):
        results = {}
        for record in SeqIO.parse(query, 'fasta'):
            length = len(record.seq)
            query = record.id
            result = BarrnapResult(query, length)
            # result.initialize_tree()
            results[query] = result
        return results

    results = get_result_dict(query)
    # print(results)
    with open(f"{prefix}_barrnap.gff") as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.strip().split('\t')
            query_id, start, end, evalue, attrs = cols[0], int(cols[3]) - 1, int(cols[4]), float(cols[5]), cols[8]
            name = attrs.split(';')[0].split('=')[1]  # Name=28S_rRNA;product=28S ribosomal RNA
            name = f"{name}:{start}..{end}"
            results[query_id].hits.add(Interval(start, end, name))


    out_file_hit_id = f"{prefix}_hit_id.txt"
    out_file_result = f"{prefix}_result.tsv"

    num_to_be_removed = 0
    with open(out_file_hit_id, 'w') as f, open(out_file_result, 'w') as f_result:
        print(f"Coverage > {rdna_cov_cutoff}% will be removed")
        print(f"Query\tTarget\tCoverage\tstatus\trDNA_length/contig_length")
        f_result.write(f"Query\tTarget\tCoverage\tstatus\trDNA_length/contig_length\n")
        for result in results.values():
            result.hits.merge_overlaps()
            rdna_length = sum([interval.length() for interval in result.hits])
            rdna_coverage = rdna_length / result.length * 100
            status = "REMOVE" if rdna_coverage > rdna_cov_cutoff else "KEEP"
            print(f"{result.query}\t{prefix}\t{rdna_coverage:.2f}\t{status}\t{rdna_length}/{result.length}")
            f_result.write(f"{result.query}\t{prefix}\t{rdna_coverage:.2f}\t{status}\t{rdna_length}/{result.length}\n")
            if rdna_coverage > rdna_cov_cutoff:
                num_to_be_removed += 1
                f.write(result.query + '\n')
    return num_to_be_removed

def remove_unwanted_sequences(query, prefix, use_singularity):
    out_file_hit_id = f"{prefix}_hit_id.txt"
    out_file_basename, extention = os.path.splitext(os.path.basename(query))
    out_file_filtered = f"{out_file_basename}.{prefix}_filtered{extention}"
    seqkit_cmd = [
        'seqkit', 'grep', '-v', '-n', '-r', '-f', out_file_hit_id, query, '-o', out_file_filtered
    ]
    if use_singularity:
        seqkit_cmd = ['singularity', 'exec', seqkit_sif] + seqkit_cmd
    print(f"Running seqkit grep with command: {' '.join(seqkit_cmd)}")
    subprocess.run(seqkit_cmd)
    print(f"Output file is written to {out_file_filtered}")

if __name__ == "__main__":
    args = parse_arg()
    query, prefix, num_threads, evalue, use_singularity, rdna_cov_cutoff = args.query, args.prefix, args.num_threads, args.evalue, args.use_singularity, args.rdna_cov_cutoff
    run_barrnap(query, prefix, num_threads, evalue, use_singularity)
    num_to_be_removed = parse_barrnap_gff(query, prefix, rdna_cov_cutoff)
    if num_to_be_removed > 0:
        remove_unwanted_sequences(query, prefix, use_singularity)
