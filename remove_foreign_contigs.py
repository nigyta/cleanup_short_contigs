import argparse
import subprocess
import os

blast_header = "qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp qcovs qcovus"

# Singularity container paths (must be updated regularly)
blast_sif = "/usr/local/biotools/b/blast:2.15.0--pl5321h6f7f691_1"
seqkit_sif = "/usr/local/biotools/s/seqkit:2.8.2--h9ee0642_0"

def parse_arg():
    parser = argparse.ArgumentParser(description='Run BLAST using Biopython.')
    parser.add_argument('--query', '-q', required=True, help='Query sequence file')
    subject_db_group = parser.add_mutually_exclusive_group(required=True)
    subject_db_group.add_argument('--subject', '-s', help='Subject sequence file')
    subject_db_group.add_argument('--db', '-d', help='Database')
    parser.add_argument('--evalue', type=float, default=1e-10, help='E-value threshold')
    parser.add_argument('--prefix', default=None, help='Output file prefix')
    parser.add_argument('--cov_cutoff', type=int, default=75, help='Qcov threshold (default=75)')
    parser.add_argument('--perc_identity', type=float, default=90, help='Percent identity (default=90)')
    parser.add_argument('--max_target_seqs', type=int, default=500, help='Maximum number of alignments (default=500)')
    parser.add_argument('--num_threads', '-n', type=int, default=1, help='Number of threads')
    parser.add_argument('--use_singularity', action='store_true', help='Use Singularity container')
    args = parser.parse_args()

    if not args.prefix:
        args.prefix = os.path.basename(args.subject).split('.')[0]

    return args

def run_blastn(query, prefix, subject, db, evalue, perc_identity, max_target_seqs, num_threads, use_singularity):

    # Construct the BLAST command
    outfmt = f'7 {blast_header}'
    blast_cmd = [
        'blastn',
        '-query', query,
        '-evalue', str(evalue),
        '-perc_identity', str(perc_identity),
        '-max_target_seqs', str(max_target_seqs),
        '-outfmt', outfmt,
        '-out', f"{prefix}_blast.tab",
        '-num_threads', str(num_threads)
    ]

    if subject:
        blast_cmd.extend(['-subject', subject])
    if db:
        blast_cmd.extend(['-db', db])
    if use_singularity:
        blast_cmd = ['singularity', 'exec', blast_sif] + blast_cmd

    # Print the BLAST command (for debugging purposes)
    print(f"Running BLAST with command: {' '.join(blast_cmd)}")

    # Run the BLAST command
    subprocess.run(blast_cmd)

def parse_blast_result(prefix, cov_cutoff):
    blast_header_split = blast_header.split()
    qaccver_idx = blast_header_split.index('qaccver')
    qcovs_idx = blast_header_split.index('qcovs')
    result = {}
    with open(f"{prefix}_blast.tab") as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.strip().split('\t')
            qaccver = cols[qaccver_idx]
            qcovs = float(cols[qcovs_idx])
            qcovs_current = result.get(qaccver, 0)
            if qcovs > qcovs_current:
                result[qaccver] = qcovs

    out_file_hit_id = f"{prefix}_hit_id.txt"
    out_file_result = f"{prefix}_result.tsv"

    num_to_be_removed = 0
    with open(out_file_hit_id, 'w') as f, open(out_file_result, 'w') as f_result:
        print(f"Coverage > {cov_cutoff}% will be removed")
        print(f"Query\tTarget\tCoverage\tstatus")
        f_result.write(f"Query\tTarget\tCoverage\tstatus\n")
        for k, v in result.items():
            status = "REMOVE" if v > cov_cutoff else "KEEP"
            print(f"{k}\t{prefix}\t{v}\t{status}")
            f_result.write(f"{k}\t{prefix}\t{v}\t{status}\n")
            if v > cov_cutoff:
                num_to_be_removed += 1
                f.write(k + '\n')
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
    query, prefix, subject, db, evalue, num_threads = args.query, args.prefix, args.subject, args.db, args.evalue, args.num_threads
    cov_cutoff, use_singularity, perc_identity, max_target_seqs = args.cov_cutoff, args.use_singularity, args.perc_identity, args.max_target_seqs
    run_blastn(query, prefix, subject, db, evalue, perc_identity, max_target_seqs, num_threads, use_singularity)
    num_to_be_removed = parse_blast_result(prefix, cov_cutoff)
    if num_to_be_removed > 0:
        remove_unwanted_sequences(query, prefix, use_singularity)


