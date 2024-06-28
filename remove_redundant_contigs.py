import argparse
import subprocess
import os
import glob
import concurrent.futures

blast_header = "qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp qcovs qcovus"

# Singularity container paths (must be updated regularly)
blast_sif = "/usr/local/biotools/b/blast:2.15.0--pl5321h6f7f691_1"
seqkit_sif = "/usr/local/biotools/s/seqkit:2.8.2--h9ee0642_0"
TIMEOUT=30

def parse_arg():
    parser = argparse.ArgumentParser(description='Run BLAST using Biopython.')
    parser.add_argument('--query', '-q', required=True, help='Query sequence file')
    subject_db_group = parser.add_mutually_exclusive_group(required=True)
    subject_db_group.add_argument('--subject', '-s', help='Subject sequence file')
    subject_db_group.add_argument('--db', '-d', help='Database')
    parser.add_argument('--evalue', type=float, default=1e-100, help='E-value threshold')
    parser.add_argument('--prefix', default="self", help='Output file prefix')
    parser.add_argument('--cov_cutoff', type=int, default=95, help='Qcov threshold (default=75)')
    parser.add_argument('--perc_identity', type=float, default=90, help='Percent identity (default=90)')
    parser.add_argument('--max_target_seqs', type=int, default=1, help='Maximum number of alignments (default=500)')
    parser.add_argument('--num_threads', '-n', type=int, default=1, help='Number of threads')
    parser.add_argument('--use_singularity', action='store_true', help='Use Singularity container')
    args = parser.parse_args()

    if not args.prefix:
        args.prefix = os.path.basename(args.subject).split('.')[0]

    return args

def make_blast_db(subject, use_singularity):
    if not os.path.exists(f"{subject}.nhr"):
        makeblastdb_cmd = ['makeblastdb', '-in', subject, '-dbtype', 'nucl']
        if use_singularity:
            makeblastdb_cmd = ['singularity', 'exec', blast_sif] + makeblastdb_cmd
        print(f"Running makeblastdb with command: {' '.join(makeblastdb_cmd)}")
        subprocess.run(makeblastdb_cmd)

def split_query(query, prefix, use_singularity):
    # split the query using seqkit split
    seqkit_cmd = ['seqkit', 'split', '-s', '1', '-O', "split", query]
    if use_singularity:
        seqkit_cmd = ['singularity', 'exec', seqkit_sif] + seqkit_cmd
    print(f"Running seqkit split with command: {' '.join(seqkit_cmd)}")
    subprocess.run(seqkit_cmd)
    split_file_list = glob.glob("split/*.part_*")
    return split_file_list

# def task_blast(query, prefix, subject, db, evalue, perc_identity, max_target_seqs, use_singularity):
#     try:
#         run_blastn(query, prefix, subject, db, evalue, perc_identity, max_target_seqs, use_singularity)
#         return f"Command: {command}\nOutput: {result.stdout.decode()}\nError: {result.stderr.decode()}"
#     except subprocess.CalledProcessError as e:
#         return f"Command '{command}' failed with error: {e.stderr.decode()}"


def run_blastn_multi(split_file_list, prefix, subject, db, evalue, perc_identity, max_target_seqs, num_threads, use_singularity):
    os.makedirs("split_blast", exist_ok=True)
    # for query in split_file_list:
    #     prefix = "split_blast/" + os.path.basename(query)
    #     run_blastn(query, prefix, subject, db, evalue, perc_identity, max_target_seqs, num_threads, use_singularity)
    blast_result_files = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
        # コマンドをスレッドプールに渡して実行
        futures = [executor.submit(run_blastn, query, "split_blast/" + os.path.basename(query), subject, db, evalue, perc_identity, max_target_seqs, use_singularity) 
        for query in split_file_list]  

        # 結果を取得して表示
        for future in concurrent.futures.as_completed(futures):
            blast_result_file, msg = future.result()
            if blast_result_file:
                blast_result_files.append(blast_result_file)
            else:
                print("!!!!! " + msg + " !!!!!")
    return blast_result_files


def run_blastn(query, prefix, subject, db, evalue, perc_identity, max_target_seqs, use_singularity):

    # Construct the BLAST command
    outfmt = f'7 {blast_header}'
    blast_cmd = [
        'blastn',
        '-query', query,
        '-evalue', str(evalue),
        '-perc_identity', str(perc_identity),
        '-max_target_seqs', str(max_target_seqs),
        '-outfmt', outfmt,
        '-out', f"{prefix}_blast.tab"
    ]

    if subject:
        blast_cmd.extend(['-subject', subject])
    if db:
        blast_cmd.extend(['-db', db])
    if use_singularity:
        blast_cmd = ['singularity', 'exec', blast_sif] + blast_cmd

    # Print the BLAST command (for debugging purposes)
    print(f"Running BLAST with command: {' '.join(blast_cmd)}")

    # # Run the BLAST command
    # subprocess.run(blast_cmd)

    try:
        result = subprocess.run(blast_cmd, shell=False, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=TIMEOUT)
        # return f"{prefix}_blast.tab", f"Command: {blast_cmd}\nOutput: {result.stdout.decode()}\nError: {result.stderr.decode()}"
        return f"{prefix}_blast.tab", f"success"

    except subprocess.CalledProcessError as e:
        return None, f"Failed: {e.stderr.decode()}"
    except subprocess.TimeoutExpired as e:
        return None, f"Timeout: '{query}' The query might be highly repetitive or the genome is too complex. Please check."

def parse_blast_result_multi(blast_result_files, prefix, cov_cutoff):
    out_file_hit_id = f"{prefix}_hit_id.txt"
    out_file_result = f"{prefix}_result.tsv"

    num_to_be_removed = 0
    with open(out_file_hit_id, 'w') as f, open(out_file_result, 'w') as f_result:
        print(f"Coverage > {cov_cutoff}% will be removed")
        print(f"Query\tTarget\tCoverage\tstatus")
        f_result.write(f"Query\tTarget\tCoverage\tstatus\n")
        for blast_result_file in blast_result_files:
            qaccver, qcovs = parse_blast_result(blast_result_file)
            status = "REMOVE" if qcovs > cov_cutoff else "KEEP"
            print(f"{qaccver}\t{prefix}\t{qcovs}\t{status}")
            f_result.write(f"{qaccver}\t{prefix}\t{qcovs}\t{status}\n")

            # print(f"{qaccver}\t{qcovs}")
            if qcovs > cov_cutoff:
                num_to_be_removed += 1
                f.write(qaccver + '\n')
    print(f"Writing seq IDs to remove: {out_file_hit_id}")
    return num_to_be_removed

def parse_blast_result(blast_result_file):
    blast_header_split = blast_header.split()
    qaccver_idx = blast_header_split.index('qaccver')
    qcovs_idx = blast_header_split.index('qcovs')

    with open(blast_result_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.strip().split('\t')
            qaccver = cols[qaccver_idx]
            qcovs = float(cols[qcovs_idx])
            return qaccver, qcovs


def remove_unwanted_sequences(query, prefix, use_singularity):
    out_file_hit_id = f"{prefix}_hit_id.txt"
    out_file_basename, extention = os.path.splitext(os.path.basename(query))
    out_file_filtered = f"{out_file_basename}.{prefix}_filtered{extention}"
    seqkit_cmd = [
        'seqkit', 'grep', '-v', '-r', '-n', '-f', out_file_hit_id, query, '-o', out_file_filtered
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

    make_blast_db(subject, use_singularity)
    split_file_list = split_query(query, prefix, use_singularity)
    # print(split_file_list)
    subject, db = None, subject
    blast_result_files = run_blastn_multi(split_file_list, prefix, subject, db, evalue, perc_identity, max_target_seqs, num_threads, use_singularity)
    num_to_be_removed = parse_blast_result_multi(blast_result_files, prefix, cov_cutoff)
    # run_blastn(query, prefix, subject, db, evalue, perc_identity, max_target_seqs, num_threads, use_singularity)
    # num_to_be_removed = parse_blast_result(prefix, cov_cutoff)
    if num_to_be_removed > 0:
        remove_unwanted_sequences(query, prefix, use_singularity)
    else:
        print("No sequences to be removed.")


