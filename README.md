# Cleanup_short_contigs

Currently, it is only internal use in NIG-SuperComputer system, since it uses singularity container in `/usr/loca/biotools/`.  

This script can remove short contigs derived from organellar genomes (chloroplast, mitochondrion), rDNA cluster regions, and bacterial contamination. It can also remove redundant short contigs overlapped with longer sequences.

__Use this script at your own risk.__

## Usage
```
cleanup_short_contigs [-n num_threads] [-c cp.fasta] [-m mt.fasta] [-r] [-p] [-s] <query.fa> <outdir>

```
The result will be written into `outdir/query.clean.fa`

The input query genomes will be split into long contigs (>1M) and short contigs (<1M). For the short contigs, the following processing is performed, and those that meet the criteria will be removed.


- Removal of chroplast sequences (`-c`)  
    Specify the pass to the FASTA file of the reference chloroplast genome with option `-c`  
    Contigs with more than 75% coverage will be removed.

- Removal of mitochondrial sequences (`-m`)  
    Specify the pass to the FASTA file of the reference mitochondrial genome with option `-m`  
    Contigs with more than 75% coverage will be removed.

- Removal of short contigs derived from rDNA cluster regionss (`-r`)  
    rDNA will be identified using Barrnap, and contigs with more than 45% of rDNA regions will be removed.

- Removal of short contigs derived from bacterial contamination (`-p`)  
    Queries will be searched against RefSeq representative prokaryotic genomes (`/home/ddbjshare/blast/db/v5/ref_prok_rep_genomes` in NIG-SC). This step may take time.  
    Contigs with more than 75% coverage will be removed.

- Removal of short contigs overlapped with longer sequences by self-BLAST (`-s`)  
	Short contigs (<1M) will be search against longer contigs (>1M), and contigs with more than 90% coverage will be removed.  
	Highly repetitive contigs may take time. If the running time exceeds 120s, the task will be killed and the query contig will be kept.  
To see which contigs have been removed, see `*_result.tsv` files for each step.

## Prerequisite
`Biopython` and `intervaltree` are required. Install them with 
```
mamba install -c bioconda biopython intervaltree
```

## Note
Since this script uses singularity, the query file must be in your home directory or somewhere that is accessible from the singularity container.
