# Cleanup_short_contigs

Currently, it is only internal use in NIG-SuperComputer system, since it uses singularity container in `/usr/loca/biotools/`  

This script can remove short contigs derived from organellar genomes (chloroplast, mitochondrion) or rDNA cluster regions.  
Reference genomes for chloroplast and mitochondrion must be specified with options `-c` and `-p`, respectively. If not specified, the step will be skipped.  

Optionally, it can remove contigs derived from bacterial contamination (with `-p` flag), which internally runs BLAST against representative genomes of Refseq Prokaryote. So it may take time. 

## Usage
```
cleanup_short_contigs [-n num_threads] [-c cp.fasta] [-m mt.fasta] [-p] <query> <outdir>
```
