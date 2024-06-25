#!/bin/bash

# 変数の初期化
NUM_THREADS=1
CPFASTA=""
MTFASTA=""
PROK_FLAG=false

SCRIPTDIR=$(cd $(dirname $0) && pwd)


# オプションの解析
while getopts "n:c:m:p" opt; do
  case $opt in
    n)
      NUM_THREADS="$OPTARG"
      ;;
    c)
      CPFASTA="$OPTARG"
      ;;
    m)
      MTFASTA="$OPTARG"
      ;;
    p)
      PROK_FLAG=true
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

# getoptsが処理したオプション引数を取り除く
shift $((OPTIND -1))

# 引数のチェック
if [ "$#" -ne 2 ]; then
    echo
    echo "Usage: $0 [-n num_threads] [-c cp.fasta] [-m mt.fasta] [-p] <query> <outdir>"
    echo 
    exit 1
fi

# unset this if you use seqkit in your PATH
SEQKIT_SIF="singularity exec /usr/local/biotools/s/seqkit:2.8.2--h9ee0642_0"


# 引数の割り当て
QUERY=$1
OUTDIR=$2
# CPFASTA=$3
# MTFASTA=$4


# クエリとして与えられたパスにファイルが存在するか確認
if [ ! -f "$QUERY" ]; then
    echo "Error: The file specified by query ($QUERY) does not exist."
    exit 1
else
    QUERY_FULLPATH=$(cd $(dirname $QUERY) && pwd)/$(basename $QUERY)
    echo INPUT QUERY FASTA: $QUERY_FULLPATH

fi

if [ -n "$CPFASTA" ] &&　[ ! -f "$CPFASTA" ]; then
    echo "Error: The file specified by query ($CPFASTA) does not exist."
    exit 1
else
    if [ -n "$CPFASTA" ]; then
        CPFASTA_FULLPATH=$(cd $(dirname $CPFASTA) && pwd)/$(basename $CPFASTA)
        echo Reference FASTA for chloroplast: $CPFASTA_FULLPATH
    fi
fi

if [ -n "$MTFASTA" ] &&　[ ! -f "$MTFASTA" ]; then
    echo "Error: The file specified by query ($MTFASTA) does not exist."
    exit 1
else
    if [ -n "$MTFASTA" ]; then
        MTFASTA_FULLPATH=$(cd $(dirname $MTFASTA) && pwd)/$(basename $MTFASTA)
        echo Reference FASTA for mitochondorion: $MTFASTA_FULLPATH
    fi
fi




# 出力ディレクトリの作成
mkdir -p "$OUTDIR"
OUTDIR_FULLPATH=$(cd $OUTDIR && pwd)


echo Number of threads: $NUM_THREADS
echo Enable BLAST against RefSeq_representative_prokaryotes: $PROK_FLAG


# ファイル名（拡張子付き）を取得
QUERY_BASENAME=$(basename "$QUERY")

# 拡張子を除いた部分と拡張子を取得
OUTPUT_BASE_NAME="${QUERY_BASENAME%.*}"
EXTENSION="${QUERY_BASENAME##*.}"

# 最終出力ファイル
OUTPUT_FASTA=${OUTDIR_FULLPATH}/${OUTPUT_BASE_NAME}.clean.${EXTENSION}

# echo $OUTPUT_FASTA


# 処理開始
cd $OUTDIR

# 1M以上とそれ未満に分ける
$SEQKIT_SIF seqkit seq -M 999999 $QUERY_FULLPATH > genome_lt1M.fa
$SEQKIT_SIF seqkit seq -m 10000000 $QUERY_FULLPATH > genome_1M.fa


cp genome_lt1M.fa tmp.genome.fa


if [ ! -z "$CPFASTA" ]; then
    echo ========== Removing short contigs derived from chloroplast ========== 
    python ${SCRIPTDIR}/remove_foreign_contigs.py -q tmp.genome.fa -s ${CPFASTA_FULLPATH} --use_singularity --prefix cp 
    mv tmp.genome.cp_filtered.fa tmp.genome.fa
fi

if [ ! -z "$MTFASTA" ]; then
    echo ========== Removing short contigs derived from mitochondorion ==========
    python ${SCRIPTDIR}/remove_foreign_contigs.py -q tmp.genome.fa -s ${MTFASTA_FULLPATH} --use_singularity --prefix mt 
    mv tmp.genome.mt_filtered.fa tmp.genome.fa
fi

echo ========== Removing short contigs derived from rDNA cluster region ==========
python ${SCRIPTDIR}/remove_rDNA_contigs.py --query tmp.genome.fa --use_singularity --num_threads $NUM_THREADS --prefix rDNA
mv tmp.genome.rDNA_filtered.fa tmp.genome.fa


if [ "$PROK_FLAG" == "true" ]; then
    echo ========== Removing short contigs derived from bacterial contamination ==========
    echo THIS STEP MAY TAKE TIME
    export SINGULARITY_BIND=/home/ddbjshare/blast/db/v5/
    python ${SCRIPTDIR}/remove_foreign_contigs.py -q tmp.genome.fa -d /home/ddbjshare/blast/db/v5/ref_prok_rep_genomes --num_threads $NUM_THREADS --prefix prok
    mv tmp.genome.prok_filtered.fa tmp.genome.fa
fi


mv tmp.genome.fa genome_lt1M.clean.fa
echo Writing output file to ${OUTPUT_FASTA}
cat genome_1M.fa genome_lt1M.clean.fa > $OUTPUT_FASTA
echo "Done."

echo ============================================================
$SEQKIT_SIF seqkit stats genome_1M.fa genome_lt1M.fa genome_lt1M.clean.fa ${OUTPUT_BASE_NAME}.clean.${EXTENSION}
