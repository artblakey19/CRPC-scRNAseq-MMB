#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# PRJNA699369 (Duke, Chen et al. Clin Cancer Res 2022) -> 10x-style FASTQ
#
# The SRA submission is a CellRanger BAM converted to cSRA. The R1 read
# (barcode + UMI) was not archived as a read, but CellRanger's *corrected*
# barcode and UMI survive in the SEQUENCE.LINKAGE_GROUP column as
#     CB:<16bp>-1|UB:<10bp>
# so R1 can be reconstructed losslessly and the run re-processed with
# CellRanger exactly like our own samples.
#
# Reads without a CB (~2.5%) are dropped; CellRanger would discard them as
# invalid barcodes anyway.
#
# Usage: 01_sra_to_fastq.sh <sra_file> <outdir> <sample_name> [zthreads] [row_range]
#   row_range is optional and only for pilot runs, e.g. 1-2000000
# ---------------------------------------------------------------------------
set -euo pipefail

SRA="${1:?usage: 01_sra_to_fastq.sh <sra_file> <outdir> <sample_name> [zthreads] [row_range]}"
OUTDIR="${2:?missing outdir}"
SAMPLE="${3:?missing sample name}"
ZTHREADS="${4:-3}"
ROWRANGE="${5:-}"

mkdir -p "$OUTDIR"
R1="$OUTDIR/${SAMPLE}_S1_L001_R1_001.fastq.gz"
R2="$OUTDIR/${SAMPLE}_S1_L001_R2_001.fastq.gz"
STATS="$OUTDIR/${SAMPLE}.recon_stats.tsv"

RANGE_ARG=()
[[ -n "$ROWRANGE" ]] && RANGE_ARG=(-R "$ROWRANGE")

echo "[$(date +%H:%M:%S)] $SAMPLE : reconstructing FASTQ from $SRA"

vdb-dump "$SRA" -T SEQUENCE \
    -C "READ,(INSDC:quality:text:phred_33)QUALITY,LINKAGE_GROUP" \
    -f tab "${RANGE_ARG[@]}" 2>/dev/null \
| mawk -F'\t' -v R1="$R1" -v R2="$R2" -v ZT="$ZTHREADS" -v STATS="$STATS" '
BEGIN {
    r1c = "bgzip -@ " ZT " -c > \"" R1 "\""
    r2c = "bgzip -@ " ZT " -c > \"" R2 "\""
    # quality string donor for the synthetic R1 (barcode+UMI are already
    # error-corrected by CellRanger, so assign max Phred)
    QB  = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
    prevL = -1
}
{
    total++
    lg = $3
    if (substr(lg, 1, 3) != "CB:") { nocb++; next }
    p = index(lg, "|UB:")
    if (p == 0)                   { nocb++; next }

    cb = substr(lg, 4, p - 4)
    d  = index(cb, "-"); if (d > 0) cb = substr(cb, 1, d - 1)   # strip "-1"
    ub = substr(lg, p + 4)

    seq = $1; qual = $2
    if (length(seq) == 0 || length(seq) != length(qual)) { badlen++; next }

    bc = cb ub
    L  = length(bc)
    if (L != prevL) { prevL = L; qs = substr(QB, 1, L) }

    kept++
    nm = "@" total
    print nm " 1:N:0:1" | r1c
    print bc            | r1c
    print "+"           | r1c
    print qs            | r1c

    print nm " 2:N:0:1" | r2c
    print seq           | r2c
    print "+"           | r2c
    print qual          | r2c
}
END {
    close(r1c); close(r2c)
    printf "sample\t%s\n",       "'"$SAMPLE"'" >  STATS
    printf "total_spots\t%d\n",  total         >> STATS
    printf "kept\t%d\n",         kept          >> STATS
    printf "no_CB\t%d\n",        nocb          >> STATS
    printf "bad_length\t%d\n",   badlen        >> STATS
    printf "barcode_len\t%d\n",  prevL         >> STATS
}'

echo "[$(date +%H:%M:%S)] $SAMPLE : done"
cat "$STATS"
