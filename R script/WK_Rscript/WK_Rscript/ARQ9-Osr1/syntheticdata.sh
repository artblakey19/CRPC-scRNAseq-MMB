#!/bin/bash
for ((i=1; i<=20;i++))
do
    python3 chipulate-master/chipulate.py --input-file example/rep$i/test.tsv --genome-file chr1.fa --chrom-size-file chr1.fa.fai -d 50 --output-dir example/rep$i
    bedops -u example/rep$i/test.chip_reads.bed syn.bed > example/syn_with_peak$i.bed
    bedops -u example/rep$i/test.control_reads.bed syn.bed > example/syn2_with_peak$i.bed
    macs2 callpeak -t example/syn_with_peak$i.bed example/syn2_with_peak$i.bed -c 1.bed -f BED -n twosample -B -q 0.5 --outdir example/rep$i
    macs2 callpeak -t example/syn_with_peak$i.bed -f BED -n experimental -B -q 1 --outdir example/rep$i
    macs2 callpeak -t example/syn2_with_peak$i.bed -f BED -n experimental2 -B -q 1 --outdir example/rep$i
    macs2 callpeak -t 1.bed -f BED -n control -B -q 1 --outdir example/rep$i
done

