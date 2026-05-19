#!/usr/bin/env python3
"""
Reference-free pileup scan of DNPC-relevant hotspot regions across the three
CRPC 10x BAMs (3' v3 chemistry — so coverage is limited to last exons / 3' UTRs).

Strategy: tally A/C/G/T per genomic position, infer REF = majority base, and
emit positions where the second-most-common base reaches a configurable
allele-fraction threshold AND minimum alt-read count.
"""
import os
import sys
import pysam
from collections import Counter

BAM_DIR = "/home/MMB/projects/Human CRPC scRNAseq/Raw_data/CRPC_bam"
SAMPLES = ["CRPC1", "CRPC2", "CRPC3"]

# (label, chrom, start, end)  — GRCh38, regions chosen to cover 3' / last exons
REGIONS = [
    ("AR_LBD_ex4-8",     "chrX",  67711000, 67730619),   # AR LBD: L702H / W742C / H875Y / F877L / T878A
    ("FOXA1_full",       "chr14", 37589551, 37595011),
    ("TP53_ex7-11",      "chr17",  7670600,  7676600),   # minus strand; ex11 ≈ 7669600–7670600, scan a bit wider
    ("PTEN_ex7-9",       "chr10", 87952000, 87971930),
    ("RB1_ex22-27",      "chr13", 48470000, 48481890),
    ("MYC_full",         "chr8", 127735434,127741434),
    ("SPOP_MATH",        "chr17", 49609000, 49619594),
]

MIN_DEPTH   = 20
MIN_ALT_FRAC = 0.10
MIN_ALT_CNT  = 4
MIN_BASE_Q   = 20

def scan_bam(bam_path, sample):
    bam = pysam.AlignmentFile(bam_path, "rb")
    hits = []
    for label, chrom, start, end in REGIONS:
        for col in bam.pileup(chrom, start, end,
                              truncate=True,
                              min_base_quality=MIN_BASE_Q,
                              min_mapping_quality=20,
                              ignore_overlaps=True,
                              stepper="samtools",
                              max_depth=100000):
            bases = []
            for read in col.pileups:
                if read.is_del or read.is_refskip:
                    continue
                qpos = read.query_position
                if qpos is None:
                    continue
                bases.append(read.alignment.query_sequence[qpos].upper())
            if not bases:
                continue
            depth = len(bases)
            if depth < MIN_DEPTH:
                continue
            cnt = Counter(bases)
            if len(cnt) < 2:
                continue
            top = cnt.most_common()
            ref_base, ref_cnt = top[0]
            alt_base, alt_cnt = top[1]
            alt_frac = alt_cnt / depth
            if alt_frac < MIN_ALT_FRAC or alt_cnt < MIN_ALT_CNT:
                continue
            hits.append({
                "sample": sample, "region": label, "chrom": chrom,
                "pos": col.reference_pos + 1,  # 1-based
                "depth": depth, "ref": ref_base, "ref_n": ref_cnt,
                "alt": alt_base, "alt_n": alt_cnt, "alt_frac": round(alt_frac, 3),
                "all": ",".join(f"{b}:{n}" for b,n in top),
            })
    bam.close()
    return hits

def main():
    out_path = "/home/MMB/projects/Human CRPC scRNAseq/Results/DNPC_hotspot_pileup.tsv"
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    header = ["sample","region","chrom","pos","depth","ref","ref_n","alt","alt_n","alt_frac","all"]
    n_hits_total = 0
    with open(out_path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for s in SAMPLES:
            bam = os.path.join(BAM_DIR, s, "possorted_genome_bam.bam")
            sys.stderr.write(f"[{s}] scanning...\n")
            hits = scan_bam(bam, s)
            sys.stderr.write(f"[{s}] {len(hits)} candidate sites\n")
            n_hits_total += len(hits)
            for h in hits:
                fh.write("\t".join(str(h[k]) for k in header) + "\n")
    sys.stderr.write(f"total candidate sites: {n_hits_total}\n")
    sys.stderr.write(f"output: {out_path}\n")

if __name__ == "__main__":
    main()
