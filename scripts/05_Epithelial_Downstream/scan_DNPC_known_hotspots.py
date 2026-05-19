#!/usr/bin/env python3
"""
Directed pileup at exact known DNPC/CRPC hotspot positions (GRCh38).
Reports base counts regardless of depth/VAF thresholds so we can see WT vs mutant calls
even when coverage is sparse.

Coordinates are GRCh38; for genes on the minus strand the REF/ALT shown
is the genomic-strand call (we still annotate the protein change for context).
"""
import os, sys
import pysam
from collections import Counter

BAM_DIR = "/home/MMB/projects/Human CRPC scRNAseq/Raw_data/CRPC_bam"
SAMPLES = ["CRPC1", "CRPC2", "CRPC3"]

# (gene, protein_change, chrom, pos_1based, ref_genomic, alt_genomic, note)
# Coordinates from ClinVar / COSMIC for canonical transcripts (GRCh38).
HOTSPOTS = [
    # AR (chrX, +) — LBD hotspots in exon 8; codon numbering is +1 vs old (T877A → T878A)
    ("AR",   "p.L702H",  "chrX",  67711613, "T", "A", "exon4 (LBD)"),
    ("AR",   "p.W742C",  "chrX",  67711734, "G", "T", "exon5"),
    ("AR",   "p.H875Y",  "chrX",  67723692, "C", "T", "exon8"),
    ("AR",   "p.F877L",  "chrX",  67723698, "T", "C", "exon8"),
    ("AR",   "p.T878A",  "chrX",  67723701, "A", "G", "exon8 — classic ARPC"),
    ("AR",   "p.T878S",  "chrX",  67723702, "C", "G", "exon8"),
    # FOXA1 (chr14, +) — forkhead hotspots
    ("FOXA1","p.R219S",  "chr14", 37592294, "C", "A", "forkhead"),
    ("FOXA1","p.F254L",  "chr14", 37592399, "T", "C", "wing2"),
    ("FOXA1","p.D226N",  "chr14", 37592315, "G", "A", "forkhead"),
    # PTEN (chr10, +)
    ("PTEN", "p.R130Q",  "chr10", 87933147, "G", "A", "phosphatase"),
    ("PTEN", "p.R130G",  "chr10", 87933147, "G", "C", "phosphatase"),
    ("PTEN", "p.R233*",  "chr10", 87957915, "C", "T", "C2 domain stop"),
    # TP53 (chr17, -) — REF/ALT shown on genomic strand
    ("TP53", "p.R175H",  "chr17",  7675088, "C", "T", "DBD"),
    ("TP53", "p.R248Q",  "chr17",  7674220, "C", "T", "DBD"),
    ("TP53", "p.R248W",  "chr17",  7674221, "G", "A", "DBD"),
    ("TP53", "p.R273H",  "chr17",  7673803, "C", "T", "DBD"),
    ("TP53", "p.R282W",  "chr17",  7673776, "G", "A", "DBD"),
    # SPOP (chr17, -)
    ("SPOP", "p.F133L",  "chr17", 49619095, "A", "G", "MATH"),
    ("SPOP", "p.Y87C",   "chr17", 49619140, "T", "C", "MATH"),
    # RB1 (chr13, +) — common nonsense in last exons
    ("RB1",  "p.R556*",  "chr13", 48381403, "C", "T", "pocket A"),
]

def query_position(bam, chrom, pos1):
    """Return Counter of base counts at 1-based pos."""
    pos0 = pos1 - 1
    cnt = Counter()
    for col in bam.pileup(chrom, pos0, pos0+1,
                          truncate=True,
                          min_base_quality=20,
                          min_mapping_quality=20,
                          ignore_overlaps=True,
                          stepper="samtools",
                          max_depth=200000):
        if col.reference_pos != pos0:
            continue
        for read in col.pileups:
            if read.is_del or read.is_refskip:
                continue
            qpos = read.query_position
            if qpos is None:
                continue
            cnt[read.alignment.query_sequence[qpos].upper()] += 1
    return cnt

def main():
    out_path = "/home/MMB/projects/Human CRPC scRNAseq/Results/DNPC_known_hotspots.tsv"
    header = ["sample","gene","protein","chrom","pos","ref","alt","note",
              "depth","ref_n","alt_n","other_n","alt_frac","base_counts"]
    rows = []
    for s in SAMPLES:
        bam = pysam.AlignmentFile(os.path.join(BAM_DIR, s, "possorted_genome_bam.bam"), "rb")
        sys.stderr.write(f"[{s}] querying {len(HOTSPOTS)} hotspots...\n")
        for gene, prot, chrom, pos, ref, alt, note in HOTSPOTS:
            cnt = query_position(bam, chrom, pos)
            depth = sum(cnt.values())
            ref_n = cnt.get(ref, 0)
            alt_n = cnt.get(alt, 0)
            other_n = depth - ref_n - alt_n
            alt_frac = round(alt_n/depth, 3) if depth else 0.0
            base_counts = ",".join(f"{b}:{n}" for b,n in cnt.most_common()) or "-"
            rows.append([s, gene, prot, chrom, pos, ref, alt, note,
                         depth, ref_n, alt_n, other_n, alt_frac, base_counts])
        bam.close()

    with open(out_path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    sys.stderr.write(f"output: {out_path}\n")

    # also pretty-print summary
    print(f"\n{'sample':6} {'gene':5} {'protein':10} {'note':25} {'depth':>5} {'ref':>4} {'alt':>4} {'altf':>5}  bases")
    print("-"*100)
    for r in rows:
        s,gene,prot,_,_,ref,alt,note,depth,ref_n,alt_n,_,af,bc = r
        flag = " *" if alt_n >= 3 and af >= 0.05 else ""
        print(f"{s:6} {gene:5} {prot:10} {note:25} {depth:5d} {ref_n:4d} {alt_n:4d} {af:5.2f}  {bc}{flag}")

if __name__ == "__main__":
    main()
