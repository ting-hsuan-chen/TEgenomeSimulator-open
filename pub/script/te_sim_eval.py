#!/usr/bin/env python3
"""
te_sim_eval.py

Comprehensive evaluation scripts for TE simulator outputs vs RepeatMasker results.

Inputs:
 - --genome         : whole sequence FASTA (single-seq)
 - --te_gff         : GFF file containing ground-truth TE intervals
 - --nonte          : optional non-TE FASTA (auto-generated if not provided)
 - --te_only        : optional TE-only FASTA (auto-generated if not provided)
 - --te_fasta       : optional multi-FASTA with one entry per TE insertion (auto-generated if not provided)
 - --rm_out         : RepeatMasker .out file (from whole sequence)

Outputs: printed summary + optional CSVs and PNGs saved to output directory.

Dependencies: biopython, numpy, scipy, scikit-learn, pandas, matplotlib
Install: pip install biopython numpy scipy scikit-learn pandas matplotlib

Usage:
python te_sim_eval.py \
   --genome genome.fa \
   --te_gff truth.gff --rm_out genome.out \
   --outdir results --k 6 --win 1000 --eval all
"""

import argparse
import os
import sys
import math
import gzip
import lzma
from collections import Counter, defaultdict
from itertools import product

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from scipy.stats import entropy
from sklearn.metrics import precision_recall_curve, auc
import matplotlib.pyplot as plt

# ---- Utility functions ----
def read_single_fasta(path):
    records = list(SeqIO.parse(path, 'fasta'))
    if len(records) != 1:
        print(f"Warning: {path} contains {len(records)} sequences (expect 1). Using first.")
    return str(records[0].seq).upper()

def read_multifasta_sequences(path):
    return {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(path, 'fasta')}

def write_fasta(records, path):
    SeqIO.write(records, path, 'fasta')

# ---- GFF parsing ----
def parse_gff_intervals(gff_file):
    intervals = []
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue
            start, end = int(fields[3]), int(fields[4])
            intervals.append({'start': start, 'end': end})
    return intervals

# ---- Generate derived FASTAs ----
def generate_te_and_non_te(genome_seq, gff_intervals, out_te, out_nonte):
    mask = np.zeros(len(genome_seq), dtype=bool)
    for iv in gff_intervals:
        mask[iv['start']-1:iv['end']] = True
    te_seq = ''.join([base if mask[i] else '' for i, base in enumerate(genome_seq)])
    nonte_seq = ''.join([base if not mask[i] else '' for i, base in enumerate(genome_seq)])
    write_fasta([SeqRecord(Seq(te_seq), id='TE_only', description='')], out_te)
    write_fasta([SeqRecord(Seq(nonte_seq), id='nonTE', description='')], out_nonte)
    return te_seq, nonte_seq

def generate_te_copies(genome_seq, gff_intervals, out_te_fasta):
    records = []
    for i, iv in enumerate(gff_intervals, 1):
        seq = genome_seq[iv['start']-1:iv['end']]
        records.append(SeqRecord(Seq(seq), id=f'TE_{i}', description=f"{iv['start']}-{iv['end']}"))
    write_fasta(records, out_te_fasta)
    return {rec.id: str(rec.seq) for rec in records}


# ---------- IO helpers ----------
def read_single_fasta(path):
    recs = list(SeqIO.parse(path, 'fasta'))
    if len(recs) != 1:
        print(f"Warning: {path} contains {len(recs)} sequences (expect 1). Using first.")
    return str(recs[0].seq).upper()

def read_multifasta_sequences(path):
    return {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(path, 'fasta')}

# ---------- k-mer utilities ----------
def canonical_kmer(kmer):
    rc = kmer.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    return kmer if kmer <= rc else rc


def count_kmers(seq, k=6, canonical=True):
    counts = Counter()
    seq = seq.replace('N','')
    n = len(seq) - k + 1
    if n <= 0:
        return counts
    for i in range(n):
        kmer = seq[i:i+k]
        if 'N' in kmer:
            continue
        if canonical:
            kmer = canonical_kmer(kmer)
        counts[kmer] += 1
    return counts


def kmer_vector(counts, all_kmers):
    return np.array([counts.get(k, 0) for k in all_kmers], dtype=float)


# ---------- similarity metrics ----------
def cosine_similarity(a, b):
    num = np.dot(a, b)
    den = np.linalg.norm(a) * np.linalg.norm(b)
    if den == 0:
        return 0.0
    return float(num/den)


def spectral_angle_from_cosine(cosine):
    c = np.clip(cosine, -1.0, 1.0)
    return float(math.acos(c))


# ---------- entropy & sliding windows ----------
def shannon_entropy_bases(seq):
    seq = seq.upper()
    counts = Counter(b for b in seq if b in 'ACGT')
    total = sum(counts.values())
    if total == 0:
        return 0.0
    ent = 0.0
    for b in 'ACGT':
        p = counts.get(b, 0) / total
        if p > 0:
            ent -= p * math.log2(p)
    return ent


def kmer_entropy(seq, k=3):
    seq = seq.upper()
    n = len(seq) - k + 1
    if n <= 0:
        return 0.0
    counts = Counter(seq[i:i+k] for i in range(n) if 'N' not in seq[i:i+k])
    total = sum(counts.values())
    ent = -sum((v/total) * math.log2(v/total) for v in counts.values())
    return ent


def sliding_entropy(seq, win=1000, step=None, k=None):
    if step is None:
        step = win
    positions = []
    values = []
    L = len(seq)
    for start in range(0, L, step):
        end = min(start + win, L)
        window = seq[start:end]
        if k is None:
            ent = shannon_entropy_bases(window)
        else:
            ent = kmer_entropy(window, k=k)
        positions.append(start)
        values.append(ent)
        if end == L:
            break
    return np.array(positions), np.array(values)


# ---- Compressibility ----
def compressed_size_bytes(text, method='gzip'):
    data = text.encode()
    if method == 'gzip':
        return len(gzip.compress(data))
    elif method == 'xz':
        return len(lzma.compress(data))
    else:
        raise ValueError('unknown compression method')

# ---- Reciprocal overlap ----
def reciprocal_overlap(start1, end1, start2, end2):
    inter_start = max(start1, start2)
    inter_end = min(end1, end2)
    if inter_start >= inter_end:
        return 0.0
    inter_len = inter_end - inter_start + 1
    len1 = end1 - start1 + 1
    len2 = end2 - start2 + 1
    return inter_len / min(len1, len2)

# ---------- RepeatMasker .out parser (simple) ----------
# This parser expects the classic RepeatMasker .out alignment table with a header. It will
# attempt to extract columns: score, perc.div, perc.del, perc.ins, sequence, begin, end, (left),
# repeat, class/family, rbegin, rend, (rleft). Coordinates converted to 0-based half-open [start,end).

def parse_repeatmasker_out(path):
    preds = []
    with open(path) as fh:
        lines = fh.readlines()
    # skip header lines until one that starts with 'score' or with a number column
    start_idx = 0
    for i,l in enumerate(lines[:50]):
        if l.strip().startswith('SW') or l.strip().startswith('score') or l.strip().startswith('scores'):
            start_idx = i+1
            break
    for l in lines[start_idx:]:
        if not l.strip():
            continue
        # columns separated by whitespace but repeat name and class/family may have spaces; typical RM uses fixed columns
        parts = l.rstrip('\n')
        if len(parts) < 20:
            continue
        # best-effort: split into at most 15 fields
        cols = parts.split()
        try:
            score = float(cols[0])
            perc_div = float(cols[1])
            # query seq name
            qname = cols[4]
            qbegin = int(cols[5])
            qend = int(cols[6])
            # matching repeat name and class are often in cols[9] and cols[10]
            # reconstruct repeat name/class by taking remaining columns
            repeat_name = cols[9]
            class_family = cols[10] if len(cols) > 10 else ''
            # Convert to 0-based half-open
            start = qbegin - 1
            end = qend
            preds.append({'score': score, 'div': perc_div, 'name': repeat_name, 'class': class_family, 'start': start, 'end': end})
        except Exception:
            continue
    return preds


# ---------- GFF parser (TE truth) ----------
def parse_gff_intervals(path, seqname=None):
    intervals = []
    with open(path) as fh:
        for l in fh:
            if l.startswith('#'):
                continue
            cols = l.rstrip('\n').split('\t')
            if len(cols) < 9:
                continue
            chrom = cols[0]
            s = int(cols[3]) - 1  # GFF is 1-based inclusive
            e = int(cols[4])      # keep half-open
            attr = cols[8]
            if seqname and chrom != seqname:
                continue
            intervals.append({'chrom': chrom, 'start': s, 'end': e, 'attr': attr})
    return intervals


# ---- RepeatMasker recovery evaluation ----
# This parser expects the classic RepeatMasker .out alignment table with a header. It will
# attempt to extract columns: score, perc.div, perc.del, perc.ins, sequence, begin, end, (left),
# repeat, class/family, rbegin, rend, (rleft). Coordinates converted to 0-based half-open [start,end).
def evaluate_recovery(truth_intervals, pred_intervals, score_field='score', thresholds=None, ro_th=0.5):
    if thresholds is None:
        scores = [p.get(score_field, 0.0) for p in pred_intervals]
        thresholds = sorted(set([0.0] + [float(s) for s in scores]), reverse=True)

    results = []
    for th in thresholds:
        TP = 0
        FP = 0
        FN = 0
        truth_used = [False] * len(truth_intervals)

        for p in pred_intervals:
            if p.get(score_field, 0.0) < th:
                continue
            matched = False
            for i, t in enumerate(truth_intervals):
                if truth_used[i]:
                    continue
                ro = reciprocal_overlap(p['start'], p['end'], t['start'], t['end'])
                if ro >= ro_th:
                    TP += 1
                    truth_used[i] = True
                    matched = True
                    break
            if not matched:
                FP += 1

        FN = truth_used.count(False)
        precision = TP / (TP + FP) if TP + FP > 0 else 0.0
        recall = TP / (TP + FN) if TP + FN > 0 else 0.0
        results.append({'threshold': th, 'TP': TP, 'FP': FP, 'FN': FN, 'precision': precision, 'recall': recall})

    return pd.DataFrame(results)

# ---------- Main pipeline ----------
def main():
    p = argparse.ArgumentParser()
    p.add_argument('--genome', required=True)
    p.add_argument('--te_gff', required=True)
    p.add_argument('--nonte')
    p.add_argument('--te_only')
    p.add_argument('--te_fasta')
    p.add_argument('--rm_out')
    p.add_argument('--outdir', default='results')
    p.add_argument('--k', type=int, default=6)
    p.add_argument('--win', type=int, default=1000)
    p.add_argument('--step', type=int, default=None)
    p.add_argument('--ro', type=float, default=0.5)
    p.add_argument('--eval', choices=['kmer','entropy','compress','rm','all'], default='all')
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Load genome and truth
    genome = read_single_fasta(args.genome)
    truth = parse_gff_intervals(args.te_gff)

    # Generate nonTE and TE-only if not provided
    nonte_path = args.nonte if args.nonte else os.path.join(args.outdir, 'nonte.fa')
    te_only_path = args.te_only if args.te_only else os.path.join(args.outdir, 'te_only.fa')
    if not os.path.exists(nonte_path) or not os.path.exists(te_only_path):
        te_only, nonte = generate_te_and_non_te(genome, truth, te_only_path, nonte_path)
    else:
        nonte = read_single_fasta(nonte_path)
        te_only = read_single_fasta(te_only_path)

    # Generate TE copies if not provided
    te_fasta_path = args.te_fasta if args.te_fasta else os.path.join(args.outdir, 'te_copies.fa')
    if not os.path.exists(te_fasta_path):
        te_copies = generate_te_copies(genome, truth, te_fasta_path)
    else:
        te_copies = read_multifasta_sequences(te_fasta_path)

    ## Load sequences only if needed
    #genome = nonte = te_only = te_copies = None
    #if args.eval in ['kmer','entropy','compress','all']:
    #    if not args.genome or not args.nonte or not args.te_only:
    #        sys.exit("Error: --genome, --nonte, and --te_only are required for kmer/entropy/compress/all")
    #    genome = read_single_fasta(args.genome)
    #    nonte = read_single_fasta(args.nonte)
    #    te_only = read_single_fasta(args.te_only)
    #if args.eval in ['entropy','all']:
    #    if not args.te_fasta:
    #        sys.exit("Error: --te_fasta required for entropy/all")
    #    te_copies = read_multifasta_sequences(args.te_fasta)
    #if args.eval in ['rm','all']:
    #    if not args.te_gff or not args.rm_out:
    #        sys.exit("Error: --te_gff and --rm_out required for rm/all")

    # ---------- k-mer spectra ----------
    if args.eval in ['kmer','all']:
        print('Computing k-mer counts...')
        k = args.k
        bases = 'ACGT'
        raw_kmers = [''.join(p) for p in product(bases, repeat=k)]
        all_kmers = sorted({canonical_kmer(km) for km in raw_kmers})

        counts_genome = count_kmers(genome, k=k, canonical=True)
        counts_nonte = count_kmers(nonte, k=k, canonical=True)
        counts_teonly = count_kmers(te_only, k=k, canonical=True)

        vec_genome = kmer_vector(counts_genome, all_kmers)
        vec_nonte = kmer_vector(counts_nonte, all_kmers)
        vec_teonly = kmer_vector(counts_teonly, all_kmers)

        cos_g_te = cosine_similarity(vec_genome, vec_teonly)
        angle_g_te = spectral_angle_from_cosine(cos_g_te)
        cos_g_non = cosine_similarity(vec_genome, vec_nonte)
        angle_g_non = spectral_angle_from_cosine(cos_g_non)

        print(f'k={k} cosine(genome,TE-only) = {cos_g_te:.6f} spectral_angle={angle_g_te:.6f} rad')
        print(f'k={k} cosine(genome,non-TE) = {cos_g_non:.6f} spectral_angle={angle_g_non:.6f} rad')

        pd.DataFrame({'kmer': all_kmers, 'genome': vec_genome, 'nonTE': vec_nonte, 'TEonly': vec_teonly}).to_csv(os.path.join(args.outdir, f'kmer_counts_k{k}.csv'), index=False)

    # ---------- entropy ----------
    if args.eval in ['entropy','all']:
        print('Computing entropy per TE...')
        te_entropies = []
        for tid, seq in te_copies.items():
            ent = shannon_entropy_bases(seq)
            te_entropies.append({'id': tid, 'len': len(seq), 'entropy': ent})
        pd.DataFrame(te_entropies).to_csv(os.path.join(args.outdir, 'te_entropies.csv'), index=False)

        print('Sliding window entropy on genome, non-TE, TE-only...')
        step = args.step if args.step is not None else args.win
        pos_g, ent_g = sliding_entropy(genome, win=args.win, step=step, k=None)
        pos_non, ent_non = sliding_entropy(nonte, win=args.win, step=step, k=None)
        pos_teonly, ent_teonly = sliding_entropy(te_only, win=args.win, step=step, k=None)

        pd.DataFrame({'pos': pos_g, 'entropy': ent_g}).to_csv(os.path.join(args.outdir, 'sliding_entropy_genome.csv'), index=False)
        pd.DataFrame({'pos': pos_non, 'entropy': ent_non}).to_csv(os.path.join(args.outdir, 'sliding_entropy_nonte.csv'), index=False)
        pd.DataFrame({'pos': pos_teonly, 'entropy': ent_teonly}).to_csv(os.path.join(args.outdir, 'sliding_entropy_teonly.csv'), index=False)

        plt.figure(figsize=(8,4))
        plt.plot(pos_g, ent_g, label='genome')
        plt.plot(pos_non, ent_non, label='non-TE', alpha=0.8)
        plt.plot(pos_teonly, ent_teonly, label='TE-only', alpha=0.8)
        plt.legend()
        plt.xlabel('position (bp)')
        plt.ylabel('Shannon entropy (bits)')
        plt.tight_layout()
        plt.savefig(os.path.join(args.outdir, 'sliding_entropy.png'), dpi=150)
        plt.close()

    # ---------- compressibility ----------
    if args.eval in ['compress','all']:
        print('Estimating compressibility (gzip)...')
        genome_size = len(genome)
        teonly_size = len(te_only)
        nonte_size = len(nonte)
        gz_genome = compressed_size_bytes(genome)
        gz_teonly = compressed_size_bytes(te_only)
        gz_nonte = compressed_size_bytes(nonte)
        comp_df = pd.DataFrame([
            {'region': 'genome', 'raw_bytes': genome_size, 'gz_bytes': gz_genome, 'ratio': gz_genome/genome_size},
            {'region': 'TE-only', 'raw_bytes': teonly_size, 'gz_bytes': gz_teonly, 'ratio': gz_teonly/teonly_size},
            {'region': 'non-TE', 'raw_bytes': nonte_size, 'gz_bytes': gz_nonte, 'ratio': gz_nonte/nonte_size},
        ])
        comp_df.to_csv(os.path.join(args.outdir, 'compression_summary.csv'), index=False)
        print(comp_df)

    # ---------- RepeatMasker recovery ----------
    if args.eval in ['rm','all']:
        print('Parsing RepeatMasker .out...')
        preds = parse_repeatmasker_out(args.rm_out)
        print(f'Parsed {len(preds)} predicted elements from RepeatMasker .out')
        truth = parse_gff_intervals(args.te_gff)
        print(f'Parsed {len(truth)} truth TE intervals from GFF')

        res_df = evaluate_recovery(truth, preds, score_field='score', ro_th=args.ro)
        res_df.to_csv(os.path.join(args.outdir, 'repeatmasker_recovery_by_threshold.csv'), index=False)

        y_true = []
        scores = []
        for p in preds:
            best_ro = 0.0
            for t in truth:
                ro = reciprocal_overlap(p['start'], p['end'], t['start'], t['end'])
                if ro > best_ro:
                    best_ro = ro
            y_true.append(1 if best_ro >= args.ro else 0)
            scores.append(p.get('score', 0.0))

        if len(set(y_true)) > 1:
            precision, recall, thr = precision_recall_curve(y_true, scores)
            pr_auc = auc(recall, precision)
            plt.figure()
            plt.plot(recall, precision, label=f'PR AUC={pr_auc:.3f}')
            plt.xlabel('Recall')
            plt.ylabel('Precision')
            plt.title('RepeatMasker: Precision-Recall')
            plt.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(args.outdir, 'repeatmasker_pr_curve.png'), dpi=150)
            plt.close()
        else:
            print('Not enough class balance in RepeatMasker hits to compute PR curve.')

        summary = res_df.iloc[0]
        print('\nRepeatMasker recovery summary (first threshold):')
        print(summary)

    print('\nDone. Results written to', args.outdir)

if __name__ == '__main__':
    main()










