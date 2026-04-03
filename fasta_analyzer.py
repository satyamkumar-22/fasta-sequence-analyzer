# ============================================================
#  FASTA Sequence Statistics Calculator
#  Author  : Satyam Kumar
#  Purpose : Analyse DNA/protein sequences from a FASTA file
#            and export a statistics report as CSV + plots
# ============================================================

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import csv
import os
import sys


# ── 1. PARSE FASTA FILE ─────────────────────────────────────

def parse_fasta(filepath):
    """Read all sequences from a FASTA file."""
    records = list(SeqIO.parse(filepath, "fasta"))
    if not records:
        print("ERROR: No sequences found in file.")
        sys.exit(1)
    print(f"  Loaded {len(records)} sequence(s) from {filepath}\n")
    return records


# ── 2. CALCULATE STATISTICS ──────────────────────────────────

def calculate_stats(record):
    """Return a dict of statistics for one sequence."""
    seq   = str(record.seq).upper()
    length = len(seq)

    a = seq.count("A")
    t = seq.count("T")
    g = seq.count("G")
    c = seq.count("C")
    n = seq.count("N")          # unknown bases

    gc = round((g + c) / length * 100, 2) if length > 0 else 0
    at = round((a + t) / length * 100, 2) if length > 0 else 0

    # Simple ORF finder: count ATG start codons
    orfs = seq.count("ATG")

    return {
        "ID"          : record.id,
        "Description" : record.description,
        "Length"      : length,
        "A_count"     : a,
        "T_count"     : t,
        "G_count"     : g,
        "C_count"     : c,
        "N_count"     : n,
        "GC_percent"  : gc,
        "AT_percent"  : at,
        "ATG_codons"  : orfs,
    }


# ── 3. PRINT REPORT TO TERMINAL ──────────────────────────────

def print_report(stats_list):
    """Pretty-print all statistics."""
    print("=" * 60)
    print("       FASTA SEQUENCE STATISTICS REPORT")
    print("=" * 60)
    for s in stats_list:
        print(f"\nSequence ID  : {s['ID']}")
        print(f"Description  : {s['Description'][:60]}")
        print(f"Length       : {s['Length']} bp")
        print(f"Nucleotides  : A={s['A_count']}  T={s['T_count']}  "
              f"G={s['G_count']}  C={s['C_count']}  N={s['N_count']}")
        print(f"GC content   : {s['GC_percent']} %")
        print(f"AT content   : {s['AT_percent']} %")
        print(f"ATG codons   : {s['ATG_codons']}")
    print("\n" + "=" * 60)


# ── 4. EXPORT CSV REPORT ─────────────────────────────────────

def export_csv(stats_list, output_path="results/statistics_report.csv"):
    """Save statistics to a CSV file."""
    os.makedirs("results", exist_ok=True)
    fields = list(stats_list[0].keys())
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(stats_list)
    print(f"  CSV report saved  →  {output_path}")


# ── 5. PLOT 1 — NUCLEOTIDE BAR CHART ─────────────────────────

def plot_nucleotide_chart(stats_list):
    """Bar chart of nucleotide counts for each sequence."""
    os.makedirs("results", exist_ok=True)

    seq_ids = [s["ID"][:15] for s in stats_list]
    a_vals  = [s["A_count"] for s in stats_list]
    t_vals  = [s["T_count"] for s in stats_list]
    g_vals  = [s["G_count"] for s in stats_list]
    c_vals  = [s["C_count"] for s in stats_list]

    x      = range(len(seq_ids))
    width  = 0.2
    colors = ["#4CAF50", "#2196F3", "#FF9800", "#E91E63"]

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar([i - 1.5*width for i in x], a_vals, width, label="A", color=colors[0])
    ax.bar([i - 0.5*width for i in x], t_vals, width, label="T", color=colors[1])
    ax.bar([i + 0.5*width for i in x], g_vals, width, label="G", color=colors[2])
    ax.bar([i + 1.5*width for i in x], c_vals, width, label="C", color=colors[3])

    ax.set_xticks(list(x))
    ax.set_xticklabels(seq_ids, rotation=30, ha="right")
    ax.set_ylabel("Nucleotide count")
    ax.set_title("Nucleotide Composition per Sequence")
    ax.legend()
    plt.tight_layout()
    path = "results/nucleotide_chart.png"
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"  Nucleotide chart  →  {path}")


# ── 6. PLOT 2 — GC CONTENT BAR CHART ────────────────────────

def plot_gc_chart(stats_list):
    """Horizontal bar chart showing GC% for each sequence."""
    os.makedirs("results", exist_ok=True)

    seq_ids = [s["ID"][:20] for s in stats_list]
    gc_vals = [s["GC_percent"] for s in stats_list]

    colors = ["#E74C3C" if g < 40 else "#27AE60" if g > 60 else "#3498DB"
              for g in gc_vals]

    fig, ax = plt.subplots(figsize=(9, max(3, len(seq_ids) * 0.6)))
    bars = ax.barh(seq_ids, gc_vals, color=colors, edgecolor="white")

    ax.set_xlabel("GC Content (%)")
    ax.set_title("GC Content per Sequence")
    ax.axvline(50, color="gray", linestyle="--", linewidth=1, label="50% line")

    for bar, val in zip(bars, gc_vals):
        ax.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height()/2,
                f"{val}%", va="center", fontsize=9)

    low   = mpatches.Patch(color="#E74C3C", label="Low GC  (<40%)")
    mid   = mpatches.Patch(color="#3498DB", label="Medium GC (40–60%)")
    high  = mpatches.Patch(color="#27AE60", label="High GC  (>60%)")
    ax.legend(handles=[low, mid, high], loc="lower right", fontsize=8)

    plt.tight_layout()
    path = "results/gc_content_chart.png"
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"  GC content chart  →  {path}")


# ── 7. PLOT 3 — SEQUENCE LENGTH CHART ───────────────────────

def plot_length_chart(stats_list):
    """Bar chart of sequence lengths."""
    os.makedirs("results", exist_ok=True)

    seq_ids = [s["ID"][:15] for s in stats_list]
    lengths = [s["Length"] for s in stats_list]

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.bar(seq_ids, lengths, color="#9B59B6", edgecolor="white")
    ax.set_ylabel("Length (bp)")
    ax.set_title("Sequence Length Comparison")
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()
    path = "results/sequence_length_chart.png"
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"  Length chart      →  {path}")


# ── 8. CREATE SAMPLE FASTA (for testing) ─────────────────────

def create_sample_fasta():
    """Write a sample FASTA file so the user can test immediately."""
    sample = """>NM_000518 Human HBB (Haemoglobin Beta) gene fragment
ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGT
GAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCT
>NM_005228 Human EGFR (Epidermal Growth Factor Receptor) fragment
ATGCGACCCTCCGGGACGGCCGGGGCAGCGCTCCTGGCGCTGCTGGCTGCGCTCTGCCCGGCGAGTCGGGCTCTGGAG
GAAAAGAAAGTTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGA
>NM_000546 Human TP53 (Tumour Suppressor) gene fragment
ATGGAGGAGCCGCAGTCAGATCCTAGCGTTGAGCCACCCGAGCCCAGGAGCCCAGAGCAGCGAGATCGAGAGAGGAGG
AGGATGAGCACTGAAGCGAAAATGGTTTTTCAGACCTATGGAAACTACTTCCTGAAAACAACGTTCTGTCCCCCTTGCC
>NM_000207 Human INS (Insulin) gene fragment
ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCCGCAGCCTTTGTG
AACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTGCGGGGAACGAGGCTTCTTCTACACACCCG
"""
    with open("sample_sequences.fasta", "w") as f:
        f.write(sample)
    print("  Sample FASTA created  →  sample_sequences.fasta\n")
    return "sample_sequences.fasta"


# ── MAIN ─────────────────────────────────────────────────────

def main():
    print("\n" + "=" * 60)
    print("   FASTA SEQUENCE STATISTICS CALCULATOR")
    print("   Bioinformatics Project | Python + Biopython")
    print("=" * 60 + "\n")

    # Use provided file or create sample
    if len(sys.argv) > 1:
        fasta_file = sys.argv[1]
        if not os.path.exists(fasta_file):
            print(f"ERROR: File '{fasta_file}' not found.")
            sys.exit(1)
    else:
        print("No input file given. Creating sample FASTA for demo...\n")
        fasta_file = create_sample_fasta()

    # Run analysis
    records    = parse_fasta(fasta_file)
    stats_list = [calculate_stats(r) for r in records]

    print_report(stats_list)

    print("\nSaving outputs to results/ folder...")
    export_csv(stats_list)
    plot_nucleotide_chart(stats_list)
    plot_gc_chart(stats_list)
    plot_length_chart(stats_list)

    print("\nDone! Open the results/ folder to see your charts and CSV.\n")


if __name__ == "__main__":
    main()
