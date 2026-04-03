# FASTA Sequence Statistics Calculator

A bioinformatics tool built with Python and Biopython that analyses DNA sequences from FASTA files and generates detailed statistics reports with visualizations.

---

## What This Project Does

This tool takes any `.fasta` file containing DNA sequences and automatically:

- Calculates nucleotide composition (A, T, G, C counts)
- Computes GC content percentage for each sequence
- Detects ATG start codons (potential ORFs)
- Compares sequence lengths
- Exports all results to a CSV report
- Generates 3 publication-ready charts

---

## Sample Output

Running the tool on 4 human gene sequences (HBB, EGFR, TP53, INS) produces:

```
====================================================
       FASTA SEQUENCE STATISTICS REPORT
====================================================

Sequence ID  : NM_000518
Description  : Human HBB (Haemoglobin Beta) gene fragment
Length       : 160 bp
Nucleotides  : A=38  T=37  G=46  C=39  N=0
GC content   : 53.12 %
AT content   : 46.88 %
ATG codons   : 2
```

**Charts generated:**
- `results/nucleotide_chart.png` — nucleotide composition per sequence
- `results/gc_content_chart.png` — GC% comparison with colour coding
- `results/sequence_length_chart.png` — length comparison bar chart
- `results/statistics_report.csv` — full data export

---

## Technologies Used

| Tool | Purpose |
|------|---------|
| Python 3.x | Core programming language |
| Biopython | FASTA parsing, sequence analysis |
| Matplotlib | Data visualization and charts |
| CSV module | Report export |

---

## How to Run

### Step 1 — Install dependencies
```bash
pip install biopython matplotlib
```

### Step 2 — Clone this repository
```bash
git clone https://github.com/YOUR-USERNAME/fasta-sequence-analyzer.git
cd fasta-sequence-analyzer
```

### Step 3 — Run with sample data (no input needed)
```bash
python fasta_analyzer.py
```

### Step 4 — Run with your own FASTA file
```bash
python fasta_analyzer.py your_sequences.fasta
```

Results will be saved in the `results/` folder automatically.

---

## Biological Significance

**GC content** is an important property of DNA sequences:
- High GC content (>60%) → more thermally stable DNA
- Low GC content (<40%) → less stable, easier to denature
- GC content varies significantly between organisms and genes

**ATG codons** mark potential start sites for protein translation (Open Reading Frames), which is fundamental to gene prediction in genomics.

---

## Project Structure

```
fasta-sequence-analyzer/
├── fasta_analyzer.py        ← main analysis script
├── sample_sequences.fasta   ← sample data (4 human genes)
├── results/
│   ├── statistics_report.csv
│   ├── nucleotide_chart.png
│   ├── gc_content_chart.png
│   └── sequence_length_chart.png
└── README.md
```

---

## Author

**Satyam Kumar**
B.Tech Biotechnology | [NIT JALANDHAR]
GitHub: [satyamkumar-22] | LinkedIn: [https://www.linkedin.com/in/satyam-kumar-96a219375]

---

## Future Improvements

- Add protein sequence support
- Implement BLAST search integration
- Add multiple sequence alignment
- Build a simple web interface using Flask
