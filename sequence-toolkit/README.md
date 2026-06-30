# Sequence Toolkit

A Python command-line tool for comprehensive DNA sequence analysis. Perform GC content calculations, protein translations, codon usage analysis, and generate publication-quality visualizations.

## Features

- **Sequence Loading:** Parse FASTA files with Biopython
- **GC Content Calculation:** Determine GC percentage for any DNA sequence
- **Reverse Complement:** Generate reverse complement for DNA sequences
- **DNA to Protein Translation:** Translate DNA sequences to proteins in any reading frame
- **Codon Usage Analysis:** Calculate frequency of each codon in your sequences
- **Summary Statistics:** Export sequence metrics to CSV
- **Data Visualization:** Generate bar charts of codon usage frequencies

## Installation

### Requirements
- Python 3.9+
- Biopython
- Pandas
- Matplotlib
- Numpy

### Setup

1. Clone the repository:
```bash
git clone https://github.com/kashihiehie/sequence-toolkit.git
cd sequence-toolkit
```

2. Create a conda environment (recommended):
```bash
conda create -n bioinfo python=3.11
conda activate bioinfo
```

3. Install dependencies:
```bash
pip install biopython pandas matplotlib numpy
```

## Usage

### Basic Usage

```bash
python sequence_toolkit.py <fasta_file>
```

### Example

```bash
python sequence_toolkit.py my_genes.fasta
```

This will generate:
- `sequence_summary.csv` — Summary statistics for each sequence
- `codon_usage.png` — Bar chart of codon frequencies

## Functions

### `load_sequences(filepath)`
Load sequences from a FASTA file.

**Parameters:**
- `filepath` (str): Path to FASTA file

**Returns:**
- list of SeqRecord objects

---

### `calculate_gc_content(sequence)`
Calculate GC content as a percentage.

**Parameters:**
- `sequence` (str or Seq): DNA sequence

**Returns:**
- float: GC percentage (0-100)

**Example:**
```python
gc = calculate_gc_content("ATGCGCAT")
print(f"GC Content: {gc}%")  # Output: GC Content: 50.0%
```

---

### `reverse_complement(sequence)`
Generate reverse complement of a DNA sequence.

**Parameters:**
- `sequence` (str or Seq): DNA sequence

**Returns:**
- str: Reverse complement sequence

**Example:**
```python
rc = reverse_complement("ATGC")
print(rc)  # Output: GCAT
```

---

### `translate_dna_to_protein(sequence, frame=0)`
Translate DNA sequence to protein.

**Parameters:**
- `sequence` (str or Seq): DNA sequence
- `frame` (int): Reading frame (0, 1, or 2)

**Returns:**
- str: Protein sequence (single letter amino acid codes)

**Example:**
```python
protein = translate_dna_to_protein("ATGAAATTT")
print(protein)  # Output: MKF
```

---

### `calculate_codon_usage(sequence)`
Count frequency of each codon.

**Parameters:**
- `sequence` (str or Seq): DNA sequence

**Returns:**
- dict: Codon → frequency mapping

**Example:**
```python
codons = calculate_codon_usage("ATGATGATG")
print(codons)  # Output: {'ATG': 3}
```

---

### `create_summary_table(records)`
Create DataFrame with sequence statistics.

**Parameters:**
- `records` (list): List of SeqRecord objects

**Returns:**
- pandas.DataFrame with columns:
  - Sequence_ID
  - Length_bp
  - GC_Content_%
  - Protein_Length_aa
  - First_50_bp

---

### `save_summary_table(df, output_file="sequence_summary.csv")`
Save summary statistics to CSV file.

**Parameters:**
- `df` (DataFrame): Summary table
- `output_file` (str): Output filename

---

### `plot_codon_usage(codon_counts, output_file="codon_usage.png")`
Create bar chart of codon frequencies.

**Parameters:**
- `codon_counts` (dict): Codon frequency dictionary
- `output_file` (str): Output filename

## Output Files

### sequence_summary.csv
Tab-separated file with:
- Sequence ID
- Length (base pairs)
- GC Content (%)
- Protein length (amino acids)
- First 50 bases

### codon_usage.png
Bar chart showing the top 15 most frequent codons. Useful for analyzing codon bias.

## Workflow Example

```python
from sequence_toolkit import *
from Bio import SeqIO

# Load sequences
records = load_sequences("genes.fasta")

# Analyze first sequence
seq = records[0]
gc = calculate_gc_content(seq.seq)
protein = translate_dna_to_protein(seq.seq)
rc = reverse_complement(seq.seq)
codons = calculate_codon_usage(seq.seq)

print(f"GC Content: {gc}%")
print(f"Protein: {protein}")
print(f"Reverse Complement: {rc}")
print(f"Codons: {codons}")

# Create outputs
df = create_summary_table(records)
save_summary_table(df)
plot_codon_usage(dict(Counter(codons)))
```

## Use Cases

- **Initial Sequence QC:** Quickly assess sequence properties
- **Codon Bias Analysis:** Understand organism-specific codon preferences
- **Batch Processing:** Analyze multiple sequences at once
- **Teaching:** Learn bioinformatics concepts and Python
- **Research:** Prepare data for downstream analysis

## Example Output

**sequence_summary.csv:**
```
Sequence_ID,Length_bp,GC_Content_%,Protein_Length_aa,First_50_bp
gene1,1500,48.5,500,ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGA
gene2,2100,52.3,700,ATGAAACCCGGGTTTATGAAACCCGGGTTTATGAAACCCGGGTTTATGAA
```

**codon_usage.png:**
Bar chart showing frequency of each codon (top 15 displayed)

## Requirements

```
biopython>=1.79
pandas>=1.3.0
matplotlib>=3.4.0
numpy>=1.21.0
```

## Tools Used

- **Biopython:** Sequence parsing and handling
- **Pandas:** Data manipulation and export
- **Matplotlib:** Data visualization
- **Python 3.11:** Core language

## Limitations

- Requires valid FASTA format input
- Translation assumes standard genetic code
- Codon chart shows top 15 codons (configurable in code)
- No support for ambiguous bases (IUPAC codes) in some functions

## Future Enhancements

- Support for multiple reading frames simultaneously
- Codon adaptation index (CAI) calculation
- GC3 (synonymous codon position) analysis
- Protein property predictions (MW, pI, hydrophobicity)
- GUI interface

## Citation

If you use this toolkit in your research, please cite:

```
Sequence Toolkit v1.0
https://github.com/kashihiehie/sequence-toolkit
```

## License

MIT License - feel free to use and modify

## Author

Kashika Vaish
- GitHub: [@kashihiehie](https://github.com/kashihiehie)
- Portfolio: [Bioinformatics Projects](https://github.com/kashihiehie?tab=repositories)

## Support

For issues, questions, or suggestions, please open an issue on GitHub or contact the author.

---

**Happy sequencing!** 🧬