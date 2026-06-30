from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter

def load_sequences(filepath):
    """Load sequences from FASTA file."""
    records = list(SeqIO.parse(filepath, "fasta"))
    if not records:
        print("Error: No sequences found in file")
        return None
    return records

def print_sequence_info(record):
    """Print basic info about a sequence."""
    print(f"\nSequence ID: {record.id}")
    print(f"Length: {len(record.seq)} bp")
    print(f"Sequence (first 50 bp): {str(record.seq)[:50]}...")

def calculate_gc_content(sequence):
    """Calculate GC content as percentage."""
    seq_str = str(sequence).upper()
    if len(seq_str) == 0:
        return 0
    
    gc_count = seq_str.count('G') + seq_str.count('C')
    gc_percent = (gc_count / len(seq_str)) * 100
    return round(gc_percent, 2)

def reverse_complement(sequence):
    """Generate reverse complement of DNA sequence."""
    complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 
                      'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
    
    # First reverse the sequence
    reversed_seq = str(sequence)[::-1]
    
    # Then complement each base
    rc = ''.join([complement_map.get(base, 'N') for base in reversed_seq])
    return rc

def translate_dna_to_protein(sequence, frame=0):
    """Translate DNA sequence to protein (single reading frame)."""
    codon_table = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
    
    seq_str = str(sequence).upper()
    protein = ""
    
    # Translate from specified reading frame
    for i in range(frame, len(seq_str) - 2, 3):
        codon = seq_str[i:i+3]
        amino_acid = codon_table.get(codon, 'X')
        protein += amino_acid
    
    return protein

def calculate_codon_usage(sequence):
    """Count frequency of each codon in sequence."""
    seq_str = str(sequence).upper()
    codon_counts = {}
    
    for i in range(0, len(seq_str) - 2, 3):
        codon = seq_str[i:i+3]
        if len(codon) == 3 and 'N' not in codon:  # Only count complete codons
            codon_counts[codon] = codon_counts.get(codon, 0) + 1
    
    return codon_counts

def create_summary_table(records):
    """Create pandas DataFrame with sequence statistics."""
    data = []
    
    for record in records:
        gc = calculate_gc_content(record.seq)
        protein = translate_dna_to_protein(record.seq, frame=0)
        
        data.append({
            'Sequence_ID': record.id,
            'Length_bp': len(record.seq),
            'GC_Content_%': gc,
            'Protein_Length_aa': len(protein),
            'First_50_bp': str(record.seq)[:50]
        })
    
    df = pd.DataFrame(data)
    return df

def save_summary_table(df, output_file="sequence_summary.csv"):
    """Save summary table to CSV."""
    df.to_csv(output_file, index=False)
    print(f"\n✅ Summary saved to {output_file}")
    print(df)

def plot_codon_usage(codon_counts, output_file="codon_usage.png"):
    """Create bar chart of codon frequencies."""
    if not codon_counts:
        print("No codons to plot")
        return
    
    # Get top 15 codons
    top_codons = dict(sorted(codon_counts.items(), 
                             key=lambda x: x[1], 
                             reverse=True)[:15])
    
    plt.figure(figsize=(12, 6))
    plt.bar(top_codons.keys(), top_codons.values(), color='steelblue')
    plt.xlabel('Codon')
    plt.ylabel('Frequency')
    plt.title('Top 15 Codon Usage')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    print(f"✅ Chart saved to {output_file}")
    plt.close()

if __name__ == "__main__":
    print("🔍 DEBUG: Starting analysis...")
    
    # Create test FASTA
    test_fasta = "test_seq.fasta"
    print(f"📝 Creating {test_fasta}...")
    with open(test_fasta, 'w') as f:
        f.write(">gene1\n")
        f.write("ATGATGATGATGATGATGATGATGATGATGATGATGATGATG\n")
        f.write(">gene2\n")
        f.write("ATGAAACCCGGGTTTATGAAACCCGGGTTTATGAAACCCGGGTTTATG\n")
    print(f"✅ {test_fasta} created")
    
    # Load sequences
    print("\n📂 Loading sequences...")
    records = load_sequences(test_fasta)
    print(f"✅ Loaded {len(records) if records else 0} sequences")
    
    if records:
        print("\n📊 Creating summary table...")
        df = create_summary_table(records)
        print("✅ Summary table created:")
        print(df)
        
        print("\n💾 Saving summary table...")
        save_summary_table(df, "sequence_summary.csv")
        print("✅ Summary saved")
        
        print("\n📈 Creating plot...")
        codons = calculate_codon_usage(records[0].seq)
        plot_codon_usage(codons, "codon_usage.png")
        print("✅ Plot created")
        
        print("\n✅ ALL DONE!")
    else:
        print("❌ No sequences loaded!")