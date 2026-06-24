from Bio import AlignIO
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
import numpy as np

# Read the aligned FASTA file
print("Reading alignment...")
alignment = AlignIO.read("brca1_aligned.fasta", "fasta")
print(f"Loaded {len(alignment)} sequences: {[record.id for record in alignment]}")

# Calculate pairwise distances (simple method: count differences)
def calculate_distance_matrix(alignment):
    sequences = [str(record.seq) for record in alignment]
    n = len(sequences)
    distance_matrix = np.zeros((n, n))
    
    for i in range(n):
        for j in range(i+1, n):
            # Count differences
            differences = sum(1 for a, b in zip(sequences[i], sequences[j]) if a != b and a != '-' and b != '-')
            distance = differences / len(sequences[i])
            distance_matrix[i, j] = distance
            distance_matrix[j, i] = distance
    
    return distance_matrix

print("Calculating distances...")
dist_matrix = calculate_distance_matrix(alignment)

# Convert to condensed distance matrix for scipy
condensed_dist = squareform(dist_matrix)

# Build tree using hierarchical clustering
print("Building tree...")
linkage_matrix = linkage(condensed_dist, method='average')

# Get species names
species_names = [record.id for record in alignment]

# Draw and save the tree
print("Drawing tree...")
plt.figure(figsize=(14, 8))
dendrogram(linkage_matrix, labels=species_names, leaf_font_size=12)
plt.title("BRCA1 Phylogenetic Tree (Hierarchical Clustering)", fontsize=14)
plt.xlabel("Species", fontsize=12)
plt.ylabel("Distance", fontsize=12)
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig("brca1_tree.png", dpi=300, bbox_inches='tight')
print("✅ Tree saved as brca1_tree.png!")
plt.show()