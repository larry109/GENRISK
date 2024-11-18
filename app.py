# Importation des bibliothèques nécessaires
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
import numpy as np

# Fonction pour lire les séquences d'un fichier FASTA
def lire_sequences(fichier):
    sequences = []
    for record in SeqIO.parse(fichier, "fasta"):
        sequences.append({"id": record.id, "sequence": str(record.seq)})
    return pd.DataFrame(sequences)

# Fonction pour calculer le contenu GC
def calculer_gc(sequence):
    return GC(sequence)

# Fonction pour trouver les motifs répétés
def trouver_motifs_repetes(sequence, longueur_motif):
    motifs = {}
    for i in range(len(sequence) - longueur_motif + 1):
        motif = sequence[i:i+longueur_motif]
        if motif in motifs:
            motifs[motif] += 1
        else:
            motifs[motif] = 1
    return {k: v for k, v in motifs.items() if v > 1}

# Fonction principale
def analyser_sequences(fichier_fasta):
    # Lecture des séquences
    df = lire_sequences(fichier_fasta)
    print(f"Nombre de séquences lues : {len(df)}")

    # Calcul du contenu GC
    df['gc_content'] = df['sequence'].apply(calculer_gc)

    # Recherche de motifs répétés
    df['motifs_repetes'] = df['sequence'].apply(lambda x: trouver_motifs_repetes(x, 3))

    # Calcul de la longueur des séquences
    df['longueur'] = df['sequence'].apply(len)

    # Visualisation de la distribution du contenu GC
    plt.figure(figsize=(10, 6))
    sns.histplot(df['gc_content'], kde=True)
    plt.title("Distribution du contenu GC")
    plt.xlabel("Contenu GC (%)")
    plt.ylabel("Fréquence")
    plt.savefig("distribution_gc.png")
    plt.close()

    # Clustering basé sur le contenu GC et la longueur des séquences
    X = df[['gc_content', 'longueur']].values
    kmeans = KMeans(n_clusters=3, random_state=42)
    df['cluster'] = kmeans.fit_predict(X)

    # Visualisation des clusters
    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=df, x='gc_content', y='longueur', hue='cluster', palette='deep')
    plt.title("Clustering des séquences")
    plt.xlabel("Contenu GC (%)")
    plt.ylabel("Longueur de la séquence")
    plt.savefig("clustering_sequences.png")
    plt.close()

    # Analyse des motifs répétés
    tous_motifs = [motif for motifs in df['motifs_repetes'] for motif in motifs]
    motifs_frequents = pd.Series(tous_motifs).value_counts().head(10)

    plt.figure(figsize=(12, 6))
    motifs_frequents.plot(kind='bar')
    plt.title("Top 10 des motifs répétés les plus fréquents")
    plt.xlabel("Motif")
    plt.ylabel("Fréquence")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig("motifs_repetes_frequents.png")
    plt.close()

    # Sauvegarde des résultats
    df.to_csv("resultats_analyse.csv", index=False)
    print("Analyse terminée. Les résultats ont été sauvegardés.")

# Exécution de l'analyse
if __name__ == "__main__":
    fichier_fasta = "sequences.fasta"  # Remplacez par le chemin de votre fichier FASTA
    analyser_sequences(fichier_fasta)
