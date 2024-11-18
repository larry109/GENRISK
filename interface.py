import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QVBoxLayout, QPushButton, QLabel, QLineEdit, QTextEdit, QFileDialog
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Interface d'interrogation de la base de connaissances génétiques")
        self.setGeometry(100, 100, 800, 600)

        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout(central_widget)

        self.load_button = QPushButton("Charger la base de connaissances")
        self.load_button.clicked.connect(self.load_knowledge_base)
        layout.addWidget(self.load_button)

        self.query_input = QLineEdit()
        self.query_input.setPlaceholderText("Entrez votre requête ici")
        layout.addWidget(self.query_input)

        self.query_button = QPushButton("Exécuter la requête")
        self.query_button.clicked.connect(self.execute_query)
        layout.addWidget(self.query_button)

        self.result_text = QTextEdit()
        self.result_text.setReadOnly(True)
        layout.addWidget(self.result_text)

        self.figure = plt.figure(figsize=(5, 4))
        self.canvas = FigureCanvas(self.figure)
        layout.addWidget(self.canvas)

        self.df = None

    def load_knowledge_base(self):
        file_name, _ = QFileDialog.getOpenFileName(self, "Sélectionner la base de connaissances", "", "Fichiers CSV (*.csv)")
        if file_name:
            self.df = pd.read_csv(file_name)
            self.result_text.setText(f"Base de connaissances chargée : {len(self.df)} séquences")

    def execute_query(self):
        if self.df is None:
            self.result_text.setText("Veuillez d'abord charger la base de connaissances.")
            return

        query = self.query_input.text().lower()
        
        if "gc content" in query:
            self.show_gc_content_distribution()
        elif "sequence length" in query:
            self.show_sequence_length_distribution()
        elif "cluster" in query:
            self.show_cluster_distribution()
        elif "motif" in query:
            self.show_motif_frequency()
        else:
            self.result_text.setText("Requête non reconnue. Essayez 'GC content', 'sequence length', 'cluster', ou 'motif'.")

    def show_gc_content_distribution(self):
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        self.df['gc_content'].hist(ax=ax, bins=20)
        ax.set_title("Distribution du contenu GC")
        ax.set_xlabel("Contenu GC (%)")
        ax.set_ylabel("Fréquence")
        self.canvas.draw()
        self.result_text.setText(f"Contenu GC moyen : {self.df['gc_content'].mean():.2f}%")

    def show_sequence_length_distribution(self):
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        self.df['longueur'].hist(ax=ax, bins=20)
        ax.set_title("Distribution de la longueur des séquences")
        ax.set_xlabel("Longueur de la séquence")
        ax.set_ylabel("Fréquence")
        self.canvas.draw()
        self.result_text.setText(f"Longueur moyenne des séquences : {self.df['longueur'].mean():.2f}")

    def show_cluster_distribution(self):
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        self.df['cluster'].value_counts().plot(kind='bar', ax=ax)
        ax.set_title("Distribution des clusters")
        ax.set_xlabel("Cluster")
        ax.set_ylabel("Nombre de séquences")
        self.canvas.draw()
        self.result_text.setText(f"Nombre de clusters : {self.df['cluster'].nunique()}")

    def show_motif_frequency(self):
        all_motifs = [motif for motifs in self.df['motifs_repetes'].apply(eval) for motif in motifs]
        motif_freq = pd.Series(all_motifs).value_counts().head(10)
        
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        motif_freq.plot(kind='bar', ax=ax)
        ax.set_title("Top 10 des motifs répétés les plus fréquents")
        ax.set_xlabel("Motif")
        ax.set_ylabel("Fréquence")
        plt.xticks(rotation=45, ha='right')
        self.canvas.draw()
        
        self.result_text.setText(f"Motif le plus fréquent : {motif_freq.index[0]} (fréquence : {motif_freq.iloc[0]})")

def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()