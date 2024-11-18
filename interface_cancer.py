import sys
import pandas as pd
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, 
                             QLabel, QLineEdit, QPushButton, QTextEdit, QFileDialog, 
                             QComboBox, QSpinBox, QCheckBox)
from PyQt5.QtCore import Qt
from Bio.SeqUtils import gc_fraction
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np
class CancerRiskCalculator(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Calculateur de Risque de Cancer basé sur l'ADN")
        self.setGeometry(100, 100, 800, 600)
        self.knowledge_base = None
        self.initUI()

    def initUI(self):
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)

        # Bouton pour charger la base de connaissances
        self.load_kb_button = QPushButton("Charger la base de connaissances")
        self.load_kb_button.clicked.connect(self.load_knowledge_base)
        main_layout.addWidget(self.load_kb_button)

        # Champ pour entrer la séquence ADN
        sequence_layout = QHBoxLayout()
        sequence_layout.addWidget(QLabel("Séquence ADN:"))
        self.sequence_input = QLineEdit()
        sequence_layout.addWidget(self.sequence_input)
        main_layout.addLayout(sequence_layout)

        # Informations supplémentaires
        info_layout = QHBoxLayout()
        info_layout.addWidget(QLabel("Sexe:"))
        self.sex_combo = QComboBox()
        self.sex_combo.addItems(["Homme", "Femme"])
        info_layout.addWidget(self.sex_combo)

        info_layout.addWidget(QLabel("Âge:"))
        self.age_spin = QSpinBox()
        self.age_spin.setRange(0, 120)
        info_layout.addWidget(self.age_spin)

        self.smoker_check = QCheckBox("Fumeur")
        info_layout.addWidget(self.smoker_check)

        main_layout.addLayout(info_layout)

        # Bouton de calcul
        self.calculate_button = QPushButton("Calculer le risque")
        self.calculate_button.clicked.connect(self.calculate_risk)
        main_layout.addWidget(self.calculate_button)

        # Zone de résultat
        self.result_text = QTextEdit()
        self.result_text.setReadOnly(True)
        main_layout.addWidget(self.result_text)

        # Zone pour afficher le graphique
        self.canvas = FigureCanvas(Figure(figsize=(5, 5)))
        main_layout.addWidget(self.canvas)

    def load_knowledge_base(self):
        file_name, _ = QFileDialog.getOpenFileName(self, "Charger la base de connaissances", "", "CSV Files (*.csv)")
        if file_name:
            try:
                self.knowledge_base = pd.read_csv(file_name)
                required_columns = ['gc_content', 'risk_score', 'risk_category']
                if not all(col in self.knowledge_base.columns for col in required_columns):
                    raise ValueError(f"Les colonnes {', '.join(required_columns)} sont requises dans le fichier CSV.")
                self.result_text.setText("Base de connaissances chargée avec succès.")
                self.plot_graphs()
            except Exception as e:
                self.result_text.setText(f"Erreur lors du chargement de la base de connaissances: {str(e)}")

    def calculate_risk(self):
        if self.knowledge_base is None:
            self.result_text.setText("Veuillez d'abord charger la base de connaissances.")
            return

        sequence = self.sequence_input.text().upper()
        if not sequence:
            self.result_text.setText("Veuillez entrer une séquence ADN.")
            return

        # Calcul du contenu GC
        gc_content = gc_fraction(sequence) * 100

        # Recherche dans la base de connaissances
        closest_match = self.knowledge_base.iloc[(self.knowledge_base['gc_content'] - gc_content).abs().argsort()[:1]]
        risk_score = closest_match['risk_score'].values[0]
        risk_category = closest_match['risk_category'].values[0]

        # Ajustement du risque en fonction des informations supplémentaires
        age = self.age_spin.value()
        is_smoker = self.smoker_check.isChecked()
        sex = self.sex_combo.currentText()

        if age > 50:
            risk_score *= 1.2
        if is_smoker:
            risk_score *= 1.5
        if sex == "Homme":
            risk_score *= 1.1

        result = f"Analyse de la séquence ADN :\n\n"
        result += f"Longueur de la séquence : {len(sequence)} nucléotides\n"
        result += f"Contenu GC : {gc_content:.2f}%\n\n"
        result += f"Score de risque estimé : {risk_score:.4f}\n"
        result += f"Catégorie de risque : {risk_category}\n\n"
        result += f"Informations supplémentaires :\n"
        result += f"Sexe : {sex}\n"
        result += f"Âge : {age} ans\n"
        result += f"Fumeur : {'Oui' if is_smoker else 'Non'}\n\n"
        result += "Note : Cette évaluation est basée sur une analyse simplifiée. "
        result += "Pour une évaluation précise des risques de cancer, consultez un professionnel de santé."

        self.result_text.setText(result)

    def plot_graphs(self):
        """Affiche des graphiques améliorés et plus lisibles."""
        self.canvas.figure.clear()
        
        # Création de deux sous-graphiques
        gs = self.canvas.figure.add_gridspec(2, 1, height_ratios=[1, 1], hspace=0.3)
        ax1 = self.canvas.figure.add_subplot(gs[0])
        ax2 = self.canvas.figure.add_subplot(gs[1])
        
        # 1. Graphe en courbe lissé : Contenu GC vs Score de Risque
        # Trier les données pour une meilleure visualisation
        data = self.knowledge_base.sort_values('gc_content')
        
        # Créer une courbe lissée
        ax1.plot(data['gc_content'], data['risk_score'], 
                'b-', alpha=0.3, label='Données brutes')
        
        # Ajouter une ligne de tendance
        z = np.polyfit(data['gc_content'], data['risk_score'], 3)
        p = np.poly1d(z)
        ax1.plot(data['gc_content'], p(data['gc_content']), 
                'r-', linewidth=2, label='Tendance')
        
        ax1.set_title('Relation entre Contenu GC et Score de Risque', 
                    fontsize=10, pad=10)
        ax1.set_xlabel('Contenu GC (%)')
        ax1.set_ylabel('Score de Risque')
        ax1.grid(True, alpha=0.3)
        ax1.legend()

        # 2. Histogramme amélioré
        ax2.hist(self.knowledge_base['gc_content'], 
                bins=30, 
                color='skyblue',
                edgecolor='black',
                alpha=0.7)
        
        ax2.set_title('Distribution du Contenu GC', 
                    fontsize=10, pad=10)
        ax2.set_xlabel('Contenu GC (%)')
        ax2.set_ylabel('Fréquence')
        ax2.grid(True, alpha=0.3)

        # Ajouter des statistiques descriptives
        mean_gc = self.knowledge_base['gc_content'].mean()
        std_gc = self.knowledge_base['gc_content'].std()
        
        stats_text = f'Moyenne: {mean_gc:.2f}%\nÉcart-type: {std_gc:.2f}%'
        ax2.text(0.95, 0.95, stats_text,
                transform=ax2.transAxes,
                verticalalignment='top',
                horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        # Ajuster la mise en page
        self.canvas.figure.set_tight_layout(True)
        self.canvas.draw()


if __name__ == '__main__':
    app = QApplication(sys.argv)
    calculator = CancerRiskCalculator()
    calculator.show()
    sys.exit(app.exec_())
