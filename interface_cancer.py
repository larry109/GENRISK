import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, QTextEdit, QFileDialog, QComboBox, QSpinBox, QCheckBox, QScrollArea)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QPixmap
from Bio.SeqUtils import gc_fraction
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

class CancerRiskCalculator(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("GENRISK - Système d'Analyse du Risque de Cancer")
        self.setGeometry(100, 100, 1200, 900)
        self.knowledge_base = None
        self.initUI()

    def initUI(self):
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        main_layout = QVBoxLayout(main_widget)

        # En-tête
        header_widget = QWidget()
        header_widget.setMinimumHeight(70)
        header_widget.setStyleSheet("background-color: white;")
        header_layout = QHBoxLayout(header_widget)

        logo_label = QLabel()
        logo_pixmap = QPixmap("logo_ispm.png")  # Logo vide si le fichier n'existe pas
        logo_pixmap = logo_pixmap.scaled(50, 50, Qt.KeepAspectRatio, Qt.SmoothTransformation)
        logo_label.setPixmap(logo_pixmap)
        header_layout.addWidget(logo_label)

        header_label = QLabel("GENRISK - ONCOLOGIA")
        header_label.setStyleSheet("""
            QLabel {
                font-size: 16px;
                font-weight: bold;
                color: #2c3e50;
                padding: 5px;
            }
        """)
        header_layout.addWidget(header_label)
        header_layout.addStretch()
        main_layout.addWidget(header_widget)

        # Zone de contrôle
        control_group = QWidget()
        control_layout = QVBoxLayout(control_group)

        self.load_kb_button = QPushButton("Charger la base de connaissances")
        self.load_kb_button.setStyleSheet("""
            QPushButton {
                padding: 5px;
                background-color: #4a90e2;
                color: white;
                font-weight: bold;
                border-radius: 5px;
                margin: 5px;
            }
            QPushButton:hover {
                background-color: #357abd;
            }
        """)
        self.load_kb_button.clicked.connect(self.load_knowledge_base)
        control_layout.addWidget(self.load_kb_button)

        sequence_layout = QHBoxLayout()
        sequence_layout.addWidget(QLabel("Séquence ADN:"))
        self.sequence_input = QLineEdit()
        self.sequence_input.setPlaceholderText("Entrez la séquence ADN ici...")
        self.sequence_input.setStyleSheet("""
            QLineEdit {
                padding: 8px;
                border: 1px solid #bdc3c7;
                border-radius: 4px;
            }
        """)
        sequence_layout.addWidget(self.sequence_input)
        control_layout.addLayout(sequence_layout)

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
        control_layout.addLayout(info_layout)

        self.calculate_button = QPushButton("Calculer le risque")
        self.calculate_button.setStyleSheet("""
            QPushButton {
                padding: 10px;
                background-color: #2ecc71;
                color: white;
                font-weight: bold;
                border-radius: 5px;
                margin: 5px;
            }
            QPushButton:hover {
                background-color: #27ae60;
            }
        """)
        self.calculate_button.clicked.connect(self.calculate_risk)
        control_layout.addWidget(self.calculate_button)

        main_layout.addWidget(control_group)

        # Nouveau layout horizontal pour diviser l'écran en deux colonnes
        split_layout = QHBoxLayout()

        # Colonne de gauche pour le résultat
        result_column = QVBoxLayout()
        self.result_text = QTextEdit()
        self.result_text.setReadOnly(True)
        self.result_text.setStyleSheet("""
            QTextEdit {
                background-color: #1e1e1e;
                color: white;
                padding: 5px;
                border: none;
                border-radius: 5px;
                font-family: Consolas, Monaco, monospace;
                margin: 5px;
            }
        """)
        self.result_text.setMinimumHeight(300)
        result_column.addWidget(self.result_text)
        split_layout.addLayout(result_column, 1)  # Ratio 1

        # Colonne de droite pour les graphiques
        graph_column = QVBoxLayout()
        self.graph_scroll = QScrollArea()
        self.graph_scroll.setWidgetResizable(True)
        self.graph_scroll.setStyleSheet("""
            QScrollArea {
                border: none;
                background-color: white;
            }
            QScrollBar:vertical {
                border: none;
                background: #f0f0f0;
                width: 10px;
                margin: 0px;
            }
            QScrollBar::handle:vertical {
                background: #c0c0c0;
                min-height: 30px;
                border-radius: 5px;
            }
            QScrollBar::handle:vertical:hover {
                background: #a0a0a0;
            }
        """)

        self.graph_widget = QWidget()
        self.graph_layout = QVBoxLayout(self.graph_widget)
        self.canvas = FigureCanvas(Figure(figsize=(10, 12)))
        self.graph_layout.addWidget(self.canvas)
        self.graph_widget.setMinimumSize(800, 1000)
        self.graph_scroll.setWidget(self.graph_widget)
        graph_column.addWidget(self.graph_scroll)
        split_layout.addLayout(graph_column, 2)  # Ratio 2 (plus grand que la colonne de résultat)

        main_layout.addLayout(split_layout)

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

        # Vérification de la séquence ADN
        valid_nucleotides = set('ATCG')
        invalid_chars = set(sequence) - valid_nucleotides
        if invalid_chars:
            self.result_text.setText(f"Erreur : La séquence contient des caractères invalides : {', '.join(invalid_chars)}. "
                                     f"Seuls les nucléotides A, T, C et G sont autorisés.")
            return

        gc_content = gc_fraction(sequence) * 100

        closest_match = self.knowledge_base.iloc[(self.knowledge_base['gc_content'] - gc_content).abs().argsort()[:1]]
        risk_score = closest_match['risk_score'].values[0]
        risk_category = closest_match['risk_category'].values[0]

        age = self.age_spin.value()
        is_smoker = self.smoker_check.isChecked()
        sex = self.sex_combo.currentText()

        if age > 50:
            risk_score *= 1.2
        if is_smoker:
            risk_score *= 1.5
        if sex == "Homme":
            risk_score *= 1.1

        self.plot_graphs(gc_content, risk_score)

        result = f"Analyse de la séquence ADN :\n\n"
        result += f"Longueur de la séquence : {len(sequence)} nucléotides\n"
        result += f"Contenu GC : {gc_content:.2f}%\n\n"
        result += f"Score de risque estimé : {risk_score:.4f}\n"
        result += f"Catégorie de risque : {risk_category}\n\n"
        result += f"Facteurs de risque :\n"
        result += f"• Sexe : {sex}\n"
        result += f"• Âge : {age} ans\n"
        result += f"• Fumeur : {'Oui' if is_smoker else 'Non'}\n\n"
        result += "Note : Cette évaluation est basée sur une analyse simplifiée. "
        result += "Pour une évaluation précise des risques de cancer, consultez un professionnel de santé."

        self.result_text.setText(result)

    def plot_graphs(self, current_gc=None, current_risk=None):
        if self.knowledge_base is None:
            return

        self.canvas.figure.clear()

        gs = self.canvas.figure.add_gridspec(2, 1, height_ratios=[1, 1], hspace=0.4)
        ax1 = self.canvas.figure.add_subplot(gs[0])
        ax2 = self.canvas.figure.add_subplot(gs[1])

        plt.rcParams.update({'font.size': 10})

        data = self.knowledge_base.sort_values('gc_content')
        ax1.scatter(data['gc_content'], data['risk_score'], alpha=0.3, color='blue', label='Données')
        z = np.polyfit(data['gc_content'], data['risk_score'], 3)
        p = np.poly1d(z)
        x_new = np.linspace(data['gc_content'].min(), data['gc_content'].max(), 100)
        ax1.plot(x_new, p(x_new), 'r-', label='Tendance', linewidth=2)

        if current_gc is not None and current_risk is not None:
            ax1.scatter([current_gc], [current_risk], color='green', s=100, marker='*', label='Séquence actuelle')

        ax1.set_title('Relation entre Contenu GC et Score de Risque', pad=20, fontsize=12, fontweight='bold')
        ax1.set_xlabel('Contenu GC (%)', labelpad=10, fontsize=10)
        ax1.set_ylabel('Score de Risque', labelpad=10, fontsize=10)
        ax1.grid(True, alpha=0.3)
        ax1.legend(fontsize=8)

        ax2.hist(self.knowledge_base['gc_content'], bins=30, color='skyblue', edgecolor='black', alpha=0.7)
        if current_gc is not None:
            ax2.axvline(x=current_gc, color='green', linestyle='--', label='Séquence actuelle')

        ax2.set_title('Distribution du Contenu GC', pad=20, fontsize=12, fontweight='bold')
        ax2.set_xlabel('Contenu GC (%)', labelpad=10, fontsize=10)
        ax2.set_ylabel('Fréquence', labelpad=10, fontsize=10)
        ax2.grid(True, alpha=0.3)

        mean_gc = self.knowledge_base['gc_content'].mean()
        std_gc = self.knowledge_base['gc_content'].std()
        stats_text = f'Moyenne: {mean_gc:.2f}%\nÉcart-type: {std_gc:.2f}%'
        ax2.text(0.95, 0.95, stats_text, transform=ax2.transAxes, verticalalignment='top', horizontalalignment='right',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8), fontsize=8)

        self.canvas.figure.tight_layout()
        self.canvas.draw()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self.plot_graphs()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    calculator = CancerRiskCalculator()
    calculator.show()
    sys.exit(app.exec_())
