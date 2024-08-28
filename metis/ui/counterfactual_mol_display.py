from PySide6.QtWidgets import (
    QLabel,
)
from PySide6.QtGui import QPixmap
from PySide6 import QtWidgets, QtCore
from rdkit.Chem import Draw
from rdkit.Chem import AllChem as Chem


class CounterfactualWindow(QtWidgets.QMainWindow):
    def __init__(self, parent=None):
        super(CounterfactualWindow, self).__init__(parent)
        self.setMinimumSize(1300, 700)
        self.setMaximumSize(1300, 700)

    def initGUI(self, smiles):
        self.globallayout = QtWidgets.QHBoxLayout()
        self.smiles = smiles
        self.gen_mol_images()
        for i in range(len(self.smiles)):
            image = QLabel(self)
            image.setPixmap(QPixmap(f"../data/molimage_{i}.png"))
            self.globallayout.addWidget(image)
        self.center = QtWidgets.QWidget()
        self.center.setLayout(self.globallayout)
        self.setCentralWidget(self.center)

        self.show()

    def gen_mol_images(self):
        for i, x in enumerate(self.smiles):
            Draw.MolToFile(
                Chem.MolFromSmiles(x), f"../data/molimage_{i}.png", size=(300, 300)
            )
