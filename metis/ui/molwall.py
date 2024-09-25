from PySide6.QtWidgets import QLabel, QPushButton, QScrollArea
from PySide6.QtGui import QPixmap, QPalette, QColor
from PySide6 import QtWidgets, QtCore
from rdkit.Chem import Draw
from rdkit.Chem import AllChem as Chem
from typing import Dict, List
import os, shutil
from metis import PKGDIR


# The `molwall` class is a subclass of `QtWidgets.QMainWindow` that creates a GUI window with a grid
# layout and displays molecular images as bricks in the grid.
class molwall(QtWidgets.QMainWindow):
    def __init__(self, parent=None):
        super(molwall, self).__init__(parent)
        self.setMinimumSize(1300, 700)
        self.setMaximumSize(1300, 700)
        self.mol_dict = dict()

    def initGUI(
        self,
        iteration: int,
        mol_info: list = None,
    ) -> None:

        self.iteration = iteration
        self.mol_dict[iteration] = mol_info
        scrollable = QScrollArea()
        self.globallayout = QtWidgets.QGridLayout()
        self.genBricks()
        self.setBricks()
        self.center = QtWidgets.QWidget()
        self.center.setLayout(self.globallayout)
        scrollable.setWidget(self.center)
        self.setCentralWidget(scrollable)
        self.show()

    def genBricks(self) -> None:
        """
        The function `genBricks` generates a list of `molbrick` objects by iterating over the files in the
        current working directory.
        """

        self.molBricks = []
        for iteration in self.mol_dict:
            for smi in self.mol_dict[iteration]:
                mol_inf = {"Iteration": iteration + 1}  # Add Iteration Number as C
                mol_inf.update(smi["mol_info"])
                self.molBricks.append(
                    molbrick(
                        smi["image"],
                        mol_inf,
                    )
                )

    def setBricks(self) -> None:
        """
        The function sets the layout of bricks in a grid pattern using a counter variable.
        """
        counter = 0
        row = 0
        numMols = sum([len(self.mol_dict[i]) for i in self.mol_dict])
        while counter < numMols:
            self.globallayout.addWidget(
                self.molBricks[counter], row, counter % 5, QtCore.Qt.AlignTop
            )
            counter += 1

            row = counter // 5


# The `molbrick` class is a subclass of `QtWidgets.QWidget` that displays an image specified by
# `image_path` and has a minimum and maximum size of 240x400 pixels.
class molbrick(QtWidgets.QWidget):
    def __init__(self, image: QPixmap, mol_info: List[Dict[str, str]] = None):
        super(molbrick, self).__init__()
        palette = self.palette()
        palette.setColor(QPalette.Window, QColor(255, 255, 255))
        self.setPalette(palette)
        self.setAutoFillBackground(True)

        self.setMinimumSize(240, 400)
        self.setMaximumSize(240, 400)

        mainLayout = QtWidgets.QVBoxLayout()
        pixmap = image.scaled(
            220,
            220,
            QtCore.Qt.KeepAspectRatio,
            mode=QtCore.Qt.SmoothTransformation,
        )

        self.molButton = QLabel()
        self.molButton.setPixmap(pixmap)
        self.molButton.resize(240, 400)
        mainLayout.addWidget(self.molButton)
        if mol_info is not None:
            for key in mol_info:
                mainLayout.addWidget(
                    QLabel(f"{key}{':' if key!= ''else ''} {mol_info[key]}")
                )

        self.setLayout(mainLayout)
