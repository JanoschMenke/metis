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
        self.cwd = f"{PKGDIR}/resources/temp_images/molImages/"
        if os.path.exists(self.cwd):
            shutil.rmtree(self.cwd)
        os.makedirs(self.cwd)
        self.image_paths = []
        self.molInfoDicts = []
        self.smilesDict = {0: {}}

    def initGUI(
        self,
        iteration: int,
        smiles: List[str] = None,
        molInfo: List[Dict[str, str]] = None,
    ) -> None:
        self.iteration = iteration
        if iteration not in self.smilesDict:
            self.smilesDict[iteration] = {}
        scrollable = QScrollArea()
        self.globallayout = QtWidgets.QGridLayout()
        if smiles is not None:
            self.loadMolecules(smiles=smiles, infos=molInfo)
        self.genBricks()
        self.setBricks()
        self.center = QtWidgets.QWidget()
        self.center.setLayout(self.globallayout)
        scrollable.setWidget(self.center)
        self.setCentralWidget(scrollable)
        self.show()

    def loadMolecules(
        self, smiles: List[str], infos: List[Dict[str, str]] = None
    ) -> None:
        """
        The function `loadMolecules` takes a list of SMILES strings and an optional list of dictionaries as
        input, and saves images of the molecules represented by the SMILES strings as PNG files.
        """

        smiles_to_remove = set(self.smilesDict[self.iteration]) - set(smiles)
        [self.smilesDict[self.iteration].pop(x) for x in smiles_to_remove]

        for i, smi in enumerate(smiles):
            self.smilesDict[self.iteration][smi] = {"molInfo": infos[i]}
            img = Draw.MolToImage(Chem.MolFromSmiles(smi), size=(240, 240), canvas=None)
            img.save(f"{self.cwd}mol_{smi}.png")
            self.smilesDict[self.iteration][smi]["imgPath"] = f"{self.cwd}mol_{smi}.png"

    def genBricks(self) -> None:
        """
        The function `genBricks` generates a list of `molbrick` objects by iterating over the files in the
        current working directory.
        """

        self.molBricks = []
        for iteration in self.smilesDict:
            for smi in self.smilesDict[iteration]:
                self.molBricks.append(
                    molbrick(
                        self.smilesDict[iteration][smi]["imgPath"],
                        self.smilesDict[iteration][smi]["molInfo"],
                    )
                )

    def clearImages(self) -> None:
        """
        The `clearImages` function deletes all images in the directory.
        """
        shutil.rmtree(self.cwd)

    def setBricks(self) -> None:
        """
        The function sets the layout of bricks in a grid pattern using a counter variable.
        """
        counter = 0
        row = 0
        numMols = sum([len(self.smilesDict[i]) for i in self.smilesDict])
        while counter < numMols:
            self.globallayout.addWidget(
                self.molBricks[counter], row, counter % 5, QtCore.Qt.AlignTop
            )
            counter += 1

            row = counter // 5


# The `molbrick` class is a subclass of `QtWidgets.QWidget` that displays an image specified by
# `image_path` and has a minimum and maximum size of 240x400 pixels.
class molbrick(QtWidgets.QWidget):
    def __init__(self, image_path: str, molInfo: List[Dict[str, str]] = None):
        super(molbrick, self).__init__()
        palette = self.palette()
        palette.setColor(QPalette.Window, QColor(255, 255, 255))
        self.setPalette(palette)
        self.setAutoFillBackground(True)

        self.setMinimumSize(240, 400)
        self.setMaximumSize(240, 400)

        mainLayout = QtWidgets.QVBoxLayout()
        pixmap = QPixmap(image_path)
        pixmap = pixmap.scaled(
            220,
            220,
            QtCore.Qt.KeepAspectRatio,
            mode=QtCore.Qt.SmoothTransformation,
        )

        self.molButton = QLabel()
        self.molButton.setPixmap(pixmap)
        self.molButton.resize(240, 400)
        mainLayout.addWidget(self.molButton)
        if molInfo is not None:
            for key in molInfo:
                mainLayout.addWidget(
                    QLabel(f"{key}{':' if key!= ''else ''} {molInfo[key]}")
                )

        self.setLayout(mainLayout)
