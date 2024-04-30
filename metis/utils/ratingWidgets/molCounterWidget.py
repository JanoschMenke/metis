from PySide2 import QtWidgets, QtCore


class molCounter(QtWidgets.QWidget):
    def __init__(self, totalNumMolecules):
        super(molCounter, self).__init__()
        pageLayout = QtWidgets.QHBoxLayout()
        self.totalNumMol = totalNumMolecules
        self.currentMol = 1
        self.text = QtWidgets.QLabel(
            f"Compound: {self.currentMol} / {self.totalNumMol}"
        )
        self.text.setAlignment(QtCore.Qt.AlignCenter)
        pageLayout.addWidget(self.text)
        self.setLayout(pageLayout)

    def updateCounter(self, newIndex):
        self.text.setText(f"Compound: {newIndex+1} / {self.totalNumMol}")
