from PySide6 import QtWidgets, QtGui, QtCore


class navigationBar(QtWidgets.QWidget):
    """
    Combines multiple Checkboxes into a single Widget
    """

    directionSignal = QtCore.Signal(int)
    openWindowSignal = QtCore.Signal(str)

    def __init__(self, settings):
        super(navigationBar, self).__init__()
        buttonStyleBlue = (
            "color: #4295f5;background-color: transparent;border: 3px solid #4295f5"
        )
        mainLayout = QtWidgets.QHBoxLayout()
        self.backButton = QtWidgets.QPushButton("< Back")

        self.finishButton = QtWidgets.QPushButton("Finish?")
        self.finishButton.setStyleSheet(buttonStyleBlue)
        self.editButton = QtWidgets.QPushButton("Edit?")
        self.editButton.setStyleSheet(buttonStyleBlue)
        self.sendButton = QtWidgets.QPushButton("Send?")
        self.sendButton.setStyleSheet(buttonStyleBlue)

        self.unratedMolecule = QtWidgets.QPushButton("View")
        self.unratedMolecule.setStyleSheet(buttonStyleBlue)
        self.nextButton = QtWidgets.QPushButton("Next >")

        self.moleculeViewerButton = QtWidgets.QPushButton("Compare")
        self.moleculeViewerButton.setStyleSheet(buttonStyleBlue)

        mainLayout.addWidget(self.backButton)
        mainLayout.addWidget(self.finishButton)
        mainLayout.addWidget(self.editButton)
        mainLayout.addWidget(self.sendButton)
        mainLayout.addWidget(self.unratedMolecule)
        mainLayout.addWidget(self.nextButton)
        if settings["sendButton"]["render"] == False:
            self.sendButton.hide()
        if settings["editButton"]["render"] == False:
            self.editButton.hide()

        self.connectButtons()

        self.setLayout(mainLayout)

    def connectButtons(self):

        self.backButton.clicked.connect(lambda: self.directionButtonClick(-1))
        self.nextButton.clicked.connect(lambda: self.directionButtonClick(+1))

        QtGui.QShortcut(
            QtGui.QKeySequence("left"),
            self.backButton,
            lambda: self.directionButtonClick(-1),
        )

        QtGui.QShortcut(
            QtGui.QKeySequence("right"),
            self.nextButton,
            lambda: self.directionButtonClick(+1),
        )

        self.finishButton.clicked.connect(lambda: self.openWindowClick("finish"))
        self.editButton.clicked.connect(lambda: self.openWindowClick("edit"))
        self.sendButton.clicked.connect(lambda: self.openWindowClick("send"))
        self.unratedMolecule.clicked.connect(lambda: self.openWindowClick("view"))

    def openWindowClick(self, window: str):
        self.openWindowSignal.emit(window)

    def directionButtonClick(self, direction: int):
        self.directionSignal.emit(direction)
