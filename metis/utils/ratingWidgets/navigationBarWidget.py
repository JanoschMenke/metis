from PySide2 import QtWidgets


class navigationBar(QtWidgets.QWidget):
    """
    Combines multiple Checkboxes into a single Widget
    """

    def __init__(self, settings):
        super(navigationBar, self).__init__()
        buttonStyleBlue = (
            "color: #4295f5;background-color: transparent;border: 3px solid #4295f5"
        )
        mainLayout = QtWidgets.QHBoxLayout()
        self.backButton = QtWidgets.QPushButton("< Back")

        self.downloadButton = QtWidgets.QPushButton("Finish?")
        self.downloadButton.setStyleSheet(buttonStyleBlue)
        self.editButton = QtWidgets.QPushButton("Edit?")
        self.editButton.setStyleSheet(buttonStyleBlue)
        self.sendButton = QtWidgets.QPushButton("Send?")
        self.sendButton.setStyleSheet(buttonStyleBlue)

        self.unratedMolecule = QtWidgets.QPushButton("Unrated Mol.")
        self.unratedMolecule.setStyleSheet(buttonStyleBlue)
        self.downloadButton.pressed.connect(self.downloadClickAnimation)
        self.nextButton = QtWidgets.QPushButton("Next >")

        self.moleculeViewerButton = QtWidgets.QPushButton("Compare")
        self.moleculeViewerButton.setStyleSheet(buttonStyleBlue)

        mainLayout.addWidget(self.backButton)
        mainLayout.addWidget(self.downloadButton)
        mainLayout.addWidget(self.moleculeViewerButton)
        mainLayout.addWidget(self.editButton)
        mainLayout.addWidget(self.sendButton)
        mainLayout.addWidget(self.unratedMolecule)
        mainLayout.addWidget(self.nextButton)
        if settings["sendButton"]["render"] == False:
            self.sendButton.hide()
        if settings["editButton"]["render"] == False:
            self.editButton.hide()
        if settings["compareButton"]["render"] == False:
            self.moleculeViewerButton.hide()

        self.setLayout(mainLayout)

    def downloadClickAnimation(self):
        self.downloadButton.setStyleSheet(
            "color: #ECECEC;background-color: #4295f5;border: 3px solid #4295f5"
        )
        self.downloadButton.setText("Completed")

    def resetDownloadText(self):
        self.downloadButton.setText("Finish?")
        self.downloadButton.setStyleSheet(
            "color: #4295f5;background-color: transparent;border: 3px solid #4295f5"
        )


# %%
