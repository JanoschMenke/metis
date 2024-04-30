from PySide2 import QtWidgets, QtCore


class selectLocalSubstructures(QtWidgets.QWidget):
    def __init__(self, settings):
        super(selectLocalSubstructures, self).__init__()
        layout = QtWidgets.QGridLayout()

        self.colorDictButton = {name: settings[name]["color"] for name in settings}
        self.colorDictButton["other"] = "#c0c0c0"

        self.buttonDict = {
            name: QtWidgets.QPushButton(settings[name]["name"]) for name in settings
        }
        self.buttonDict["other"] = QtWidgets.QLineEdit("Other?")

        self.buttonDict["other"].setAlignment(QtCore.Qt.AlignCenter)
        self.buttonDict["other"].returnPressed.connect(self.onEnterHit)

        self.transparentButton = QtWidgets.QPushButton()
        self.transparentButton.setStyleSheet(
            f"border: none;background-color:transparent"
        )
        self.transparentButton.pressed.connect(lambda: self.activateButton("other"))

        for button in self.buttonDict:
            self.buttonDict[button].setStyleSheet(
                f"border: 3px solid {self.colorDictButton[button]};color:{self.colorDictButton[button]}"
            )
            if button != "other":
                self.buttonDict[button].pressed.connect(
                    lambda color=button: self.activateButton(color)
                )

        nameList = [name for name in self.colorDictButton]
        self._initialLiability = nameList[0]

        for i, button in enumerate(self.buttonDict):
            layout.addWidget(self.buttonDict[button], i, 0)
        layout.addWidget(self.transparentButton, i, 0)
        layout.setAlignment(QtCore.Qt.AlignCenter)
        self.setLayout(layout)
        self.activateButton(self._initialLiability)
        self.buttonDict["other"].clearFocus()
        self.transparentButton.clearFocus()

    def onEnterHit(self):
        "Move outside of textbox when 'Enter' is hit."
        self.buttonDict["other"].clearFocus()

    def activateButton(self, buttonClicked: str):
        "Changes Color when button is clicked"
        self.buttonDict[buttonClicked].setStyleSheet(
            f"border: 3px solid {self.colorDictButton[buttonClicked]};background-color: {self.colorDictButton[buttonClicked]};color: #ECECEC"
        )

        # all other buttons are resetted to their original color scheme
        for button in self.buttonDict:
            if button != buttonClicked:
                self.buttonDict[button].setStyleSheet(
                    f"border: 3px solid {self.colorDictButton[button]};color:{self.colorDictButton[button]}; background:transparent"
                )
        if buttonClicked == "other":
            self.buttonDict["other"].setFocus()
        else:
            self.buttonDict["other"].clearFocus()
