import logging
from PySide6 import QtWidgets, QtCore

# Configure logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class selectLocalSubstructures(QtWidgets.QWidget):
    sanitizeSignal = QtCore.Signal(str)
    other_text_changed = QtCore.Signal(str)

    def __init__(self, settings):
        super(selectLocalSubstructures, self).__init__()
        logger.info("Initializing selectLocalSubstructures")

        layout = QtWidgets.QGridLayout()

        self.colorDictButton = {name: settings[name]["color"] for name in settings}
        self.colorDictButton["other"] = "#c0c0c0"

        self.buttonDict = {
            name: QtWidgets.QPushButton(settings[name]["name"]) for name in settings
        }
        self.buttonDict["other"] = CustomLineEdit("Other?")
        self.buttonDict["other"].editingFinished.connect(self.on_other_editing_finished)
        self.buttonDict["other"].focusLost.connect(self.on_other_editing_finished)

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
        self.set_active_button(self._initialLiability)
        self.buttonDict["other"].clearFocus()
        self.transparentButton.clearFocus()

        logger.debug(
            f"Initialization complete. Initial button: {self._initialLiability}"
        )

    def onEnterHit(self):
        "Move outside of textbox when 'Enter' is hit."
        logger.debug("Enter hit in 'other' field")
        self.buttonDict["other"].clearFocus()

    def activateButton(self, buttonClicked: str):
        "Changes Color when button is clicked"
        logger.debug(f"Substructure Button clicked: {buttonClicked}")
        logger.info(f"Switching to Substructure: {buttonClicked}")

        self.sanitizeSignal.emit(buttonClicked)
        self.set_active_button(buttonClicked)

    def change_button_style(self, buttonClicked: str):
        logger.debug(f"Changing button style for: {buttonClicked}")
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

    def get_active_button(self):
        """
        Returns the name of the currently activated button and the text in 'other' if it's selected.
        """
        logger.info(f"Getting active button: {self.currentlyActivated}")
        if self.currentlyActivated == "other":
            other_text = self.buttonDict["other"].text()
            logger.debug(f"'Other' is active with text: {other_text}")
            return self.currentlyActivated, other_text
        else:
            return self.currentlyActivated, None

    def set_active_button(self, button: str, other_text: str = None):
        """
        Sets the active button and updates the UI accordingly.
        If 'other' is selected, you can also set its text.
        """
        logger.debug(f"Setting active button to: {button}")
        if button not in self.buttonDict:
            logger.error(f"Invalid button name: {button}")
            raise ValueError(f"Invalid button name: {button}")

        self.currentlyActivated = button
        self.change_button_style(button)

        if button == "other" and other_text is not None:
            logger.debug(f"Setting 'other' text to: {other_text}")
            self.buttonDict["other"].setText(other_text)

        logger.debug(f"Active button set to: {button}")

    def on_other_editing_finished(self):
        # This slot will be called when focus is cleared from the QLineEdit
        current_text = self.buttonDict["other"].text()
        self.other_text_changed.emit(current_text)


class CustomLineEdit(QtWidgets.QLineEdit):
    focusLost = QtCore.Signal()

    def focusOutEvent(self, event):
        super().focusOutEvent(event)
        self.focusLost.emit()
