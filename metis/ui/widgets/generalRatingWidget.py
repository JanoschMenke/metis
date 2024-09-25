from PySide6 import QtWidgets, QtCore, QtGui


class generalRating(QtWidgets.QWidget):
    """
    Combines Slider or Buttons into a single Widget for general rating.
    """

    valueChanged = QtCore.Signal(int)  # New signal for value changes

    def __init__(self, settings):
        super(generalRating, self).__init__()
        self.layout = QtWidgets.QVBoxLayout()
        self.layout.setAlignment(QtCore.Qt.AlignTop)
        if settings["slider"] == False:
            self.createButtons()
            self.layout.addLayout(self.buttonlayout)
        else:
            self.createSlider()
            self.layout.addLayout(self.buttonlayout)

        self.setLayout(self.layout)

    def onClick(self, value):
        """
        Handles button clicks for the button-based rating.

        Args:
            value (int): The index of the clicked button.
        """
        for i in range(len(self.buttonList)):
            if (i == value) and (value != self.currentPressed):
                self.buttonList[i].setStyleSheet(
                    f"background-color: {self.colors[i]} ;color: #ECECEC; border: 3px solid {self.colors[i]}"
                )
                self.currentPressed = i
                self.valueChanged.emit(i)  # Emit the new value
            elif (i == value) and (value == self.currentPressed):
                self.buttonList[i].setStyleSheet(
                    f"background-color: transparent;color: {self.colors[i]}; border: 3px solid {self.colors[i]}"
                )
                self.currentPressed = None
                self.valueChanged.emit(-1)  # Emit -1 for deselection
            else:
                self.buttonList[i].setStyleSheet(
                    f"background-color: transparent;color: {self.colors[i]}; border: 3px solid {self.colors[i]}"
                )

    def createButtons(self):
        """Creates the button-based rating interface."""
        self.buttonlayout = QtWidgets.QHBoxLayout()
        self.colors = ["#FF605C", "#F29F05", "#00CA4E"]
        ratings = ["Not at All", "Somewhat", "Very Much"]
        self.currentPressed = None
        self.buttonList = [QtWidgets.QPushButton(x) for x in ratings]

        for i in range(len(self.buttonList)):
            self.buttonList[i].setStyleSheet(
                f"background-color: transparent;color: {self.colors[i]};border: 3px solid {self.colors[i]}"
            )
            self.buttonList[i].pressed.connect(lambda value=i: self.onClick(value))
            self.buttonlayout.addWidget(self.buttonList[i])

    def createSlider(self):
        """Creates the slider-based rating interface."""
        self.currentPressed = None
        self.buttonList = None
        self.buttonlayout = QtWidgets.QHBoxLayout()
        self.slider = QtWidgets.QSlider()
        self.slider.setMaximumWidth(300)
        self.slider.setOrientation(QtCore.Qt.Horizontal)
        self.slider.setTickPosition(QtWidgets.QSlider.TicksBelow)
        self.slider.setTickInterval(10)
        self.slider.setMinimum(0)
        self.slider.setMaximum(100)

        self.slider.valueChanged.connect(self.changedValue)

        self.label = QtWidgets.QLabel("0")
        self.label.setFont(QtGui.QFont("Sanserif", 15))
        self.buttonlayout.addWidget(self.slider)
        self.buttonlayout.addWidget(self.label)

    def changedValue(self):
        """Handles value changes for the slider-based rating."""
        self.currentPressed = self.slider.value()
        self.label.setText(str(self.currentPressed))
        self.valueChanged.emit(self.currentPressed)  # Emit the new value

    def setSavedButton(self, value):
        """
        Sets the saved button state or slider value.

        Args:
            value (int): The value to set.
        """
        if self.buttonList is not None:
            for i in range(len(self.buttonList)):
                if i == value:
                    self.buttonList[i].setStyleSheet(
                        f"background-color: {self.colors[i]} ;color: #ECECEC; border: 3px solid {self.colors[i]}"
                    )
                else:
                    self.buttonList[i].setStyleSheet(
                        f"background-color: transparent;color: {self.colors[i]}; border: 3px solid {self.colors[i]}"
                    )
        else:
            if value is None:
                self.slider.setValue(0)
            else:
                self.slider.setValue(value)
        # self.valueChanged.emit(value if value is not None else 0)  # Emit the set value

    @property
    def eval_score(self):
        """
        Gets the current evaluation score.

        Returns:
            int: The current evaluation score.
        """
        return self.currentPressed

    @eval_score.setter
    def eval_score(self, value):
        """
        Sets the evaluation score.

        Args:
            value (int): The evaluation score to set.
        """
        self.setSavedButton(value)
        self.currentPressed = value
