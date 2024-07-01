from PySide6 import QtWidgets, QtCore, QtGui


class generalRating(QtWidgets.QWidget):
    """
    Combines Slider into a single Widget
    """

    def __init__(self, settings):
        super(generalRating, self).__init__()
        layout = QtWidgets.QVBoxLayout()
        layout.setAlignment(QtCore.Qt.AlignTop)
        if settings["slider"] == False:
            self.createButtons()
            layout.addLayout(self.buttonlayout)
        else:
            self.createSlider()
            layout.addLayout(self.buttonlayout)

        self.setLayout(layout)

    def onClick(self, value):
        for i in range(len(self.buttonList)):
            if (i == value) and (value != self.currentPressed):
                self.buttonList[i].setStyleSheet(
                    f"background-color: {self.colors[i]} ;color: #ECECEC; border: 3px solid {self.colors[i]}"
                )
                self.currentPressed = i
            elif (i == value) and (value == self.currentPressed):
                self.buttonList[i].setStyleSheet(
                    f"background-color: transparent;color: {self.colors[i]}; border: 3px solid {self.colors[i]}"
                )
                self.currentPressed = None
            else:
                self.buttonList[i].setStyleSheet(
                    f"background-color: transparent;color: {self.colors[i]}; border: 3px solid {self.colors[i]}"
                )

    def createButtons(self):
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
        self.currentPressed = self.slider.value()
        self.label.setText(str(self.currentPressed))

    def setSavedButton(self, value):
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
