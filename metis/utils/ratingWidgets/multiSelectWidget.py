from PySide6 import QtWidgets
from typing import List
import math


class multiSelect(QtWidgets.QWidget):
    """
    Combines multiple Checkboxes into a single Widget
    """

    def __init__(self, liabilities):
        super(multiSelect, self).__init__()

        mainLayout = QtWidgets.QVBoxLayout()
        self.labelsCheckbox = liabilities

        self.rowCheckboxes = QtWidgets.QGridLayout()
        layout_value = len(liabilities) / math.floor(len(liabilities) ** 0.5)
        if layout_value > 6:
            layout_value = 6
        for i in range(len(self.labelsCheckbox)):
            temp = QtWidgets.QCheckBox(self.labelsCheckbox[i])
            temp.setTristate(True)
            self.rowCheckboxes.addWidget(
                temp,
                i // layout_value,
                i % layout_value,
            )

        mainLayout.addLayout(self.rowCheckboxes)

        self.setLayout(mainLayout)


class customCheckBox(QtWidgets.QCheckBox):
    def __init__(self, text, parent=None):
        super().__init__(text, parent)
        self.current_state = 0

    def setChecked(self, value=None):
        if value is None:
            self.current_state += 1
        else:
            self.current_state = value
