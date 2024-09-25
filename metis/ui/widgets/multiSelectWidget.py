from PySide6 import QtWidgets, QtCore
from PySide6.QtCore import Qt, Signal
import math
from typing import Dict


class MultiSelect(QtWidgets.QWidget):
    """
    Combines multiple Checkboxes into a single Widget with a valueChanged signal.
    """

    valueChanged = Signal(dict)

    def __init__(self, liabilities):
        super(MultiSelect, self).__init__()

        mainLayout = QtWidgets.QVBoxLayout()
        self.labelsCheckbox = liabilities

        self.rowCheckboxes = QtWidgets.QGridLayout()
        layout_value = len(liabilities) / math.floor(len(liabilities) ** 0.5)
        if layout_value > 6:
            layout_value = 6
        self.checkbox_dict = {}  # Dictionary to store checkboxes by label
        for i, label in enumerate(self.labelsCheckbox):
            temp = QtWidgets.QCheckBox(label)
            temp.setTristate(True)
            temp.stateChanged.connect(self._on_checkbox_changed)
            self.rowCheckboxes.addWidget(
                temp,
                i // layout_value,
                i % layout_value,
            )
            self.checkbox_dict[label] = temp  # Store checkbox in dictionary

        mainLayout.addLayout(self.rowCheckboxes)

        self.setLayout(mainLayout)

    def __getQStates(self, value):
        if value == -1:
            return Qt.CheckState.Unchecked
        elif value == 0:
            return Qt.CheckState.PartiallyChecked
        else:
            return Qt.CheckState.Checked

    def _on_checkbox_changed(self):
        """
        Emits the valueChanged signal with the current checkbox states.
        """
        self.valueChanged.emit(self.eval_score)

    @property
    def eval_score(self) -> Dict[str, int]:
        """
        Returns the current state of all checkboxes.

        Returns:
            Dict[str, int]: A dictionary mapping checkbox labels to their states.
        """
        return {
            label: int(checkbox.checkState().value) - 1
            for label, checkbox in self.checkbox_dict.items()
        }

    @eval_score.setter
    def eval_score(self, score: Dict[str, int]) -> None:
        """
        Sets the state of checkboxes based on the provided score.

        Args:
            score (Dict[str, int]): A dictionary mapping checkbox labels to their desired states.
        """

        if not score:
            score = {label: -1 for label in self.labelsCheckbox}

        for label, value in score.items():
            if label in self.checkbox_dict:
                self.checkbox_dict[label].setCheckState(self.__getQStates(value))
