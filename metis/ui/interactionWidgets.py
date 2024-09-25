from PySide6.QtWidgets import QLabel
from PySide6 import QtWidgets, QtCore
from metis.ui.widgets import (
    generalRating,
    TPProfile,
    multiSelectWidget,
    selectLocalSubstructures,
)
from typing import Dict


class evaluationWindow(QtWidgets.QWidget):
    def __init__(self, settings):
        super(evaluationWindow, self).__init__()
        self.layout = QtWidgets.QVBoxLayout()
        self.layout.setSpacing(0)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.TPPDispay = TPProfile(settings)
        self.uses_tab = settings["ui"]["tab"]["render"]

        self.evaluationWidget = stackedWidget(settings)
        self.evaluationWidget.setMaximumSize(1000, 1000)

        self.layout.addWidget(self.TPPDispay)
        self.layout.addWidget(self.evaluationWidget)
        self.setLayout(self.layout)


class stackedWidget(QtWidgets.QWidget):
    tab_changed = QtCore.Signal(str)
    local_liability_changed = QtCore.Signal(str)
    text_in_other_changed = QtCore.Signal(str)

    def __init__(self, settings):
        super(stackedWidget, self).__init__()

        pageLayout = QtWidgets.QGridLayout()

        self.substructureLayout = QtWidgets.QVBoxLayout()

        pageLayout.setSpacing(0)
        pageLayout.setContentsMargins(0, 0, 0, 0)
        self.globalInteractionWindow = generalEvaluation(settings)

        self.text1 = QLabel("How much would you prioritize this molecule?")
        self.text1.setWordWrap(True)
        self.text1.setAlignment(QtCore.Qt.AlignTop)
        self.text2 = QLabel("What properties of the molecule raise concerns?")
        self.text2.setWordWrap(True)

        pageLayout.addWidget(self.text1, 0, 0)
        pageLayout.addWidget(self.globalInteractionWindow.sliderModule, 1, 0)
        pageLayout.addWidget(self.text2, 2, 0)
        pageLayout.addWidget(self.globalInteractionWindow.listwidget, 3, 0)

        if settings["ui"]["tab"]["render"] == True:
            self.tab = QtWidgets.QTabWidget()
            self.localEvaluation = {
                name: selectLocalSubstructures(
                    settings["ui"]["substructures"]["liabilities"]
                )
                for name in settings["ui"]["tab"]["tab_names"]
            }

            for name in self.localEvaluation:
                self.tab.addTab(self.localEvaluation[name], name)
            self.tab.setCurrentIndex(0)
        else:
            self.localEvaluation = {
                "general": selectLocalSubstructures(
                    settings["ui"]["substructures"]["liabilities"]
                )
            }

        for name in self.localEvaluation:
            self.localEvaluation[name].sanitizeSignal.connect(self.on_liability_change)
            self.localEvaluation[name].other_text_changed.connect(
                self.on_other_text_change
            )
            self.local_liability_changed.connect(
                self.localEvaluation[name].set_active_button
            )

        self.tab.currentChanged.connect(self.onTabChanged)

        self.text3 = QLabel("Any feedback on specific substructures?")
        self.text3.setWordWrap(True)
        self.text3.setAlignment(QtCore.Qt.AlignTop | QtCore.Qt.AlignCenter)
        pageLayout.addWidget(self.text3, 0, 1)
        if settings["ui"]["tab"]["render"] == True:
            pageLayout.addWidget(self.tab, 1, 1, 3, 1, QtCore.Qt.AlignTop)
        else:
            pageLayout.addWidget(
                self.localEvaluation["general"], 1, 1, 3, 1, QtCore.Qt.AlignTop
            )

        self.checkRender(settings)

        pageLayout.setRowStretch(0, 0)
        pageLayout.setRowStretch(1, 0)
        pageLayout.setRowStretch(2, 0)
        pageLayout.setRowStretch(3, 0)

        self.setLayout(pageLayout)

    def onTabChanged(self, index):
        new_tab_name = self.tab.tabText(index)
        print(f"Tab changed to index {new_tab_name}")
        # Emit our custom signal with the new tab name
        self.tab_changed.emit(new_tab_name)

    def on_liability_change(self, new_liability: str):
        self.local_liability_changed.emit(new_liability)

    def on_other_text_change(self, new_text: str):
        self.text_in_other_changed.emit(new_text)

    def checkRender(self, settings):
        if settings["ui"]["general"]["render"] == False:
            self.text1.hide()
            self.globalInteractionWindow.sliderModule.hide()
        if settings["ui"]["global_properties"]["render"] == False:
            self.text2.hide()
            self.globalInteractionWindow.listwidget.hide()
        if settings["ui"]["substructures"]["render"] == False:
            self.text3.hide()
            if settings["ui"]["tab"]["render"] == True:
                self.tab.hide()
            else:
                self.localEvaluation["general"].hide()


class generalEvaluation(QtWidgets.QWidget):
    """
    Custom Widget that combines all elements that are used to let the user evaluate compounds.
    """

    def __init__(self, settings, *args, **kwargs):
        super(generalEvaluation, self).__init__(*args, **kwargs)

        layout = QtWidgets.QVBoxLayout()
        self.sliderModule = generalRating(settings["ui"]["general"])
        self.sliderModule.setMaximumSize(800, 200)

        layout.addWidget(self.sliderModule)
        self.listwidget = multiSelectWidget.MultiSelect(
            settings["ui"]["global_properties"]["liabilities"]
        )
        self.listwidget.setMaximumSize(900, 100)

        layout.addWidget(self.listwidget, stretch=3)

        self.setLayout(layout)

    @property
    def eval_score(self) -> Dict[str, int]:
        "returns the users rating of the molecule"
        ratingMolecule = self.sliderModule.eval_score
        concernsMolecule = self.listwidget.eval_score
        concernsMolecule["general"] = ratingMolecule
        return concernsMolecule

    def setRatings(self, ratingMolecule: int, concernsMolecule: Dict[str, int]) -> None:
        self.sliderModule.eval_score = ratingMolecule
        self.listwidget.eval_score = concernsMolecule

        # for loop to set the checkbox to the saved state, default: off/0
