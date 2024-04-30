from PySide2.QtWidgets import QLabel
from PySide2 import QtWidgets, QtCore
from .ratingWidgets import (
    generalRating,
    TPProfile,
    QHLine,
    QVLine,
    multiSelect,
    navigationBar,
    selectLocalSubstructures,
    molCounter,
)
from PySide2.QtCore import Qt


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

    @property
    def initialLiability(self):
        name = self.getTabName()
        return self.evaluationWidget.localEvaluation[name]._initialLiability

    @property
    def substructButtonDict(self):
        name = self.getTabName()
        return self.evaluationWidget.localEvaluation[name].buttonDict

    @property
    def substructTransButton(self):
        name = self.getTabName()
        return self.evaluationWidget.localEvaluation[name].transparentButton

    @property
    def globalRatings(self):
        return self.evaluationWidget.globalInteractionWindow.getRatings()

    @property
    def textInOther(self):
        return self.substructButtonDict["other"].text()

    def getTabIndex(self) -> int:
        if self.uses_tab:
            idx = self.evaluationWidget.tab.currentIndex()
        else:
            idx = 0
        return idx

    @property
    def localEvaluations(self):
        return self.evaluationWidget.localEvaluation

    def getTabName(self) -> str:
        idx = self.getTabIndex()
        name = list(self.evaluationWidget.localEvaluation)[idx]
        return name

    def setRatings(self, ratingMolecule, concernsMolecule):
        self.evaluationWidget.globalInteractionWindow.setRatings(
            ratingMolecule, concernsMolecule
        )

    def resetOtherText(self):
        self.setOtherText("Other?")

    def setOtherText(self, text):
        self.substructButtonDict["other"].setText(text)


class stackedWidget(QtWidgets.QWidget):
    def __init__(self, settings):
        super(stackedWidget, self).__init__()

        pageLayout = QtWidgets.QGridLayout()

        self.substructureLayout = QtWidgets.QVBoxLayout()

        pageLayout.setSpacing(0)
        pageLayout.setContentsMargins(0, 0, 0, 0)
        self.globalInteractionWindow = generalEvaluation(settings)

        self.text1 = QLabel(
            "How much would you prioritize this molecule as a DRD2 binder?"
        )
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

        # self.tabwidget.currentChanged.connect(self.resetColorOnTabChange)
        self.setLayout(pageLayout)

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
        self.listwidget = multiSelect(
            settings["ui"]["global_properties"]["liabilities"]
        )
        self.listwidget.setMaximumSize(900, 100)

        layout.addWidget(self.listwidget, stretch=3)

        self.setLayout(layout)

    def getRatings(self):
        "returns the users rating of the molecule"
        ratingMolecule = self.sliderModule.currentPressed
        concernsMolecule = [
            int(self.listwidget.rowCheckboxes.itemAt(i).widget().checkState())
            for i in range(len(self.listwidget.labelsCheckbox))
        ]
        return ratingMolecule, concernsMolecule

    def setRatings(self, ratingMolecule, concernsMolecule):
        self.sliderModule.setSavedButton(ratingMolecule)
        self.sliderModule.currentPressed = ratingMolecule

        # for loop to set the checkbox to the saved state, default: off/0
        [
            self.listwidget.rowCheckboxes.itemAt(i)
            .widget()
            .setCheckState(self.__getQStates(concernsMolecule[i]))
            for i in range(len(concernsMolecule))
        ]

    def __getQStates(self, value):
        if value == 0:
            return Qt.CheckState.Unchecked
        elif value == 1:
            return Qt.CheckState.PartiallyChecked
        else:
            return Qt.CheckState.Checked
