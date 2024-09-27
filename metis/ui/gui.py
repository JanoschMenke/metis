from PySide6 import QtWidgets, QtCore
from rdeditor import molEditWidget
from PySide6.QtCore import Signal, Slot
from PySide6.QtGui import *
from metis.config.settings import BaseConfig
from metis.ui import widgets, molwall, interactionWidgets
from metis.utils import rdeditor_wrapper, helper
from rdkit import Chem
from typing import Optional


class MetisGUI(QtWidgets.QWidget):
    alternative_molecule_saved = Signal(str)

    def __init__(self, settings: BaseConfig):
        super(MetisGUI, self).__init__()
        self.settings = settings
        self.init_components()
        self.set_main_layout()
        self.connect_signals()

    def init_components(self):
        """Initialize all GUI components."""
        self.molWall = molwall.molwall()
        self.set_main_editor()
        self.set_molecular_editor()
        self.init_navigationbar()
        self.init_evaluation_window()
        self.create_color_dict()

    def set_main_editor(self):
        """Set up the main editor widget."""
        self.testTab1 = QtWidgets.QTabWidget()
        self.editor = molEditWidget.MolEditWidget()
        self.editor.action = "Select"

        if (
            self.settings.ui.show_atom_contributions.render
            or self.settings.ui.show_reference_molecules.render
        ):
            self.testTab1.addTab(self.editor, "Editor")
            self.testTab1.setFixedSize(600, 600)

            if self.settings.ui.show_atom_contributions.render:
                self.svgWidget1 = QtWidgets.QLabel()
                self.testTab1.addTab(self.svgWidget1, "Explanation")

            if self.settings.ui.show_reference_molecules.render:
                self.svgWidget2 = QtWidgets.QLabel()
                self.testTab1.addTab(self.svgWidget2, "Similar Actives")
        else:
            self.editor.setFixedSize(600, 600)

    def connect_signals(self):
        self.evaluationWindow.evaluationWidget.local_liability_changed.connect(
            self.set_highlight_colors
        )
        self.second_editor.saveAction.triggered.connect(self.alternative_moelcule)

    def set_main_layout(self):
        self.mainlayout = QtWidgets.QHBoxLayout()
        self.mainlayout.setContentsMargins(0, 0, 0, 0)
        if (
            self.settings.ui.show_atom_contributions.render
            or self.settings.ui.show_reference_molecules.render
        ):
            self.mainlayout.addWidget(self.testTab1)  # self.editor
        else:
            self.mainlayout.addWidget(self.editor)

        self.mainlayout.addWidget(self.evaluationWindow)

    def create_color_dict(self):
        self.color_dict = helper.create_color_dictionary(
            self.settings.ui.substructures.liabilities
        )
        self.set_highlight_colors(
            next(iter(self.settings.ui.substructures.liabilities))
        )

    def setGlobalLayout(self):
        globallayout = QtWidgets.QVBoxLayout()
        globallayout.addLayout(self.mainlayout)
        globallayout.addWidget(widgets.QHLine())
        globallayout.addWidget(self.navigation)
        globallayout.addWidget(self.molCounter)
        return globallayout

    def init_navigationbar(self):
        self.molCounter = widgets.molCounter(self.settings.data.num_molecules)
        self.navigation = widgets.navigationBar(self.settings.ui.navigationbar.dict())
        # if no reward model is specified in settings, remove sendButton
        if self.settings.reward_model is None:
            self.navigation.sendButton.hide()

    def init_evaluation_window(self):
        self.evaluationWindow = interactionWidgets.evaluationWindow(
            self.settings.dict()
        )
        self.evaluationWindow.setMinimumSize(600, 600)

    def set_molecular_editor(self):
        self.second_editor = rdeditor_wrapper.MetisEditor(loglevel="Warning")

    def update_molecule(self, mol: Chem.Mol):
        self.editor.mol = mol

    def update_atom_selection(self, atoms: list):
        self.editor.selectedAtoms = atoms

    def set_highlight_colors(self, liability: str):
        self.editor.last_selected_highlight_colour = self.color_dict[liability]
        self.editor.selected_highlight_colour = self.color_dict[liability]

    def update_text_in_other(self, tab_name: str, text: str) -> None:
        self.evaluationWindow.evaluationWidget.localEvaluation[tab_name].buttonDict[
            "other"
        ].setText(text)

    def msg_app(self, title, msg):
        userInfo = QtWidgets.QMessageBox.question(
            self, title, msg, QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No
        )
        if userInfo == QtWidgets.QMessageBox.Yes:
            return "Y"
        if userInfo == QtWidgets.QMessageBox.No:
            return "N"
        self.close()

    def exit_msg(self, num_not_evaluated: Optional[int] = None) -> str:
        if num_not_evaluated:
            text = f"""Thanks for participating. This will quit the application.\nYou haven't rated {num_not_evaluated} compounds.\nDo you want to finish the experiment and quit the application?
            """
        else:
            text = f"""Thanks for participating. This will quit the application.\nYou have rated all compounds.\nDo you want to finish the experiment and quit the application?
            """
        response = self.msg_app(
            "Confirmation",
            text,
        )
        return response

    def train_confirmation(self):

        response = self.msg_app(
            "Confirmation",
            "Are you sure?",
        )
        return response

    def load_similarity_map(self, pixmap):
        self.svgWidget2.setPixmap(pixmap)

    def load_explanation(self, pixmap):
        self.svgWidget1.setPixmap(pixmap)

    def disable_inputs(self):
        self.setEnabled(False)

    def enable_inputs(self):
        QtCore.QCoreApplication.processEvents()
        self.setEnabled(True)

    def open_second_editor(self):
        self.second_editor.editor.mol = self.editor.mol
        self.second_editor.original_mol = self.editor.mol
        self.second_editor.show_and_raise()

    def alternative_moelcule(self):
        self.alternative_molecule_saved.emit(
            Chem.MolToSmiles(self.second_editor.editor.mol)
        )
