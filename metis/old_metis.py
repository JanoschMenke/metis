#!/usr/bin/env python
from __future__ import print_function
from typing import List, Dict, Optional


# Import required modules
import sys, time, os, shutil
from PySide6.QtGui import *


from PySide6.QtWidgets import *
from PySide6.QtWidgets import QApplication
from PySide6 import QtCore, QtWidgets

import argparse
from metis.core import Backend
from metis.config.settings import load_settings, BaseConfig
from metis.ui import molwall, gui

from metis.core.data import sample_training_data
from utils import helper

from rdkit.Chem import AllChem as Chem
from . import PKGDIR
import os
from sys import platform


# TODO: File TRANSFER FOR SSH WITH SUBPROCESS
# TODO: Structure De Novo modeling better (reward model vs reward function.)


# The main window class
class Old_Metis(QtWidgets.QMainWindow):
    # Constructor function
    def __init__(
        self,
        settings: BaseConfig,
        output_folder: str,
        fileName=None,
        loglevel="WARNING",
    ):
        super(Old_Metis, self).__init__()
        # self.backend = Backend(settings, output_folder)
        self.settings = settings
        self.setSize()  # set size of main window
        self.loglevels = []

        self.initGUI(fileName=fileName)
        # ! ADAPT THIS
        self.setHighlightColors(self.backend.current_liability)
        self.editor.logger.setLevel(loglevel)

    def setSize(self):
        self.setMinimumSize(1500, 800)
        self.setMaximumSize(1920, 1080)

    def setHighlightColors(self, liability: str):
        #! CHECK THIS
        self.editor.last_selected_highlight_colour = self.backend.color_dict[liability]
        self.editor.selected_highlight_colour = self.backend.color_dict[liability]

    def initGUI(self, fileName=None):
        """
        Layout:
        ┌───────────┬────────────┐xx           xx
        │           │            │ x            x
        │   Editor  │ Evaluation │ x Mainlayout x
        │           │  Window    │ x            x
        │           │            │ x            x Globallayout
        ├───────────┴────────────┤xx            x
        │     Navigationbar      │              x
        └────────────────────────┘             xx
        """

        self.setWindowTitle("Metis")

        self.loadEvaluations()
        self.loadMolFile()
        if self.settings.ui.show_atom_contributions.render:
            self.loadExplanation()
        if self.settings.ui.show_reference_molecules.render:
            self.loadTrainMolecules()
        self.show()
        self.setFocus()

    def setCentralLayout(self):
        self.ui = gui.MetisGUI()
        self.center = QtWidgets.QWidget()
        self.center.setLayout(self.ui.setGlobalLayout())
        self.setCentralWidget(self.center)

        # self.fileName = fileName

    # INITIALIZE COMPONENTS
    # connects widgets and their buttons with functions

    # Function to setup status bar, central widget, menu bar, tool bar

    def loadMolFile(self):
        mol = Chem.MolFromSmiles(self.backend.current_smiles)
        self.editor.mol = mol

    def closeEvent(self, event):
        self.editor.logger.debug("closeEvent triggered")
        self.exitFile()
        event.ignore()

    def exitFile(self):
        self.saveEvaluations()
        numNotEvaluated = sum(self.backend.df.evaluated == 0)
        if numNotEvaluated > 0:
            text = f"""Thanks for participating. This will quit the application.\nYou haven't rated {numNotEvaluated} compounds.\nDo you want to finish the experiment and quit the application?
            """
        else:
            text = f"""Thanks for participating. This will quit the application.\nYou have rated all compounds.\nDo you want to finish the experiment and quit the application?
            """
        response = self.msgApp(
            "Confirmation",
            text,
        )
        if response == "Y":
            self.molWall.clearImages()
            self.backend.clear_temp_images()
            self.backend.save_substructure_dict()
            self.backend.save_final_dataset()
            exit(0)

    def changeTab(self):
        self.saveEvaluations()
        self.editor.selectedAtoms = []
        self.backend.currentTab = self.evaluationWindow.getTabName()
        print(f"New Tab {self.backend.currentTab}")
        # changes the color of the button to mathc the previous selection
        self.evaluationWindow.localEvaluations[
            self.backend.currentTab
        ].change_button_style(self.backend.current_liability)
        self.loadEvaluations()

    def saveSubstructDict(self):
        self.backend.save_substructure_dict()

    # Function to show Diaglog box with provided Title and Message
    def msgApp(self, title, msg):
        userInfo = QMessageBox.question(
            self, title, msg, QMessageBox.Yes | QMessageBox.No
        )
        if userInfo == QMessageBox.Yes:
            return "Y"
        if userInfo == QMessageBox.No:
            return "N"
        self.close()

    def openWindowFromNavBar(self, window: str):
        if window == "view":
            self.openMolViewer(molInfo=None)
        elif window == "unrated":
            self.findNextUnratedMol()
        elif window == "edit":
            self.openMolEditor()
        elif window == "send":
            self.train_and_update_model()
        elif window == "finish":
            self.exitFile()

    def findNextUnratedMol(self):
        stepSize = self.backend.next_unrated_mol
        if stepSize != 0:
            self.onMyToolBarButtonClick(stepSize)

    def openMolEditor(self) -> None:
        self.secondEditor.editor.mol = self.editor.mol
        self.secondEditor.original_mol = self.editor.mol
        self.secondEditor.show_and_raise()

    def openMolViewer(
        self,
        smiles: Optional[List] = None,
        molInfo: Optional[Dict] = None,
    ) -> None:
        self.saveEvaluations()
        if smiles is None:
            smiles, rating = self.backend.get_smiles_and_ratings()
        if (molInfo is None) & (smiles is not None):
            molInfo = [{"Evaluation": f"{i}"} for i in rating]
        self.molWall.initGUI(self.backend.iteration, smiles=smiles, molInfo=molInfo)

    def loadTrainMolecules(self):
        pixmap = self.backend.similar_actives
        self.svgWidget2.setPixmap(pixmap)

    def loadExplanation(self):
        pixmap = self.backend.explanation
        self.svgWidget1.setPixmap(pixmap)

    def train_and_update_model(self):
        self.saveEvaluations()
        response = self._checkMoleculesEval()
        if response == "Y":
            self.disable_inputs()
            smiles, dict_list = self.backend.update_models()
            if smiles is not None:
                self.openMolViewer(smiles, dict_list)
            else:
                self._updateInnerLoop()
                self.enable_inputs()

    def _checkMoleculesEval(self):
        response = "Y"
        if self.backend.all_molecules_evaluated:
            response = self.msgApp(
                "Confirmation",
                "Are you sure?",
            )
        return response

    def _updateInnerLoop(self):
        self.backend._update_inner_loop()
        self.initGUI(fileName=None)

    def startNextIteration(self):
        self.initGUI(fileName=None)
        self.enable_inputs()

        if self.backend.reached_max_iterations:
            self.navigation.sendButton.hide()

    def getEditedMolecule(self):
        self.backend.alternative_molecule = self.secondEditor.editor.mol

    def onMyToolBarButtonClick(self, direction: int):
        "Behavior when clicking 'Back' or 'Next' Button"
        self.saveEvaluations()
        self.backend.update_current_mol(direction)
        self.loadEvaluations()
        self.evaluationWindow.updateProperties(self.backend.property_labels)
        self.molCounter.updateCounter(self.backend.currentMolIndex)
        self.loadMolFile()
        if self.backend._atom_contributions:
            self.loadExplanation()
        if self.backend._reference_mols:
            self.loadTrainMolecules()

    def saveDataset(self, path: str = "testResults.csv") -> None:
        self.saveEvaluations()
        self.backend._save_dataset(path)

    def loadEvaluations(self):
        self.editor.selectedAtoms = self.backend.atom_selection
        ratingMolecule = self.backend.rating
        concernsMolecule = self.backend.concerns
        self.evaluationWindow.setRatings(ratingMolecule, concernsMolecule)
        self.loadTextProvidedInOther()

    def loadTextProvidedInOther(self):
        if self.backend.other_text != "":
            self.evaluationWindow.setOtherText(self.backend.other_text)
        else:
            self.evaluationWindow.resetOtherText()

    def saveEvaluations(self):
        """
        The `saveEvaluations` function updates the rating and concerns of a molecule in a dataframe, saves
        any additional text provided, and marks the molecule as evaluated.
        """
        self.gatherSelectedAtoms()
        self.backend.rating, self.backend.concerns = self.evaluationWindow.globalRatings
        self.backend.check_if_evaluated()
        self.saveTextProvidedInOther()

    def saveTextProvidedInOther(self):
        """
        The function saves the text provided in the "Other" field if it is not equal to "Other?".
        """
        if self.evaluationWindow.getTextInOther(self.backend.currentTab) != "Other?":
            self.backend.other_text = self.evaluationWindow.getTextInOther(
                self.backend.currentTab
            )

    def gatherSelectedAtoms(self):
        if len(self.editor.selectedAtoms) > 0:
            self.backend.atom_selection = self.editor.selectedAtoms
        else:
            self.backend.atom_selection = []

    def switchLiability(self, liability: str):
        print(f"Switching Liability: {self.backend.current_liability} -> {liability}")
        self.evaluationWindow.substructButtonDict["other"].clearFocus()

        self.gatherSelectedAtoms()
        self.backend.current_liability = liability
        self.setHighlightColors(liability=liability)

        self.editor.selectedAtoms = self.backend.atom_selection

    def disable_inputs(self):
        self.setEnabled(False)

    def enable_inputs(self):
        QtCore.QCoreApplication.processEvents()
        self.setEnabled(True)

    def clear_temp_images(self):
        self.backend.clear_temp_images()


def launch(loglevel="WARNING"):
    if (platform == "darwin") & (helper.is_faulty_pyside_version()):
        os.environ["QT_MAC_WANTS_LAYER"] = "1"
    "Function that launches the mainWindow Application"
    # Exception Handling
    parser = argparse.ArgumentParser("Metis")
    parser.add_argument("-f", "--file", default="setting.yml")
    parser.add_argument("-o", "--output", default="../results/")

    args = parser.parse_args()

    if not os.path.isfile(args.file):
        print(f"Settings not found or set: {args.file}")
        print()
        parser.print_help()
        sys.exit(1)

    myApp = QApplication(sys.argv)
    settings: BaseConfig = load_settings(args.file)
    mainWindow = Metis(
        settings,
        output_folder=args.output,
        fileName=args.file,
        loglevel=loglevel,
    )

    stylesheet = helper.read_and_replace_css(
        f"{PKGDIR}/resources/design/style.css", PKGDIR
    )
    myApp.setStyleSheet(stylesheet)
    myApp.exec_()

    sys.exit(0)
    """
    except NameError:
        print("Name Error:", sys.exc_info()[1])
    except SystemExit:
        print("Closing Window...")
    except Exception:
        print(sys.exc_info())
    """


if __name__ == "__main__":
    launch(loglevel="WARNING")
