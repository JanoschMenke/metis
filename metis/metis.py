#!/usr/bin/env python
from __future__ import print_function
from typing import List, Dict, Optional


# Import required modules
import sys, time, os, shutil
from PySide2.QtGui import *
from PySide2.QtWidgets import *
from PySide2 import QtCore, QtGui, QtWidgets, QtSvg
from PySide2.QtSvg import QSvgWidget
from PySide2.QtCore import Qt, QThread, QThreadPool
import argparse
import pickle

from utils import (
    interactionWidgets,
    settings_validator,
    helper,
    counterfactual_mol_display,
    data,
    tutorial,
    molwall,
    analysis,
    tracker,
    draw,
)


import pandas as pd
import numpy as np
import yaml
from pathlib import Path
from rdeditor import rdeditor
from reinvent_connect import train_rf as trf
from reinvent_connect import to_reinvent as tr
import json
import os

# Import model
from utils import molEditWidget
from rdkit.Chem import Draw
from rdkit.Chem import AllChem as Chem
from utils.data import sample_training_data


# TODO: TESTS
# TODO: save final evaluations whe closing the window into results folder

parser = argparse.ArgumentParser("Metis")
parser.add_argument("-f", "--file", default="setting.yml")
args = parser.parse_args()

# Load settings as a dictionary from the YAML file
def load_settings(file_path: str) -> dict:
    with open(file_path, 'r') as f:
        settings = yaml.safe_load(f)
    return settings


# The main window class
class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, fileName=None, loglevel="WARNING"):
        super(MainWindow, self).__init__()

        try:
            [f.unlink() for f in Path("reinvent_connect/input_files/current_run/").glob("*") if f.is_file()]
        except:
            pass

        self.setSize()  # set size of main window
        self.loadFiles()  # load settings and set files paths
        self.initRFTrainer(self.settings['interactive_model'])  # Pass the model as a dict key
        self.setVariables()  # set properties used throughout the main application

        self.iteration = 0
        self.inner_iteration = 0
        self.track_substructure = tracker.metis_tracker()

        self.setMainEditor()  # setup main editor in which only atoms can be highlighted

        self.setMolecularEditor()  # the molecular Editor can be used to actually modify the molecule
        self.setHighlightColors()
        self.molWall = molwall.molwall()

        if self.settings['use_human_component']:  # Access settings through dict keys
            self.analyzer = analysis.analyser(
                path_initial_data=self.settings['data']['initial_path'],
                path_all_scaffolds="../data/all_topological_scaffolds.csv",
            )

        self._fileName = None
        self.initGUI(fileName=fileName)

        self.editor.logger.setLevel(loglevel)

        if self.settings['tutorial']:  # Check for tutorial setting
            self.tutorial = tutorial.Second(self)

        self.counterfactual_window = counterfactual_mol_display.CounterfactualWindow()

    def loadFiles(self):
        self.designPath = os.path.abspath(os.path.dirname(__file__)) + "/design/"
        self.settings = load_settings(args.file)  # Load settings as a dictionary
        self.load_settings(self.settings)

    def initRFTrainer(self, settings):
        self.RFTrainer = trf.trainRF(**settings)  # Pass settings as a dictionary

    def load_settings(self, settings):
        self._tutorial = settings['tutorial']
        self._atom_contributions = settings['ui']['show_atom_contributions']
        self._reference_mols = settings['ui']['show_reference_molecules']
        self._human_component = settings['use_human_component']
        self._debug = settings['debug']
        self._slider = settings['ui']['general']['slider']
        self._max_iter = settings['max_iterations']
        self._max_inner_loop_iter = settings['innerloop_iterations']
        self._activity_label = settings['activity_label']
        self._training_data_path = settings['interactive_model']['training_data_path']
        self._initial_qsar_model = pickle.load(open(settings['interactive_model']['qsar_model_path'], "rb"))
        self._num_molecules = settings['data']['num_molecules']

    def setSize(self):
        self.setMinimumSize(1500, 800)
        self.setMaximumSize(1920, 1080)

    def setVariables(self, initial: bool = True):
        self.loglevels = []
        self.df_result = []

        np.random.seed(23234)
        self.df, self.unrated_sample = data.loadData(
            self.settings, model=self.RFTrainer, initial=initial
        )

        if initial == True:
            self.results_path = data.createResultsFolder(
                f"../results/{self.settings['data']['run_name']}",
                debug=self._debug,
            )
        self.numProperties = len(self.unrated_sample)
        self.numEvaluatedMols = sum(self.df.evaluated)
        self.currentMolIndex = 0
        self.numMols = len(self.df)
        self.currentLiability = None
        self.df_list = []

    def setMainEditor(self):
        """
        The function sets the main editor widget to a fixed size of 600x600 pixels.
        """
        self.testTab1 = QtWidgets.QTabWidget()

        self.editor = molEditWidget.MolEditWidget()
        if self._atom_contributions or self._reference_mols:
            self.testTab1.addTab(self.editor, "Editor")
            self.testTab1.setFixedSize(600, 600)

        if self._atom_contributions:
            self.svgWidget1 = QLabel()
            self.testTab1.addTab(self.svgWidget1, "Explanation")
        if self._reference_mols:
            self.svgWidget2 = QLabel()
            self.testTab1.addTab(self.svgWidget2, "Similar Actives")

        else:
            self.editor.setFixedSize(600, 600)

    def setMolecularEditor(self):
        self.secondEditor = rdeditor()

    def setHighlightColors(self):
        self.editor.colorDict = data.createRGBColorDict(
            self.settings['ui']['substructures']['liabilities']
        )

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
        self.setWindowIcon(QIcon(self.designPath + "logo.png"))

        self.setMainLayout()  # responsible for holding molviewer and evaluation
        self.molCounter = interactionWidgets.molCounter(self.numMols)

        self.initNavigationBar()

        self.setCentralLayout()
        # Combine Mainlayout and Navbar

        self.filters = "MOL Files (*.mol *.mol);;Any File (*)"
        self.SetupComponents()
        self.infobar = QLabel("")
        self.myStatusBar.addPermanentWidget(self.infobar, 0)
        self.loadEvaluations()
        self.loadMolFile()
        self.clear_temp_images()
        if self._atom_contributions:
            self.loadExplanation()
        if self._reference_mols:
            self.loadTrainMolecules()

        self.show()
        self.setFocus()

    def setMainLayout(self):
        self.mainlayout = QtWidgets.QHBoxLayout()
        self.mainlayout.setContentsMargins(0, 0, 0, 0)
        if self._atom_contributions or self._reference_mols:
            self.mainlayout.addWidget(self.testTab1)  # self.editor
        else:
            self.mainlayout.addWidget(self.editor)
        self.initEvaluationWindow()
        self.currentLiability = self.evaluationWindow.initialLiability
        self.currentTab = "General"
        self.mainlayout.addWidget(self.evaluationWindow)

    def setGlobalLayout(self):
        globallayout = QtWidgets.QVBoxLayout()
        globallayout.addLayout(self.mainlayout)
        globallayout.addWidget(interactionWidgets.QHLine())
        globallayout.addWidget(self.navigation)
        globallayout.addWidget(self.molCounter)
        return globallayout

    def setCentralLayout(self):
        self.center = QtWidgets.QWidget()
        self.center.setLayout(self.setGlobalLayout())
        self.setCentralWidget(self.center)
        # self.fileName = fileName

    # INITIALIZE COMPONENTS
    # connects widgets and their buttons with functions

    def initNavigationBar(self):
        """
        The `initNavigationBar` function initializes the navigation bar and connects various buttons and
        shortcuts to their respective functions.
        """

        self.navigation = interactionWidgets.navigationBar(
            self.settings['ui']['navigationbar']
        )
        self.navigation.backButton.clicked.connect(
            lambda: self.onMyToolBarButtonClick(-1)
        )
        QtWidgets.QShortcut(
            QtGui.QKeySequence("left"),
            self.navigation.backButton,
            lambda: self.onMyToolBarButtonClick(-1),
        )

        QtWidgets.QShortcut(
            QtGui.QKeySequence("right"),
            self.navigation.nextButton,
            lambda: self.onMyToolBarButtonClick(+1),
        )
        self.navigation.nextButton.clicked.connect(
            lambda: self.onMyToolBarButtonClick(+1)
        )

        self.navigation.downloadButton.clicked.connect(self.download)
        self.navigation.unratedMolecule.clicked.connect(self.findNextUnratedMol)
        self.navigation.editButton.clicked.connect(
            lambda: self.openMolEditor(molInfo=None)
        )

        self.navigation.sendButton.clicked.connect(self.getCounterfactual)
        self.secondEditor.saveAction.triggered.connect(self.getEditedMolecule)

        # Add the new button for molecule viewer
        self.navigation.moleculeViewerButton.clicked.connect(self.openMolViewer)

    def initEvaluationWindow(self):
        """
        The function initializes an evaluation window and sets its minimum size, updates molecular
        properties, and connects buttons to a switchLiability function.
        """
        self.evaluationWindow = interactionWidgets.evaluationWindow(
            self.settings
        )
        self.evaluationWindow.setMinimumSize(600, 600)
        self.evaluationWindow.TPPDispay.updateMolProperties(
            self.df.getPropertyLabels(self.currentMolIndex)
        )
        for name in self.evaluationWindow.localEvaluations:
            for button in self.evaluationWindow.localEvaluations[name].buttonDict:
                if button != "other":
                    self.evaluationWindow.localEvaluations[name].buttonDict[
                        button
                    ].pressed.connect(
                        lambda color=button: self.switchLiability(liability=color)
                    )
                else:
                    self.evaluationWindow.localEvaluations[
                        name
                    ].transparentButton.pressed.connect(
                        lambda color=button: self.switchLiability(liability=color)
                    )
        if self.evaluationWindow.uses_tab:
            self.evaluationWindow.evaluationWidget.tab.currentChanged.connect(
                self.changeTab
            )

    # Function to setup status bar, central widget, menu bar, tool bar
    def SetupComponents(self):
        self.myStatusBar = QStatusBar()

        self.setStatusBar(self.myStatusBar)
        self.myStatusBar.showMessage("Ready", 10000)
        self.CreateActions()

    def loadMolFile(self):
        mol = Chem.MolFromSmiles(self.df.SMILES.iloc[self.currentMolIndex])
        self.editor.mol = mol

    def openFile(self):
        self.fileName, self.filterName = QFileDialog.getOpenFileName(
            self, caption="Open MOL file", filter=self.filters
        )
        self.loadMolFile(self.fileName)

    def closeEvent(self, event):
        self.editor.logger.debug("closeEvent triggered")
        self.exitFile()
        event.ignore()

    def exitFile(self):
        self.saveEvaluations()
        numNotEvaluated = sum(self.df.evaluated == 0)
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
            self.saveDataset()
            self.molWall.clearImages()
            self.clear_temp_images()
            self.saveSubstructDict()

            if len(self.df_list) > 0:
                self.df = pd.concat(
                    self.df_list, axis=0, ignore_index=True
                ).reset_index(drop=True)
            self.saveDataset(f"{self.results_path}/final_evaluation_data.csv")

            if os.path.isfile("../data/scaffold_memory.csv"):
                shutil.copy2(
                    "../data/scaffold_memory.csv",
                    f"{self.results_path}/final_scaffold_memory.csv",
                )
            exit(0)

    def changeTab(self):
        self.saveEvaluations()
        self.currentTab = self.evaluationWindow.getTabName()
        self.evaluationWindow.localEvaluations[self.currentTab].activateButton(
            self.currentLiability
        )
        self.loadEvaluations()
        self.loadMolFile()

    def mousePressEvent(self, event):
        self.evaluationWindow.substructButtonDict["other"].clearFocus()

    def saveSubstructDict(self):
        self.track_substructure.save_substruct(f"{self.results_path}/substructure.json")

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

    def setLogLevel(self):
        loglevel = self.sender().objectName().split(":")[-1].upper()
        self.editor.logger.setLevel(loglevel)

    def findNextUnratedMol(self):
        indexNotEvaluated = self.df[self.df.evaluated == False].index.values
        if len(indexNotEvaluated) != 0:
            differences = indexNotEvaluated - self.currentMolIndex
            print(differences)

            if np.sum(differences > 0) > 0:
                stepSize = np.min(differences[differences > 0])
            else:
                stepSize = np.min(differences[differences < 0])
            self.onMyToolBarButtonClick(stepSize)

    def openMolEditor(
        self,
        smiles: Optional[List] = None,
        molInfo: Optional[Dict] = None,
    ) -> None:
        # self.secondEditor.editor.mol = self.editor.mol
        # self.secondEditor.editor.originalMol = self.editor.mol
        # self.secondEditor.show()
        self.saveEvaluations()
        if smiles is None:
            if self.settings['ui']['general']['slider'] == True:
                smiles = self.df.SMILES[self.df.rating != 50].values.tolist()
                rating = self.df.rating[self.df.rating != 50].values.tolist()
            else:
                smiles = self.df.SMILES[~self.df.rating.isna()].values.tolist()
                rating = self.df.rating[~self.df.rating.isna()].values.tolist()
        if (molInfo is None) & (smiles is not None):
            molInfo = [{"Evaluation": f"{i}"} for i in rating]
        self.molWall.initGUI(self.iteration, smiles=smiles, molInfo=molInfo)

    def openMolViewer(self):
        # self.secondEditor.editor.mol = self.editor.mol
        # self.secondEditor.editor.originalMol = self.editor.mol
        # self.secondEditor.show()
        self.saveEvaluations()
        smiles, molInfo = sample_training_data(
            path_training_data=self._training_data_path,
            current_smiles_query=self.df.SMILES.iloc[self.currentMolIndex],
        )
        self.molWall.initGUI(self.iteration, smiles=smiles, molInfo=molInfo)

    def loadTrainMolecules(
        self, smiles: Optional[List] = None, molInfo: Optional[Dict] = None
    ):
        self.saveEvaluations()
        if not os.path.isfile(f"mostSimilarActives{self.currentMolIndex}.png"):
            draw.save_most_similar_active_map(
                self._training_data_path,
                self.df.SMILES.iloc[self.currentMolIndex],
                save_name=f"mostSimilarActives{self.currentMolIndex}.png",
            )
        # Display the saved image using the path
        pixmap = QPixmap(
            f"{self._temp_image_folder}mostSimilarActives{self.currentMolIndex}.png"
        )
        pixmap.scaled(
            600,
            600,
            QtCore.Qt.KeepAspectRatio,
            mode=QtCore.Qt.SmoothTransformation,
        )
        self.svgWidget2.setPixmap(pixmap)

    def loadExplanation(self):
        if not os.path.isfile(
            f"utils/temp_images/atomContribution{self.currentMolIndex}.png"
        ):
            svg = draw.save_similarity_map(
                self.editor._mol,
                self._initial_qsar_model,
                save_name=f"atomContribution{self.currentMolIndex}.png",
                ecfp_settings=self.settings['interactive_model']['ECFP'],
            )
        pixmap = QPixmap(
            f"{os.getcwd()}/utils/temp_images/atomContribution{self.currentMolIndex}.png"
        )
        pixmap = pixmap.scaled(
            600,
            600,
            QtCore.Qt.KeepAspectRatio,
            mode=QtCore.Qt.SmoothTransformation,
        )
        self.svgWidget1.setPixmap(pixmap)

    def getCounterfactual(self):
        self.saveEvaluations()
        response = "Y"
        print(sum(self.df.evaluated))
        if sum(self.df.evaluated) != self._num_molecules:
            response = self.msgApp(
                "Confirmation",
                "Are you sure?",
            )
        if response == "Y":
            self.disable_inputs()

            self.df_list.append(self.df)
            if self.RFTrainer.use_oracle_for_prediction:
                smiles = self._getSmilesAndYScoresOracle()
                print(smiles)
            else:
                smiles, yScores = self._getSmilesAndYScores()

            # selectedAtoms = self.df.at[self.currentMolIndex, "uglyAtomSelection"]
            if len(smiles) > 0:  # check if atoms are selected
                if self.RFTrainer.use_oracle_for_prediction:
                    complete_mcc, new_mcc = self.RFTrainer.update_model(smiles)
                else:
                    complete_mcc, new_mcc = self.RFTrainer.update_model(smiles, yScores)

                print(f"F1 Complete Data: {complete_mcc} - F1 New Molecules: {new_mcc}")
                print(self.RFTrainer.df)
                self.df_result.append(self.RFTrainer.df)

                if self._shouldProcessInnerLoop():
                    self._updateModelAndSend(smiles=smiles)
                else:
                    self._updateInnerLoop()
                    self.enable_inputs()

    def _updateModelAndSend(self, smiles):
        self.df = pd.concat(self.df_list, axis=0, ignore_index=True).reset_index(
            drop=True
        )
        data.createResultsFolder(f"{self.results_path}/iteration_{self.iteration}")
        self.saveDataset(
            f"{self.results_path}/iteration_{self.iteration}/evaluation_data.csv"
        )

        self._copyScaffoldMemory()

        results_df = self._concatResultsDF()
        self._saveResultsCSV(results_df)

        dict_list = self._createDictList(results_df)

        self.RFTrainer.save_model("reinvent_connect/input_files/current_run/Model.pkl")
        self.RFTrainer.save_model(
            f"{self.results_path}/iteration_{self.iteration}/Model.pkl"
        )

        user_scoring_function = self._getUserScoringFunction(results_df)
        self.sent2Reinvent(user_scoring_function)
        self.openMolEditor(results_df.smiles, dict_list)

    def _copyScaffoldMemory(self):
        shutil.copy2(
            self.settings['data']['path'],
            f"{self.results_path}/iteration_{self.iteration}/scaffold_memory.csv",
        )

    def _concatResultsDF(self):
        results_df = pd.concat(self.df_result, axis=0, ignore_index=True)
        results_df["true_initial_pred"] = self.df.loc[
            self.df.SMILES.isin(results_df.smiles), self._activity_label
        ].values.tolist()

        return results_df

    def _saveResultsCSV(self, results_df):
        results_df.to_csv(
            f"{self.results_path}/iteration_{self.iteration}/oracle_results.csv",
            index=False,
        )

    def _getUserScoringFunction(self, results_df):
        if self._human_component:
            self.track_substructure.load_data(self.df)
            return self.track_substructure.create_scoring_function()
        else:
            return None

    def _getSmilesAndYScoresOracle(self):
        if self._slider == True:
            smiles = self.df.SMILES.values
        else:
            smiles = self.df.SMILES[self.df.rating == 2].values

        return smiles

    def _getSmilesAndYScores(self):
        if self._slider == True:
            smiles = self.df.SMILES.values
            y_scores = self.df.rating.values / 100
        else:
            smiles = self.df.SMILES[self.df.rating.isin([0, 2])].values
            y_scores = np.maximum(
                0, self.df.rating[self.df.rating.isin([0, 2])].values - 1
            )

        return smiles, y_scores

    def _createDictList(self, results_df):
        dict_list = []
        for i in range(results_df.shape[0]):
            temp_dict = {"Iteration": self.iteration + 1}
            temp_dict["Target"] = (
                "Active" if results_df.target.iloc[i] == 1 else "Inactive"
            )
            temp_dict["Model Prediction"] = (
                "Active" if results_df.iloc[i, 3] == 1 else "Inactive"
            )
            dict_list.append(temp_dict)
        return dict_list

    def _shouldProcessInnerLoop(self):
        return (
            self._max_inner_loop_iter is None
            or (self._max_inner_loop_iter - 1) == self.inner_iteration
        )

    def _updateInnerLoop(self):
        print(f"Inner Loop: {self.inner_iteration}")
        self.inner_iteration += 1

        np.random.seed(23234)
        self.df, self.unrated_sample = data.loadData(
            self.settings, model=self.RFTrainer, initial=False
        )
        self.numProperties = len(self.unrated_sample)
        self.numEvaluatedMols = sum(self.df.evaluated)
        self.currentMolIndex = 0
        self.numMols = len(self.df)
        self.currentLiability = None

        self.initGUI(fileName=None)

    def startNextIteration(self):
        shutil.copy2(
            "reinvent_connect/input_files/current_run/Agent.ckpt",
            f"{self.results_path}/iteration_{self.iteration}/Agent.ckpt",
        )
        self.inner_iteration = 0
        self.enable_inputs()
        self.setVariables(initial=False)
        self.initGUI(fileName=None)
        self.iteration += 1

        if self._max_iter is not None and self.iteration == self._max_iter:
            self.navigation.sendButton.hide()

    def trackExperiment(self):
        if self.iteration == 0:
            entire_df = pd.read_csv(self.settings['data']['initial_path'])
        else:
            entire_df = pd.read_csv(self.settings['data']['path'])
        _, out = self.analyzer.trackFoundScaffolds(self.df, entire_df)

    def sent2Reinvent(
        self,
        user_scoring_function: List = [],
    ):
        """
        The function starts a worker thread using QThreadPool which starts the
        Reinvent run on the server.
        """
        self.thread = QThreadPool()
        self.worker = tr.Worker(self._activity_label, user_scoring_function)
        self.worker.signals.finished.connect(self.startNextIteration)
        self.thread.start(self.worker)

        # toKeep = tr.gen_counter_factual(smiles, selectedAtoms)
        # return toKeep

    def getEditedMolecule(self):
        try:
            self.df.loc[self.currentMolIndex, "alternativeMolecule"] = Chem.MolToSmiles(
                self.secondEditor.editor.mol
            )
        except:
            print("Invalid Molecule!")

    def onMyToolBarButtonClick(self, direction: int):
        "Behavior when clicking 'Back' or 'Next' Button"

        self.saveEvaluations()
        self.currentMolIndex = (self.currentMolIndex + direction) % self.df.shape[0]
        self.loadEvaluations()
        self.evaluationWindow.TPPDispay.updateMolProperties(
            self.df.getPropertyLabels(self.currentMolIndex)
        )
        self.molCounter.updateCounter(self.currentMolIndex)
        self.loadMolFile()
        if self._atom_contributions:
            self.loadExplanation()
        if self._reference_mols:
            self.loadTrainMolecules()

    def saveDataset(self, path: str = "../data/testResults.csv") -> None:
        self.saveEvaluations()
        self.df.to_csv(path, index=False)

    def download(self):
        "Behavior when clicking 'Download' Button"
        self.exitFile()

    def loadEvaluations(self):
        self.loadSelectedMolecules()
        ratingMolecule = self.df.getRating(self.currentMolIndex)
        concernsMolecule = self.df.getGlobalLiab(self.currentMolIndex)
        self.evaluationWindow.setRatings(ratingMolecule, concernsMolecule)
        self.loadTextProvidedInOther()

    def loadSelectedMolecules(self):
        """
        The function "loadSelectedMolecules" loads selected atoms from a dataframe and resets the download
        button.
        """
        self.editor.selectedAtoms = self.df.loc[
            self.currentMolIndex,
            f"{self.currentLiability}AtomSelection_{self.currentTab}",
        ]
        # We use this function to also reset the download button
        # as it is called always when another button is pushed
        self.navigation.resetDownloadText()

    def gatherSelectedAtoms(self):
        # check if atoms were selected to be stored. It should work without it
        # but its is buggy for some reason

        if len(self.editor.selectedAtoms) > 0:
            self.df.setSelectedAtoms(
                self.currentMolIndex,
                [self.currentLiability, self.currentTab],
                self.editor.selectedAtoms,
                self.editor.mol,
            )
        self.editor.selectedAtoms = []

    def saveTextProvidedInOther(self):
        """
        The function saves the text provided in the "Other" field if it is not equal to "Other?".
        """
        if self.evaluationWindow.textInOther != "Other?":
            self.df.saveOtherText(
                self.currentMolIndex, self.evaluationWindow.textInOther
            )

    def loadTextProvidedInOther(self):
        if self.df.getOtherText(self.currentMolIndex, self.currentTab) != "":
            self.evaluationWindow.setOtherText(
                self.df.getOtherText(self.currentMolIndex, self.currentTab)
            )
        else:
            self.evaluationWindow.resetOtherText()

    def saveEvaluations(self):
        """
        The `saveEvaluations` function updates the rating and concerns of a molecule in a dataframe, saves
        any additional text provided, and marks the molecule as evaluated.
        """
        self.gatherSelectedAtoms()
        ratingMolecule, concernsMolecule = self.evaluationWindow.globalRatings
        self.df.at[self.currentMolIndex, "rating"] = ratingMolecule
        self.df.loc[
            self.currentMolIndex,
            self.settings['ui']['global_properties']['liabilities'],
        ] = concernsMolecule
        self.saveTextProvidedInOther()

        currentProperties = self.df.iloc[
            self.currentMolIndex, -(self.numProperties + 1) : -1
        ]

        if self.unrated_sample.equals(currentProperties) == False:
            self.df.at[self.currentMolIndex, "evaluated"] = 1

    # Function to create actions for menus and toolbars
    def CreateActions(self):
        self.loglevelactions = {}
        for key in self.loglevels:
            self.loglevelactions[key] = QAction(
                key,
                self,
                statusTip="Set logging level to %s" % key,
                triggered=self.setLogLevel,
                objectName="loglevel:%s" % key,
            )

    def switchLiability(self, liability: str):
        "moves cursor out of 'other' box."
        "This is somewhat redundent as it also happens within the localEvaluation Widget"
        "But in case we switch to global we also want to move out of the textbox"
        print("switching liability from", self.currentLiability, "to", liability)
        if liability != "other":
            self.evaluationWindow.substructButtonDict["other"].clearFocus()

        self.gatherSelectedAtoms()
        self.currentLiability = liability
        self.editor.drawColor = liability
        self.loadSelectedMolecules()

    def disable_inputs(self):
        self.setEnabled(False)

    def enable_inputs(self):
        QtCore.QCoreApplication.processEvents()
        self.setEnabled(True)

    def clear_temp_images(self):
        self._temp_image_folder = f"{os.getcwd()}/utils/temp_images/"
        if os.path.exists(self._temp_image_folder):
            shutil.rmtree(self._temp_image_folder)
        os.mkdir(self._temp_image_folder)


def launch(loglevel="WARNING"):
    "Function that launches the mainWindow Application"
    # Exception Handling
    try:
        myApp = QApplication(sys.argv)
        try:
            mainWindow = MainWindow(fileName=sys.argv[1], loglevel=loglevel)
        except:
            mainWindow = MainWindow(loglevel=loglevel)

        stylesheet = open("design/style.css").read()
        myApp.setStyleSheet(stylesheet)
        myApp.exec_()

        sys.exit(0)
    except NameError:
        print("Name Error:", sys.exc_info()[1])
    except SystemExit:
        print("Closing Window...")
    except Exception:
        print(sys.exc_info()[1])


if __name__ == "__main__":
    launch(loglevel="WARNING")
