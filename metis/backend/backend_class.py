from metis.utils import (
    helper,
    data,
    counterfactual_mol_display,
    tutorial,
    molwall,
    analysis,
    tracker,
    draw,
)
from metis.reinvent_connect import to_reinvent as tr
from metis.backend import settings_validator

from PySide6.QtCore import Qt, QThread, QThreadPool, Signal
import os
from os.path import join
import yaml
import pickle
from metis.reinvent_connect import train_rf as trf
import numpy as np
from rdkit import Chem
import shutil
import pandas as pd
from typing import List, Dict
from metis import PKGDIR


class Backend:
    def __init__(self, settings, results_folder: str):
        helper.clear_current_files(
            f"{PKGDIR}/reinvent_connect/input_files/current_run/"
        )
        self.loadFiles(settings)
        if self.settings.reward_model is not None:
            self.initRFTrainer(self.settings.reward_model)
        if self.settings.de_novo_model is not None:
            self.denovo_runner = tr.DeNovoRunner(self.settings.de_novo_model)
        else:
            self.denovo_runner = None

        self.setVariables(results_folder=results_folder, initial=True)
        self.setGlobalVariables()
        self.clear_temp_images()

    def loadFiles(self, settings):

        self.designPath = f"{PKGDIR}/design/"
        self.settings = settings_validator.BaseConfig(**yaml.safe_load(open(settings)))
        self.load_settings(self.settings)

    def load_settings(self, settings):
        self._tutorial = settings.tutorial
        self._human_component = (
            settings.de_novo_model.use_human_scoring_func
            if settings.de_novo_model is not None
            else False
        )
        self._debug = settings.debug
        self._slider = settings.ui.general.slider
        self._max_iter = settings.max_iterations
        self._max_inner_loop_iter = settings.innerloop_iterations
        self._activity_label = settings.activity_label
        self._num_molecules = settings.data.num_molecules

        self._reference_mols = settings.ui.show_reference_molecules.render
        self._training_data_path = settings.ui.show_reference_molecules.path
        self._atom_contributions = settings.ui.show_atom_contributions.render
        if self._atom_contributions:
            self._initial_qsar_model = pickle.load(
                open(settings.ui.show_atom_contributions.path, "rb")
            )

    def initRFTrainer(self, settings):
        self.RFTrainer = trf.trainRF(settings)

    def setVariables(self, results_folder: str = "", initial: bool = True):
        self.df_result = []
        self.df_list = []
        self.color_dict = data.createRGBColorDict(
            self.settings.ui.substructures.liabilities
        )
        if initial == True:
            self.results_path = data.createResultsFolder(
                join(results_folder, self.settings.data.run_name),
                debug=self._debug,
            )

        self.set_per_round_variables(initial=initial)
        self.clear_temp_images(files_only=True)

    def set_per_round_variables(self, initial: bool = True):
        np.random.seed(self.settings.seed)
        self.df, self.unrated_sample = data.loadData(
            self.settings.dict(), initial=initial
        )

        self.numProperties = len(self.unrated_sample)
        self.numEvaluatedMols = sum(self.df.evaluated)
        self.currentMolIndex = 0
        self.numMols = len(self.df)
        self._current_liability = None

    def setGlobalVariables(self):
        self.next_iteration_signal = tr.WorkerSignals()
        self.currentTab = "General"
        self.iteration = 0
        self.inner_iteration = 0
        self.track_substructure = tracker.metis_tracker()
        self._fileName = None
        if self._tutorial == True:
            self.tutorial = tutorial.Second(self)
        self.counterfactual_window = counterfactual_mol_display.CounterfactualWindow()

    def init_next_iteration(self):
        if self.settings.de_novo_model is not None:
            shutil.copy2(
                f"{PKGDIR}/reinvent_connect/input_files/current_run/Agent.ckpt",
                f"{self.results_path}/iteration_{self.iteration}/Agent.ckpt",
            )
        self.inner_iteration = 0
        self.setVariables(initial=False)
        self.iteration += 1
        self.next_iteration_signal.finished.emit()

    def clear_temp_images(self, files_only: bool = False):
        self._temp_image_folder = f"{PKGDIR}/utils/temp_images/"

        if os.path.exists(self._temp_image_folder):
            if files_only:
                self.remove_files_only(self._temp_image_folder)
            else:
                shutil.rmtree(self._temp_image_folder)
                os.mkdir(self._temp_image_folder)

    def remove_files_only(self, directory):
        for root, dirs, files in os.walk(directory):
            for file in files:
                file_path = os.path.join(root, file)
                os.remove(file_path)
            break

    def save_final_dataset(self):
        if len(self.df_list) > 0:
            self.df = pd.concat(self.df_list, axis=0, ignore_index=True).reset_index(
                drop=True
            )
        self._save_dataset(f"{self.results_path}/final_evaluation_data.csv")
        if os.path.isfile(self.settings.data.path):
            shutil.copy2(
                self.settings.data.path,
                f"{self.results_path}/final_scaffold_memory.csv",
            )

    def save_substructure_dict(self):
        self.track_substructure.save_substruct(f"{self.results_path}/substructure.json")

    def update_current_mol(self, direction):
        self.currentMolIndex = (self.currentMolIndex + direction) % self.df.shape[0]

    def check_if_evaluated(self):
        currentProperties = self.df.iloc[
            self.currentMolIndex, -(self.numProperties + 1) : -1
        ]

        if self.unrated_sample.equals(currentProperties) == False:
            self.df.at[self.currentMolIndex, "evaluated"] = 1

    def get_smiles_and_ratings(self):
        return self.df.get_rated_smiles()

    @property
    def reached_max_iterations(self):
        return self._max_iter is not None and self.iteration == self._max_iter

    @property
    def next_unrated_mol(self):
        indexNotEvaluated = self.df[self.df.evaluated == False].index.values
        stepSize = 0
        if len(indexNotEvaluated) != 0:
            differences = indexNotEvaluated - self.currentMolIndex
            if np.sum(differences > 0) > 0:
                stepSize = np.min(differences[differences > 0])
            else:
                stepSize = np.min(differences[differences < 0])
        return stepSize

    @property
    def other_text(self):
        return self.df.getOtherText(self.currentMolIndex, self.currentTab)

    @other_text.setter
    def other_text(self, text: str):
        self.df.saveOtherText(self.currentMolIndex, self.currentTab, text)

    @property
    def current_smiles(self):
        return self.df.SMILES.iloc[self.currentMolIndex]

    @property
    def atom_selection(self):
        return self.df.loc[
            self.currentMolIndex,
            f"{self._current_liability}AtomSelection_{self.currentTab}",
        ]

    @atom_selection.setter
    def atom_selection(self, atoms):
        self.df.setSelectedAtoms(
            self.currentMolIndex,
            [self._current_liability, self.currentTab],
            atoms,
            Chem.MolFromSmiles(self.current_smiles),
        )

    @property
    def rating(self):
        return self.df.getRating(self.currentMolIndex)

    @rating.setter
    def rating(self, value):
        self.df.at[self.currentMolIndex, "rating"] = value

    @property
    def concerns(self):
        return self.df.getGlobalLiab(self.currentMolIndex)

    @concerns.setter
    def concerns(self, values):
        self.df.loc[
            self.currentMolIndex,
            self.settings.ui.global_properties.liabilities,
        ] = values

    @property
    def property_labels(self):
        return self.df.getPropertyLabels(self.currentMolIndex)

    @property
    def all_molecules_evaluated(self):
        return sum(self.df.evaluated) != self._num_molecules

    @property
    def current_liability(self):
        return self._current_liability

    @current_liability.setter
    def current_liability(self, liability):
        self._current_liability = liability

    @property
    def explanation(self):
        pixmap = draw.set_image(
            self._temp_image_folder,
            "atomContribution",
            self.currentMolIndex,
            self.df.SMILES.iloc[self.currentMolIndex],
            self._initial_qsar_model,
            self.settings.ui.show_atom_contributions.ECFP,
        )
        return pixmap

    @property
    def similar_actives(self):
        pixmap = draw.set_image(
            self._temp_image_folder,
            "mostSimilarActives",
            self.currentMolIndex,
            self.df.SMILES.iloc[self.currentMolIndex],
            self._training_data_path,
        )
        return pixmap

    @property
    def alternative_molecule(self):
        return self.df.loc[self.currentMolIndex, "alternativeMolecule"]

    @alternative_molecule.setter
    def alternative_molecule(self, mol):
        try:
            smiles = Chem.MolToSmiles(mol)
        except:
            smiles = None
            print("Invalid Molecule. Cannot be saved")
        if smiles is not None:
            self.df.loc[self.currentMolIndex, "alternativeMolecule"] = smiles

    def update_models(self):
        self.df_list.append(self.df)

        smiles, y_scores = self.df.get_rated_smiles()
        if (self.settings.reward_model is not None) & (len(smiles) > 0):
            complete_mcc, new_mcc = self.RFTrainer.update_model(smiles, y_scores)

        if len(smiles) > 0:
            self.df_result.append(self.RFTrainer.df)

            if self._shouldProcessInnerLoop():
                smiles, dict_list = self._updateModelAndSend()
                return smiles, dict_list
            else:
                return None, None

    def _shouldProcessInnerLoop(self):
        return (
            self._max_inner_loop_iter is None
            or (self._max_inner_loop_iter - 1) == self.inner_iteration
        )

    def _updateModelAndSend(self):

        self.df = self.df_list[0]
        if len(self.df_list) > 1:
            [self.df.custom_append(x) for x in self.df_list[1:]]
        data.createResultsFolder(f"{self.results_path}/iteration_{self.iteration}")
        self._save_dataset(
            f"{self.results_path}/iteration_{self.iteration}/evaluation_data.csv"
        )
        self._copyScaffoldMemory()
        dict_list = None

        if self.settings.reward_model is not None:
            results_df = self._concatResultsDF()
            self._saveResultsCSV(results_df)
            dict_list = self._createDictList(results_df)
            self.RFTrainer.save_model(
                f"{PKGDIR}/reinvent_connect/input_files/current_run/Model.pkl"
            )
            self.RFTrainer.save_model(
                f"{self.results_path}/iteration_{self.iteration}/Model.pkl"
            )

        user_scoring_function = self._getUserScoringFunction()
        self.sent2Reinvent(user_scoring_function)
        return results_df.smiles, dict_list

    def _save_dataset(self, path: str):
        self.df.to_csv(path, index=False)

    def _copyScaffoldMemory(self):
        shutil.copy2(
            self.settings.data.path,
            f"{self.results_path}/iteration_{self.iteration}/scaffold_memory.csv",
        )

    def _concatResultsDF(self):
        results_df = pd.concat(self.df_result, axis=0, ignore_index=True)
        results_df["true_initial_pred"] = self.df.loc[
            self.df.SMILES.isin(results_df.smiles), self._activity_label
        ].values.tolist()

        return results_df

    def _getUserScoringFunction(self):
        if self._human_component:
            self.track_substructure.load_data(self.df)
            return self.track_substructure.create_scoring_function()
        else:
            return None

    def _saveResultsCSV(self, results_df):
        results_df.to_csv(
            f"{self.results_path}/iteration_{self.iteration}/oracle_results.csv",
            index=False,
        )

    def sent2Reinvent(
        self,
        user_scoring_function=[],
    ):
        """
        The function starts a worker thread using QThreadPool which starts the
        Reinvent run on the server.
        """

        self.worker = tr.Worker(
            self.denovo_runner,
            self._activity_label,
            self.settings.data.path,
            user_scoring_function,
        )
        self.worker.signals.finished.connect(self.init_next_iteration)
        if self.settings.de_novo_model is not None:
            self.thread = QThreadPool()
            self.thread.start(self.worker)
        else:
            self.worker.signals.finished.emit()

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

    def _update_inner_loop(self):
        print(f"Inner Loop: {self.inner_iteration}")
        self.inner_iteration += 1
        self.set_per_round_variables(initial=False)
        self.clear_temp_images(files_only=True)
