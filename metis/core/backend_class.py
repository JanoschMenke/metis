from metis.ui import counterfactual_mol_display, tutorial
from metis.core import data
from metis.utils import (
    helper,
    tracker,
    draw,
)
from metis.config.settings import BaseConfig
from metis.models.random_forest import trainRF
from metis.models.denovo_generator import DeNovoRunner, WorkerSignals, Worker
from metis.utils.tracker import MetisTracker
from metis import PKGDIR

from PySide6.QtCore import QThreadPool, QObject, Signal
import os
from os.path import join
import pickle
import numpy as np
from rdkit import Chem
import shutil
import pandas as pd
from typing import List, Dict, Optional, Any, Tuple
from dataclasses import dataclass


@dataclass(frozen=True)
class MoleculeData:
    smiles: str
    properties: Dict[str, float]
    evaluation: Optional[Dict[str, Any]] = None
    additional_info: Optional[Dict[str, Any]] = None


@dataclass(frozen=True)
class AppState:
    current_mol_index: int = 0
    iteration: int = 0
    inner_iteration: int = 0
    current_tab: str = "General"
    current_liability: Optional[str] = None


class MoleculeDataHandler:
    def __init__(self, settings: Dict):
        self.settings = settings
        self._molecules = self.load_initial_data()

    def get_molecule(self, index: int) -> MoleculeData:
        return self._molecules[index]

    def update_evaluation(self, index: int, evaluation: Dict) -> None:
        old_mol = self._molecules[index]
        new_mol = MoleculeData(
            smiles=old_mol.smiles, properties=old_mol.properties, evaluation=evaluation
        )
        self._molecules[index] = new_mol

    def get_evaluated_molecules(self) -> Tuple[List[str], List[float]]:
        evaluated = [mol for mol in self._molecules if mol.evaluation is not None]
        return (
            [mol.smiles for mol in evaluated],
            [mol.evaluation["score"] for mol in evaluated],
        )

    def load_initial_data(self) -> List[MoleculeData]:
        """
        Load initial data and convert it to a list of MoleculeData objects.

        Returns:
            List[MoleculeData]: A list of MoleculeData objects.
        """
        df, _ = data.loadData(self.settings.dict(), initial=True)

        molecules = []
        for _, row in df.iterrows():
            molecule = MoleculeData(
                smiles=row["SMILES"],
                properties={
                    name: row[self.settings.dict()["propertyLabels"][name]]
                    for name in self.settings.dict()["propertyLabels"]
                },
                evaluation=None,  # Initially, no evaluation
            )
            molecules.append(molecule)

        return molecules

    def get_all_data(self) -> pd.DataFrame:
        return pd.DataFrame([mol.__dict__ for mol in self._molecules])

    def clear_temp_images(self, files_only: bool = False) -> None:

        temp_image_folder = f"{PKGDIR}/resources/temp_images/"

        if os.path.exists(temp_image_folder):
            if files_only:
                self.remove_files_only(temp_image_folder)
            else:
                shutil.rmtree(temp_image_folder)
                os.mkdir(temp_image_folder)

    def remove_files_only(self, directory):
        for root, dirs, files in os.walk(directory):
            for file in files:
                file_path = os.path.join(root, file)
                os.remove(file_path)
            break

    @property
    def num_molecules(self) -> int:
        return len(self._molecules)


class BackendSignals(QObject):
    state_changed = Signal(AppState)
    molecule_evaluated = Signal(int)  # Emits index of evaluated molecule


class Backend(QObject):
    def __init__(self, settings: BaseConfig, results_folder: str):
        super().__init__()
        self.settings = settings
        self.load_settings(self.settings)
        self.signals = BackendSignals()
        self.results_folder = results_folder
        self.data_handler = MoleculeDataHandler(settings)
        self._state = AppState()
        self._setup_initial_state()

        self.init_models()
        self.tracker = MetisTracker()

        self.set_variables(results_folder=results_folder, initial=True)
        self.set_global_variables()
        self.data_handler.clear_temp_images()

    def _setup_initial_state(self):
        """Initialize the application state."""
        self.color_dict = self.color_dict = data.createRGBColorDict(
            self.settings.ui.substructures.liabilities
        )
        self._state = AppState()
        self.signals.state_changed.emit(self._state)

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

    def init_models(self):
        """Initialize RF and De Novo models."""
        self.rf_model = (
            trainRF(self.settings.reward_model)
            if self.settings.reward_model is not None
            else None
        )
        self.denovo_model = (
            DeNovoRunner(self.settings.de_novo_model)
            if self.settings.de_novo_model is not None
            else None
        )

    def set_variables(self, results_folder: str = "", initial: bool = True):
        if initial:
            self.results_path = data.createResultsFolder(
                join(results_folder, self.settings.data.run_name),
                debug=self._debug,
            )

        self.set_per_round_variables(initial=initial)
        self.data_handler.clear_temp_images(files_only=True)

    def set_per_round_variables(self, initial: bool = True):
        np.random.seed(self.settings.seed)
        self.df, self.unrated_sample = data.loadData(
            self.settings.dict(), initial=initial
        )
        new_state = AppState(
            current_mol_index=0,
            iteration=self._state.iteration,
            inner_iteration=self._state.inner_iteration,
        )
        self._update_state(new_state)

    def update_current_mol(self, direction: int) -> None:
        """Update the current molecule index."""
        new_index = (
            self._state.current_mol_index + direction
        ) % self.data_handler.num_molecules
        self._update_state(current_mol_index=new_index)

    def _update_state(self, **kwargs) -> None:
        """Update the application state."""
        new_state = AppState(**{**self._state.__dict__, **kwargs})
        self._state = new_state
        self.signals.state_changed.emit(new_state)

    def set_global_variables(self):
        """Set up global variables and objects."""
        self.next_iteration_signal = WorkerSignals()
        if self._tutorial:
            self.tutorial = tutorial.Second(self)
        self.counterfactual_window = counterfactual_mol_display.CounterfactualWindow()

    def init_next_iteration(self):
        if self.settings.de_novo_model is not None:
            shutil.copy2(
                f"{PKGDIR}/resources/input_files/current_run/Agent.ckpt",
                f"{self.results_path}/iteration_{self.iteration}/Agent.ckpt",
            )
        self.inner_iteration = 0
        self.setVariables(initial=False)
        self.iteration += 1
        self.next_iteration_signal.finished.emit()

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
        self.tracker.save_substruct(f"{self.results_path}/substructure.json")

    def check_if_evaluated(self):
        """Check if the current molecule has been evaluated."""
        current_mol = self.data_handler.get_molecule(self._state.current_mol_index)
        if current_mol.evaluation is not None:
            self.signals.molecule_evaluated.emit(self._state.current_mol_index)

    def get_smiles_and_ratings(self):
        """Get SMILES and ratings for evaluated molecules."""
        return self.data_handler.get_evaluated_molecules()

    @property
    def reached_max_iterations(self):
        return self._max_iter is not None and self._state.iteration == self._max_iter

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
                f"{PKGDIR}/resources/input_files/current_run/Model.pkl"
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

        self.worker = Worker(
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
