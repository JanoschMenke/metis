from typing import List, Dict, Optional, Any, Tuple
from pathlib import Path
from PySide6.QtCore import QObject, Signal, Slot
from PySide6.QtGui import QPixmap
from rdkit import Chem
import pickle
import numpy as np

from metis.core import data
from metis.core.storing import AppState, MoleculeDataHandler
from metis.config.settings import BaseConfig
from metis.utils import helper, draw


class BackendSignals(QObject):
    state_changed = Signal(AppState)
    molecule_updated = Signal(Chem.Mol)
    evaluation_updated = Signal(object)
    get_atoms = Signal(str)
    atom_selection_updated = Signal(list)
    text_in_other = Signal(str, str)
    properties = Signal(dict)
    similarity_image = Signal(QPixmap)
    explanation_image = Signal(QPixmap)


class Backend(QObject):
    def __init__(self, settings: BaseConfig, results_folder: Path):
        super().__init__()
        self.settings = settings
        self.molecule_handler = MoleculeDataHandler(settings)
        self.state = AppState(current_liability=self._initial_liability)
        self.signals = BackendSignals()
        if self.settings.ui.show_atom_contributions.render:
            self._initial_qsar_model = pickle.load(
                open(settings.ui.show_atom_contributions.path, "rb")
            )
        self.set_additional_properties(results_folder=results_folder)

    def set_additional_properties(self, results_folder: Path):
        self.results_path = data.createResultsFolder(
            results_folder / self.settings.data.run_name,
            debug=self.settings.debug,
        )

    @property
    def _initial_liability(self) -> str:
        return next(iter(self.settings.ui.substructures.liabilities))

    def load_initial_molecule(self) -> None:
        if len(self.molecule_handler) > 0:
            self.state = AppState(
                current_mol_index=0,
                current_liability=self._initial_liability,
                current_tab="General",
            )
            self._update_current_molecule()
        else:
            self._emit_empty_state()

    def _emit_empty_state(self):
        self.signals.state_changed.emit(self.state)
        self.signals.molecule_updated.emit(None)
        self.signals.evaluation_updated.emit({})

    @Slot(int)
    def handle_navigation(self, direction: int) -> None:
        self.signals.get_atoms.emit(self.state.current_liability)
        new_index = (self.state.current_mol_index + direction) % len(
            self.molecule_handler
        )
        self._update_state(current_mol_index=new_index)
        self._update_current_molecule()

    def _update_state(self, **kwargs):
        self.state = AppState(**{**self.state.__dict__, **kwargs})

    def _update_current_molecule(self) -> None:
        current_mol_data = self.molecule_handler.get_molecule(
            self.state.current_mol_index
        )
        mol = Chem.MolFromSmiles(current_mol_data.smiles)

        self.signals.state_changed.emit(self.state)
        self.signals.molecule_updated.emit(mol)
        self.signals.evaluation_updated.emit(current_mol_data.evaluation or {})
        self._retrieve_atom_selection()
        self._retrieve_text_in_other()
        self._retrieve_properties()
        self.explanation() if self.settings.ui.show_atom_contributions.render else None
        self.similar_actives() if self.settings.ui.show_reference_molecules else None
        self.gen_mol_image()

    @Slot(dict)
    def update_evaluation(self, evaluation: Dict[str, Any]):
        self.molecule_handler.update_evaluation(
            self.state.current_mol_index, evaluation
        )
        self.signals.evaluation_updated.emit(evaluation)

    @Slot(dict)
    def update_files(self, files: Dict[str, Any]):
        self.molecule_handler.update_files(self.state.current_mol_index, files)

    @Slot(int)
    def handle_general_rating_change(self, value: int):
        self._update_evaluation_field("general_rating", value)
        self._update_state(
            evaluated_molecules=self.state.evaluated_molecules
            | {self.state.current_mol_index}
        )

    @Slot(dict)
    def handle_multi_select_change(self, values: Dict[str, int]):
        self._update_evaluation_field("multi_select", values)

    def _update_evaluation_field(self, field: str, value: Any):
        current_mol_data = self.molecule_handler.get_molecule(
            self.state.current_mol_index
        )
        evaluation = current_mol_data.evaluation or {}
        evaluation[field] = value
        self.update_evaluation(evaluation)

    @Slot(str)
    def handle_liability_change(self, liability: str) -> None:
        self.signals.get_atoms.emit(self.state.current_liability)
        self._update_state(current_liability=liability)
        self.signals.state_changed.emit(self.state)
        self._retrieve_atom_selection()

    @Slot(str)
    def handle_tab_change(self, tab: str) -> None:
        self.signals.get_atoms.emit(self.state.current_liability)
        self._update_state(current_tab=tab)
        self.signals.state_changed.emit(self.state)
        self._retrieve_atom_selection()
        self._retrieve_text_in_other()

    def _retrieve_text_in_other(self) -> None:
        current_mol_data = self.molecule_handler.get_molecule(
            self.state.current_mol_index
        )
        evaluation = current_mol_data.evaluation or {}
        key = f"TextInOther_{self.state.current_tab}"
        text_in_other = evaluation.get(key, "Other?")
        if text_in_other:
            self.signals.text_in_other.emit(self.state.current_tab, text_in_other)

    def _retrieve_atom_selection(self) -> None:
        current_mol_data = self.molecule_handler.get_molecule(
            self.state.current_mol_index
        )
        evaluation = current_mol_data.evaluation or {}
        key = f"{self.state.current_liability}AtomSelection_{self.state.current_tab}"
        atoms = evaluation.get(key, [])
        self.signals.atom_selection_updated.emit(atoms)

    def store_atom_selection(
        self, atoms: List[int], text_in_other: Optional[str] = None
    ) -> None:
        current_mol_data = self.molecule_handler.get_molecule(
            self.state.current_mol_index
        )
        evaluation = current_mol_data.evaluation or {}
        key = f"{self.state.current_liability}AtomSelection_{self.state.current_tab}"
        evaluation[key] = atoms

        # SMARTS Pattern
        key = f"{self.state.current_liability}SMARTS_{self.state.current_tab}"
        evaluation[key] = self._get_smarts_pattern(current_mol_data, atoms)

        if (
            text_in_other
            and self.state.current_liability == "other"
            and text_in_other != "Other?"
        ):
            evaluation[f"TextInOther_{self.state.current_tab}"] = text_in_other

        self.update_evaluation(evaluation)

    def _retrieve_properties(self) -> None:
        current_mol_data = self.molecule_handler.get_molecule(
            self.state.current_mol_index
        )
        properties = current_mol_data.properties or {}
        self.signals.properties.emit(properties)

    def _get_smarts_pattern(self, molecular_data, atom_selection: List[int]) -> str:
        """Convert the selected atoms to a smarts pattern"""
        mol = Chem.MolFromSmiles(molecular_data.smiles)
        smarts_pattern = helper.MolFragmentToSmartsEnlarged(mol, atom_selection)
        return smarts_pattern

    def save_molecular_data(self, path: Optional[str] = None) -> None:
        df = self.molecule_handler.get_all_data()
        if path:
            df.to_csv(path)
        else:
            df.to_csv(self.results_path / "final_evaluation.csv")

    def explanation(self):
        current_mol_data = self.molecule_handler.get_molecule(
            self.state.current_mol_index
        )
        files = current_mol_data.files or {}
        if "atomContribution" in files:
            pixmap = files["atomContribution"]
        else:
            pixmap = draw.set_image(
                "atomContribution",
                self.state.current_mol_index,
                current_mol_data.smiles,
                self._initial_qsar_model,
                self.settings.ui.show_atom_contributions.ECFP,
            )
            files["atomContribution"] = pixmap
            self.update_files(files)
        self.signals.explanation_image.emit(pixmap)

    def similar_actives(self):
        current_mol_data = self.molecule_handler.get_molecule(
            self.state.current_mol_index
        )
        files = current_mol_data.files or {}
        if "mostSimilarActives" in files:
            pixmap = files["mostSimilarActives"]
        else:
            pixmap = draw.set_image(
                "mostSimilarActives",
                self.state.current_mol_index,
                current_mol_data.smiles,
                self.settings.ui.show_reference_molecules.path,
            )
            files["mostSimilarActives"] = pixmap
            self.update_files(files)
        self.signals.similarity_image.emit(pixmap)

    @property
    def num_molecules_not_evaluated(self):
        return len(self.molecule_handler) - len(self.state.evaluated_molecules)

    def gen_mol_image(self):
        current_mol_data = self.molecule_handler.get_molecule(
            self.state.current_mol_index
        )
        files = current_mol_data.files
        if "mol_image" not in files:
            print("Generated Mol Image")
            files["mol_image"] = draw.create_mol_image(current_mol_data.smiles)
            self.update_files(files)

    @property
    def mol_wall_molecules(self):
        molecules = self.molecule_handler.rated_molecules
        print(molecules)
        mol_list = []
        for mol in molecules:
            print(mol.files)
            d = {
                "image": mol.files["mol_image"],
                "mol_info": {"Rating": mol.evaluation["general_rating"]},
            }
            mol_list.append(d)
        return mol_list

    def prepare_training_data(self):
        self.save_current_molecules()
        smiles = []
        ratings = []
        for mol in self.molecule_handler.rated_molecules:
            if self.settings.ui.general.slider == False:
                if mol.evaluation.get("general_rating") in [0, 2]:
                    smiles.append(mol.smiles)
                    ratings.append(np.maximum(0, mol.evaluation["general_rating"] - 1))
            else:
                if mol.evaluation.get("general_rating") > -1:
                    smiles.append(mol.smiles)
                    ratings.append(mol.evaluation["general_rating"] / 100)

        return smiles, ratings

    def save_current_molecules(self):
        mol_memory = self.state.molecule_memory
        mol_memory += self.molecule_handler._molecules
        self._update_state(molecule_memory=mol_memory)
