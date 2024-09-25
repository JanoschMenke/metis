from dataclasses import dataclass, field
from typing import List, Dict, Optional, Any, Tuple
from metis.config.settings import BaseConfig
from PySide6.QtCore import QThreadPool, QObject, Signal
from metis.core import data
import os
import pandas as pd
import shutil
from metis import PKGDIR


@dataclass
class MoleculeData:
    smiles: str
    properties: Dict[str, float]
    evaluation: Optional[Dict[str, Any]] = None
    additional_info: Optional[Dict[str, Any]] = None
    files: Optional[Dict[str, Any]] = None


@dataclass(frozen=True)
class AppState:
    current_mol_index: int = 0
    iteration: int = 0
    inner_iteration: int = 0
    current_tab: str = "General"
    current_liability: Optional[str] = None
    evaluated_molecules: set = field(default_factory=set)
    molecule_memory: list = field(default_factory=list)


class MoleculeDataHandler:
    def __init__(self, settings: Dict):
        self.settings = settings
        self._molecules = self.load_initial_data()

    def get_molecule(self, index: int) -> MoleculeData:
        return self._molecules[index]

    def update_evaluation(self, index: int, evaluation: Dict[str, Any]) -> None:
        self._molecules[index].evaluation = evaluation

    def update_files(self, index: int, files: Dict[str, Any]) -> None:
        self._molecules[index].files = files

    def load_initial_data(self) -> List[MoleculeData]:
        """
        Load initial data and convert it to a list of MoleculeData objects.

        Returns:
            List[MoleculeData]: A list of MoleculeData objects.
        """
        df, _ = data.loadData(self.settings.dict(), initial=True)
        columns_to_use = ["SMILES"] + [
            self.settings.dict()["propertyLabels"][name]
            for name in self.settings.dict()["propertyLabels"]
        ]
        molecules = []
        for _, row in df.iterrows():
            molecule = MoleculeData(
                smiles=row["SMILES"],
                properties={col: row[col] for col in columns_to_use[1:]},
                evaluation=None,  # Initially, no evaluation
                additional_info={
                    col: row[col] for col in df.columns if col not in columns_to_use
                },
            )
            molecules.append(molecule)

        return molecules

    def get_all_data(self) -> pd.DataFrame:
        flattened_data = []
        for mol in self._molecules:
            flat_dict = {}
            flat_dict["smiles"] = mol.smiles

            # Flatten properties
            for key, value in mol.properties.items():
                flat_dict[key] = value

            # Flatten evaluation if it exists
            if mol.evaluation:
                for key, value in mol.evaluation.items():
                    if key == "multi_select":
                        for key_2, value_2 in value.items():
                            flat_dict[key_2] = value_2
                    else:
                        flat_dict[key] = value

            # Flatten additional_info if it exists
            if mol.additional_info:
                for key, value in mol.additional_info.items():
                    flat_dict[key] = value

            flattened_data.append(flat_dict)

        return pd.DataFrame(flattened_data)

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
    def evaluated_molecules(self) -> Tuple[List[str], List[float]]:
        evaluated = [mol for mol in self._molecules if mol.evaluation is not None]
        return evaluated

    @property
    def rated_molecules(self):
        eval_mols = self.evaluated_molecules
        rated_molecules = [
            mol
            for mol in eval_mols
            if ("general_rating" in mol.evaluation)
            and (mol.evaluation["general_rating"]) > -1
        ]
        return rated_molecules or []

    def __len__(self):
        return len(self._molecules)
