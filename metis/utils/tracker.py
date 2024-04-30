import pandas as pd
import numpy as np
from rdkit.Chem import AllChem as Chem
from metis.utils import helper, data
import json
import copy
from typing import List, Dict


class metis_tracker:
    def __init__(self):
        self.substructure_dict = {"desired": {}, "alerts": {}}
        self.smiles_dict = {"desired": [], "alerts": [], "somewhat": []}
        self.upper_limit_num_atoms = 50  # original REINVENT Limits
        self.lower_limit_num_atoms = 10  # original REINVENT Limits

    def load_data(self, dataframe: pd.DataFrame) -> None:
        data.extract_and_process_liabilities(dataframe, self.substructure_dict)
        self.process_molecular_size(dataframe)
        self.get_smiles(dataframe)

    def save_substruct(self, path: str):
        out_dict = {}
        for direction in ["desired", "alerts"]:
            out_dict[direction] = {
                name: list(self.substructure_dict[direction][name])
                for name in self.substructure_dict[direction]
            }
        with open(path, "w") as outfile:
            json.dump(out_dict, outfile)

    def get_smiles(self, dataframe: pd.DataFrame):
        self.smiles_dict["desired"] += dataframe.get_good_smiles()
        self.smiles_dict["alerts"] += dataframe.get_bad_smiles()
        self.smiles_dict["somewhat"] += dataframe.get_mediocre_smiles()

    def create_component_alerts(self):
        alerts = self.substructure_dict.get("alerts", [])

        if alerts:
            smiles_list = sum([list(smiles) for name, smiles in alerts.items()], [])

            component_custom_alerts = {
                "component_type": "custom_alerts",
                "name": "alerts",
                "weight": 1,
                "specific_parameters": {"smiles": smiles_list},
            }

            return copy.deepcopy(component_custom_alerts)

        return None

    def process_molecular_size(self, df: pd.DataFrame) -> None:
        size_liked_mol = [
            Chem.MolFromSmiles(x).GetNumAtoms() for x in df.get_good_smiles()
        ]
        too_big_smiles = df[df["Too Big"] > 0].SMILES.values.tolist()
        too_small_smiles = df[df["Too Small"] > 0].SMILES.values.tolist()
        if too_big_smiles:
            smallest_disliked_mol = np.min(
                [Chem.MolFromSmiles(x).GetNumAtoms() for x in too_big_smiles]
            )
            if self.upper_limit_num_atoms < smallest_disliked_mol:
                smallest_disliked_mol = self.upper_limit_num_atoms
            self.upper_limit_num_atoms -= (
                self.upper_limit_num_atoms - smallest_disliked_mol
            ) / 2

        if max(size_liked_mol) > self.upper_limit_num_atoms:
            self.upper_limit_num_atoms -= (
                self.upper_limit_num_atoms - max(size_liked_mol)
            ) / 2

        if too_small_smiles:
            largest_disliked_mol = np.max(
                [Chem.MolFromSmiles(x).GetNumAtoms() for x in too_small_smiles]
            )
            if self.lower_limit_num_atoms > largest_disliked_mol:
                largest_disliked_mol = self.lower_limit_num_atoms
            self.lower_limit_num_atoms += (
                largest_disliked_mol - self.lower_limit_num_atoms
            ) / 2
        if min(size_liked_mol) < self.lower_limit_num_atoms:
            self.lower_limit_num_atoms += (
                min(size_liked_mol) - self.lower_limit_num_atoms
            ) / 2

    def create_scoring_function(self) -> List:
        scoring_components = []
        scoring_components.append(self.create_component_alerts())
        scoring_components.append(self.create_component_desired())
        scoring_components.append(self.create_component_heavy_atoms())
        scoring_components += self.create_component_similarity()
        return scoring_components

    def create_component_heavy_atoms(self) -> Dict:
        component = {
            "component_type": "num_heavy_atoms",
            "name": "Heavy Atom Count",
            "weight": 1,
            "specific_parameters": {
                "transformation": {
                    "transformation_type": "double_sigmoid",
                    "high": self.upper_limit_num_atoms,
                    "low": self.lower_limit_num_atoms,
                    "coef_div": 1200,
                    "coef_si": 150,
                    "coef_se": 150,
                }
            },
        }
        return component

    def create_component_desired(self) -> Dict:
        desired = self.substructure_dict.get("desired", [])

        if desired:
            smiles_list = sum([list(smiles) for name, smiles in desired.items()], [])

            component_custom_desired = {
                "component_type": "matching_substructure",
                "name": "desired",
                "specific_parameters": {"smiles": smiles_list},
                "weight": 1,
            }

            return copy.deepcopy(component_custom_desired)

        return None

    def create_component_similarity(self) -> List:
        component_list = []
        weights = [1, -1, 0.5]
        transformations = [
            {
                "high": 0.5,
                "transformation": True,
                "transformation_type": "relu_min",
            },
            {
                "high": 0.5,
                "transformation": True,
                "transformation_type": "relu_min",
            },
            {
                "high": 0.5,
                "transformation": True,
                "transformation_type": "relu_min",
            },
        ]

        for i, direction in enumerate(self.smiles_dict):
            smiles_list = self.smiles_dict[direction]
            if smiles_list:
                component = {
                    "component_type": "tanimoto_similarity",
                    "name": f"similarity_{direction}",
                    "specific_parameters": {
                        "smiles": smiles_list,
                        "radius": 3,
                        "use_features": False,
                        "count": False,
                        "transformation": transformations[i],
                    },
                    "weight": weights[i],
                }
                component_list.append(component)

        return component_list
