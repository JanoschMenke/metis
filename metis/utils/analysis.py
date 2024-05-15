# %%
import pandas as pd
import numpy as np
from rdkit.Chem import AllChem as Chem
from rdkit.Chem.Scaffolds import MurckoScaffold as ms
from typing import List, Dict, Tuple


class analyser:
    def __init__(self, path_initial_data: str, path_all_scaffolds: str) -> None:
        data = pd.read_csv(path_initial_data)
        self.logging = list()
        self.training_scaffolds = data.iloc[:, 0].unique().tolist()
        self.all_scaffolds = (
            pd.read_csv(path_all_scaffolds, header=None).iloc[:, 0].values.tolist()
        )
        self.found_scaffolds_evaluated = set(
            self.training_scaffolds
        )  # tracks scaffolds evaluated
        self.found_scaffolds_generated = set(
            self.training_scaffolds
        )  # tracks scaffolds in scaffold_memory.csv

    def newTopologicalScaffolds(self, smiles: List[str], label: str = "") -> Dict:
        scaffolds = set([self.genericMurckoScaffold(smi) for smi in smiles])
        out_dict = {
            f"{label}_scaffolds_in_training": scaffolds.intersection(
                self.training_scaffolds
            )
        }
        out_dict[f"{label}_scaffolds_already_found"] = scaffolds.intersection(
            self.found_scaffolds_evaluated
        )
        out_dict[f"{label}_novel_scaffolds"] = (
            scaffolds - out_dict[f"{label}_scaffolds_already_found"]
        ).intersection(self.all_scaffolds)

        return out_dict

    def trackFoundScaffolds(
        self, evaluatedDataframe: pd.DataFrame, generatedDataframe: pd.DataFrame
    ) -> Tuple[Dict, Dict]:
        evaluatedOutDict = self.newTopologicalScaffolds(
            evaluatedDataframe.SMILES, label="evaluated"
        )
        self.found_scaffolds_evaluated.update(
            evaluatedOutDict["evaluated_novel_scaffolds"]
        )

        generatedOutDict = self.newTopologicalScaffolds(
            generatedDataframe.SMILES, label="generated"
        )
        self.found_scaffolds_generated.update(
            generatedOutDict["generated_novel_scaffolds"]
        )
        self.logging.append(
            [evaluatedOutDict[key] for key in evaluatedOutDict]
            + [generatedOutDict[key] for key in generatedOutDict]
        )
        evaluatedOutDict.update(generatedOutDict)
        valueDict = {name: len(evaluatedOutDict[name]) for name in generatedOutDict}

        return evaluatedOutDict, valueDict

    @staticmethod
    def genericMurckoScaffold(smi: str) -> str:
        scaff = ms.GetScaffoldForMol(Chem.MolFromSmiles(smi))
        scaff = ms.MakeScaffoldGeneric(scaff)
        return Chem.MolToSmiles(scaff)


# %%
