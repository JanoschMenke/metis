import pickle
import pandas as pd
import numpy as np
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdFingerprintGenerator
from sklearn.metrics import f1_score
from typing import List, Dict

# from torch import tensor


class trainRF:
    def __init__(self, reward_model_settings):
        self.settings = reward_model_settings

        self.oracle_available = False
        if self.settings.oracle_path is not None:
            self.oracle_available = True
            self.load_oracle(self.settings.oracle_path)
        self.load_model(self.settings.qsar_model_path)

        self.ecfpGenerator = ecfp_generator(
            bitSize=self.settings.ECFP.bitSize,
            radius=self.settings.ECFP.radius,
            useCounts=self.settings.ECFP.useCounts,
        )
        self.training_data = self.create_data_dict(self.settings.training_data_path)

    def create_data_dict(self, path: str) -> dict:
        """
        Creates a data dictionary containing input features (x) and corresponding labels (y)
        from a CSV file located at the specified path.

        Args:
            path (str): The file path to the CSV file containing the data.

        Returns:
            dict: A dictionary containing 'x' for input features and 'y' for labels.

        """

        data_dict = {}
        data = pd.read_csv(path)
        data_dict["y"] = data.target.values
        fps = self.ecfpGenerator.get_fingerprints(data.smiles.values)
        data_dict["x"] = fps
        if "weight" not in data.columns:
            data_dict["weight"] = np.array([1] * len(data_dict["x"]))
        else:
            data_dict["weight"] = data.weight.values

        return data_dict

    def update_model(self, smiles: List[str], y: List = None):
        """
        Updates the model using the provided SMILES representations and corresponding labels.

        Args:
            smiles (List[str]): A list of SMILES representations for chemical compounds.
            y (List, optional): A list of labels corresponding to the input SMILES. If None,
                the labels will be retrieved from the model's oracle (if available). Default is None.


        Returns:
            None: The function updates the internal model and does not return any value.

        Note:
            If the oracle is available, it retrieves the labels for the given SMILES using the oracle.
            If the oracle is not available, the function expects y to be provided explicitly.
            The model is then trained using the input SMILES and labels, and performance metrics
            such as complete_f1 and new_f1 are calculated.
        """

        fps, y, weights = self.get_model_inputs(smiles, y)
        self.train_rf(fps, y, weights)

        complete_f1, new_f1 = self.calc_f1(fps, y)

        print(f"Trainingsdata Size: {self.training_data['x'].shape[0]}")
        return complete_f1, new_f1

    def save_model(self, path: str) -> None:
        pickle.dump(self.model, open(path, "wb"))

    def calc_f1(self, x: np.ndarray, y: np.ndarray):
        """
        The function `calc_f1` calculates the F1 score for a given model's predictions on training data and
        new data.
        """

        y_hat = self.model.predict(self.training_data["x"])
        complete_f1 = f1_score(self.training_data["y"], y_hat)
        y_hat = self.model.predict(x)
        self.df["prediction_by_updated_model"] = y_hat
        new_f1 = f1_score(y, y_hat)

        return complete_f1, new_f1

    def load_oracle(self, oracle_path: str) -> None:
        self.oracle = pickle.load(open(oracle_path, "rb"))

    def load_model(self, qsar_model_path: str):
        """
        The `load_model` function loads a model from a given file path and enables warm start.
        Args:
            qsar_model_path (str): path to RF model that should be updated
        """
        self.model = pickle.load(open(qsar_model_path, "rb"))

        self.model.warm_start = False
        # self.model.n_estimators += 100
        self.model.class_weight = None

    def create_trackDF(self, size: int) -> pd.DataFrame:
        df = pd.DataFrame(np.zeros([size, 4]))
        df.columns = [
            "smiles",
            "prediction_by_current_model",
            "oracle_pred",
            "prediction_by_updated_model",
        ]
        return df

    def get_model_inputs(self, smiles: List[str], y: List = None):
        """
        Generate inputs for the model.

        Args:
            smiles: A list of SMILES strings.
            y: A list of target values. If None, oracle predictions will be used.

        Returns:
            A tuple containing:
                - fps: An array of fingerprints.
                - targets: An array of target values.
                - weights: An array of weights.
        """
        if self.settings.use_oracle_score:
            smiles = np.array(smiles)[np.array(y) == 1]
            y = np.array(y)[np.array(y) == 1]

        self.df = self.create_trackDF(len(smiles))
        fps = self.ecfpGenerator.get_fingerprints(smiles)
        y = self.oracle.predict(fps)

        self.df["smiles"] = smiles
        self.df["target"] = np.round(y).astype(int)
        self.df["prediction_by_current_model"] = self.model.predict(fps)

        if self.oracle_available:
            self.df["oracle_pred"] = self.oracle.predict(fps)

        weights = self._calculate_weights(y)

        return fps, self.df["target"].values, weights

    def train_rf(self, fps, y, weights):
        """
        The `train_rf` function appends new data to the existing training data and fits the model with the
        updated training data.

        """
        self.training_data["x"] = np.concatenate((self.training_data["x"], fps), axis=0)
        self.training_data["y"] = np.concatenate((self.training_data["y"], y), axis=0)
        self.training_data["weight"] = np.concatenate(
            (self.training_data["weight"], weights), axis=0
        )
        self.model.fit(
            self.training_data["x"],
            self.training_data["y"],
            sample_weight=self.training_data["weight"],
        )

    def get_prob_distribution(self, smiles):
        fps = self.ecfpGenerator.get_fingerprints(smiles)
        prob_dist = [
            estimator.predict_proba(fps) for estimator in self.model.estimators_
        ]
        prob_dist = np.stack(prob_dist, axis=1)
        # prob_dist = tensor(prob_dist)
        return prob_dist

    def _calculate_weights(self, target_values):
        if self.settings.weight == "pseudo_confidence":
            weights = 0.5 + np.abs(0.5 - np.array(target_values))
        elif self.settings.weight is None:
            weights = np.array([1] * len(target_values))
        return weights

    @property
    def RFModel(self):
        return self.model


# The `ecfp_generator` class is a Python class that generates ECFP fingerprints for a given list of
# SMILES strings.
class ecfp_generator:
    def __init__(self, bitSize: int = 2048, radius: int = 2, useCounts: bool = False):
        self.bitSize = bitSize
        self.radius = radius
        self.fpgen = rdFingerprintGenerator.GetMorganGenerator(
            radius=radius, fpSize=bitSize
        )
        self.generate = (
            self.fpgen.GetCountFingerprintAsNumPy
            if useCounts
            else self.fpgen.GetFingerprintAsNumPy
        )

    def get_fingerprints(self, smiles: List[str]) -> np.ndarray:
        fps = np.stack([self.generate(Chem.MolFromSmiles(smile)) for smile in smiles])
        return fps
