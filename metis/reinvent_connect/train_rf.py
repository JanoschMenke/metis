import pickle
import pandas as pd
import numpy as np
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdFingerprintGenerator
from sklearn.metrics import f1_score
from typing import List, Dict

# from torch import tensor
from utils.epig import epig_from_probs


class trainRF:
    def __init__(
        self,
        qsar_model_path: str,
        training_data_path: str,
        ECFP: Dict,
        oracle_path: str = None,
        oracle_score: bool = False,
        weight: str = None,
    ):
        self.weight = weight
        self.oracle_available = False
        self.use_oracle_for_prediction = oracle_score
        if oracle_path is not None:
            self.oracle_available = True
            self.load_oracle(oracle_path)
        self.load_model(qsar_model_path)

        self.ecfpGenerator = ecfp_generator(
            bitSize=ECFP["bitSize"], radius=ECFP["radius"], useCounts=ECFP["useCounts"]
        )
        self.training_data = self.create_data_dict(training_data_path)

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

        Raises:
            AssertionError: If the oracle is not available and y is not provided.

        Returns:
            None: The function updates the internal model and does not return any value.

        Note:
            If the oracle is available, it retrieves the labels for the given SMILES using the oracle.
            If the oracle is not available, the function expects y to be provided explicitly.
            The model is then trained using the input SMILES and labels, and performance metrics
            such as complete_f1 and new_f1 are calculated.
        """
        if self.use_oracle_for_prediction == True:
            assert self.oracle_available == True
            fps, y, weights = self.get_model_inputs(smiles)
        else:
            assert y is not None
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

    def get_model_inputs(self, smiles, y: List = None):
        """
        Retrieves model inputs, including fingerprints (fps) and labels (y), based on the provided SMILES.

        Args:
            smiles (List): A list of SMILES representations for chemical compounds.
            y (List, optional): A list of labels corresponding to the input SMILES. If None,
                the labels will be predicted using the model's oracle (if available). Default is None.

        Returns:
            Tuple[np.ndarray, List]: A tuple containing fingerprints (fps) and labels (y).

        Raises:
            AssertionError: If y is provided when the oracle is available.

        Note:
            - If the oracle is available, it predicts labels using the model's oracle and sets the
            'y' values accordingly.
            - If 'y' is provided explicitly, it uses the provided labels.
            - The function also creates a tracking DataFrame with information about the predictions.
        """
        self.df = self.create_trackDF(len(smiles))
        fps = self.ecfpGenerator.get_fingerprints(smiles)
        if self.use_oracle_for_prediction:
            assert y is None
            assert self.oracle_available == True
            y = self.oracle.predict(fps)
        self.df["smiles"] = smiles
        self.df["target"] = np.round(y).astype(int)
        self.df["prediction_by_current_model"] = self.model.predict(fps)

        if self.oracle_available == True:
            self.df["oracle_pred"] = self.oracle.predict(fps)
        if self.weight == "pseudo_confidence":
            weights = 0.5 + np.abs(0.5 - np.array(y))

        elif self.weight is None:
            weights = np.array([1] * len(y))

        return fps, self.df["target"], weights

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
        if len(smiles) < 2:
            fps = np.array(fps).reshape(1,-1)
        prob_dist = [
            estimator.predict_proba(fps) for estimator in self.model.estimators_
        ]
        prob_dist = np.stack(prob_dist, axis=1)
        # prob_dist = tensor(prob_dist)
        return prob_dist

    def estimate_epig(self, prob_pool, prob_target):
        return epig_from_probs(prob_pool, prob_target)

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

    def get_fingerprints(self, smiles: List[str]) -> np.ndarray[np.int32]:
        """
        The function `get_fingerprints` takes a list of SMILES strings as input and returns a numpy array of
        fingerprints.

        Args:
            smiles (List[str]): A list of SMILES strings representing chemical compounds

        Returns:
            np.ndarray[np.int32]: an array of fingerprints, where each fingerprint is represented as an array of integers.
        """
        if len(smiles) > 1:
            fps = np.stack([self.generate(Chem.MolFromSmiles(smile)) for smile in smiles])
        else:
            print("WARNING: Too few high-scoring molecules.")
            fps = np.array([self.generate(Chem.MolFromSmiles(smile)) for smile in smiles])
        return fps


# %%
