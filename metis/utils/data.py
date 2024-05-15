import pandas as pd
import numpy as np
from metis.utils import helper
from rdkit.Chem import AllChem as Chem
from rdkit import DataStructs
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprint
from rdkit.Chem import DataStructs, MACCSkeys
from datetime import datetime
from typing import List, Dict
from metis.reinvent_connect import train_rf as trf
import os
import shutil
from pathlib import Path
import warnings
from typing import Tuple


class evalDataFrame(pd.DataFrame):
    """A subclass of pd.DataFrame for handling evaluation data with additional attributes and settings.

    Args:
        data (dict or array-like): Data to be passed to the DataFrame constructor.
        settings (dict): Dictionary containing settings for the evaluation.
    """

    def __init__(self, data, settings: Dict):
        super().__init__(data=data)
        self.attrs["settings"] = settings
        self.attrs["globalLiab"] = settings["ui"]["global_properties"]["liabilities"]
        self.attrs["slider"] = settings["ui"]["general"]["slider"]
        self.attrs["numColsBefore"] = self.shape[0]
        self.attrs["liability_ll"] = [
            name for name in settings["ui"]["substructures"]["liabilities"]
        ] + ["other"]

        if settings["ui"]["tab"]["render"] == True:
            tab_names = settings["ui"]["tab"]["tab_names"]
        else:
            tab_names = ["General"]
        for i in self.attrs["liability_ll"]:
            for name in tab_names:
                self.loc[:, f"{i}AtomSelection_{name}"] = np.empty(
                    (len(self), 0)
                ).tolist()
                self.loc[:, f"{i}SMARTS_{name}"] = ""
        for name in tab_names:
            self.loc[:, f"otherTextEntry_{name}"] = ""

        if settings["ui"]["general"]["slider"] == True:
            self.loc[:, "rating"] = 50
        else:
            self.loc[:, "rating"] = None
        self.loc[:, self.attrs["globalLiab"]] = 0
        self.loc[:, "alternativeMolecule"] = ""
        self.loc[:, "evaluated"] = 0

    def getPropertyLabels(self, index: int):
        out_dict = {
            name: self[self.attrs["settings"]["propertyLabels"][name]].iloc[index]
            for name in self.attrs["settings"]["propertyLabels"]
        }
        return out_dict

    def getRating(self, index: int):
        return self.at[index, "rating"]

    def getGlobalLiab(self, index: int):
        return self.loc[index, self.attrs["globalLiab"]].values.tolist()

    def setSelectedAtoms(self, index, liab, selection, mol):
        self.at[index, f"{liab[0]}AtomSelection_{liab[1]}"] = selection
        self.at[index, f"{liab[0]}SMARTS_{liab[1]}"] = (
            helper.MolFragmentToSmartsEnlarged(mol, selection)
        )

    def saveOtherText(self, index: int, tab_name: str, value):
        self.at[index, f"otherTextEntry_{tab_name}"] = value

    def getOtherText(self, index: int, tab_name: str):
        return self.at[index, f"otherTextEntry_{tab_name}"]

    def get_rated_smiles(self):
        if self.attrs["slider"] == True:
            smiles = self.SMILES.values
            rating = self.rating.values / 100
        else:
            smiles = self.SMILES[self.rating.isin([0, 2])].values
            rating = np.maximum(0, self.rating[self.rating.isin([0, 2])].values - 1)
        return smiles, rating

    def get_good_smiles(self):
        if self.attrs["slider"] == True:
            smiles = self.SMILES.loc[self.rating > 66].values.tolist()
        else:
            smiles = self.SMILES.loc[self.rating == 2].values.tolist()
        return smiles

    def get_bad_smiles(self):
        if self.attrs["slider"] == True:
            smiles = self.SMILES.loc[self.rating <= 33].values.tolist()
        else:
            smiles = self.SMILES.loc[self.rating == 0].values.tolist()
        return smiles

    def get_mediocre_smiles(self):
        if self.attrs["slider"] == True:
            smiles = self.SMILES.loc[
                (self.rating > 33) & (self.rating <= 66)
            ].values.tolist()
        else:
            smiles = self.SMILES.loc[self.rating == 0].values.tolist()
        return smiles

    def custom_append(self, df: pd.DataFrame):
        self.reset_index(inplace=True)
        df.reset_index(inplace=True)
        for i in range(df.shape[0]):
            self.loc[(self.shape[0]) + i] = df.loc[i]
        self.reset_index(inplace=True, drop=True)
        self.drop(self.columns[0], axis=1, inplace=True)


def loadData(settings: Dict, initial: bool = True) -> Tuple[pd.DataFrame, int]:
    """
    The function `loadData` loads data from a CSV file based on the provided settings, applies a
    selection strategy to the data if specified, and returns the resulting DataFrame along with a subset
    of its columns.
    """

    if initial:
        df = pd.read_csv(settings["data"]["initial_path"])
    else:
        df = pd.read_csv(settings["data"]["path"])

    if "sampled" not in df.columns:
        df["sampled"] = 0

    keep_df = df[df.sampled == 0].reset_index()

    if settings["data"]["selection_strategy"] == "random":
        selector = randomActiveSelector(
            keep_df, settings["activity_label"], settings["data"]["num_molecules"]
        )
    if settings["data"]["selection_strategy"] == "qbc":
        selector = qbcActiveSelector(
            keep_df, settings["activity_label"], settings["data"]["num_molecules"]
        )
    if settings["data"]["selection_strategy"] == "greedy":
        selector = greedyActiveSelector(
            keep_df, settings["activity_label"], settings["data"]["num_molecules"]
        )
    if settings["data"]["selection_strategy"] == "margin":
        selector = marginActiveSelector(
            keep_df, settings["activity_label"], settings["data"]["num_molecules"]
        )
    if settings["data"]["selection_strategy"] == "ucb":
        selector = UCBSelector(
            keep_df, settings["activity_label"], settings["data"]["num_molecules"]
        )

    selectIdx = selector.get_selection()
    old_idx = keep_df.loc[selectIdx, "index"]

    df.loc[old_idx, "sampled"] = 1
    df.to_csv(settings["data"]["path"], index=False)

    keep_df = keep_df.iloc[selectIdx, :]
    if type(settings["data"]["selection_strategy"]) == list:
        keep_df = keep_df.iloc[settings["data"]["selection_strategy"], :]

    keep_df.reset_index(inplace=True, drop=True)
    df = evalDataFrame(keep_df, settings)

    return df, df.iloc[0, -(df.shape[1] - df.attrs["numColsBefore"]) : -1]


def createResultsFolder(directory: str, debug: bool = False) -> str:
    """
    Creates a results folder at the specified directory or returns the existing directory path.

    Args:
        directory (str): The path to the directory where the results folder will be created.
        debug (bool, optional): If True, no additional folder with a timestamp will be created
                               in case the directory already exists. Default is False.

    Returns:
        str: The path to the created or existing results folder.

    If the specified directory does not exist, a new folder is created.
    If the directory already exists and debug is False, a new folder with a timestamp is created.
    Returns the original directory path if no new folder is created.
    """

    if not os.path.exists(directory):
        os.makedirs(directory)

    elif debug == False:
        now = datetime.now()
        dt_string = now.strftime("%d%m%Y_%H:%M:%S")
        os.makedirs(f"{directory}_{dt_string}")
        return f"{directory}_{dt_string}"
    else:
        shutil.rmtree(directory)
        os.makedirs(directory)

    return directory


def createRGBColorDict(settings):
    """
    The function `createRGBColorDict` takes a dictionary of settings and returns a dictionary where the
    keys are names and the values are RGB color values.
    """

    rgbColorDict = {}
    for name in settings:
        rgbColorDict[name] = helper.hex2RGB(settings[name].color)

    rgbColorDict["other"] = (0.753, 0.753, 0.753)
    return rgbColorDict


def select_counterfactual(originalSmiles, toKeep):
    test = pd.read_csv("../data/scaffold_memory.csv")
    filter_col = [col for col in test.columns if col.startswith("to_")]
    matching = test[(test[filter_col] == 1).sum(axis=1) == len(filter_col)]
    if matching.shape[0] == 0:
        return []
    toKeep = [Chem.MolFromSmarts(smarts) for smarts in toKeep]
    toAvoid = helper.getMolecularRest(originalSmiles, toKeep)

    newSubstruct = [
        helper.getMolecularRest(smiles, toKeep) for smiles in matching.SMILES
    ]
    descriptor_df, original_values = helper.getDescriptorDataFrame(
        newSubstruct, toAvoid
    )
    smiles_dict = helper.select_molecules_to_show(
        descriptor_df, matching, original_values
    )
    return smiles_dict


# The `diverseSelector` class is used to select a diverse set of molecules based on their similarity
# using the MaxMinPicker algorithm.
class diverseSelector:
    def __init__(
        self, smiles: List[str], num_to_select: int = 20, seed: int = -1
    ) -> None:
        self.picker = MaxMinPicker()
        self.smiles = np.array(smiles)
        self.num_to_select = num_to_select
        self.seed = seed
        self.mols = [Chem.MolFromSmiles(smi) for smi in smiles]
        self.fps = [GetMorganFingerprint(m, 2) for m in self.mols]

    def distij(self, i: int, j: int, fps=None) -> float:
        """
        The function `distij` calculates the Tanimoto distance between two fingerprints stored in a list.
        """
        return 1 - DataStructs.TanimotoSimilarity(self.fps[i], self.fps[j])

    def get_selection(self) -> List[str]:
        """
        The function `get_selection` returns a list of indices selected using the `LazyPick` method from the
        `picker` object.
        """
        idx_to_select = self.picker.LazyPick(
            self.distij, len(self.mols), self.num_to_select, seed=self.seed
        )
        return list(idx_to_select)


class randomActiveSelector:
    def __init__(
        self,
        df: pd.DataFrame,
        activity_label,
        num_to_select: int = 20,
        seed: int = None,
    ) -> None:
        self.df = df
        self.num_to_select = num_to_select
        self.seed = seed
        self.activity_label = activity_label

    def get_selection(self) -> List[str]:
        """
        The function `get_selection` returns a list of indices selected using the `LazyPick` method from the
        `picker` object.
        """
        if len(self.df[self.df[self.activity_label] > 0.5]) < self.num_to_select:
            selection = self.df[self.df[self.activity_label] > 0.5]
        else:
            selection = self.df[self.df[self.activity_label] > 0.5].sample(
                n=self.num_to_select,
                replace=False,
            )
        idx_to_select = selection.index.values.tolist()
        return idx_to_select


class qbcActiveSelector:
    def __init__(
        self,
        df: pd.DataFrame,
        activity_label,
        num_to_select: int = 20,
        seed: int = None,
    ) -> None:
        self.df = df
        self.num_to_select = num_to_select
        self.seed = seed
        self.activity_label = activity_label

    def get_selection(self) -> List[str]:
        """
        The function `get_selection` returns a list of indices selected using the `LazyPick` method from the
        `picker` object.
        """
        prop_votes = self.df[self.df[self.activity_label] > 0.5][
            self.activity_label
        ].values
        vote_entropy = []
        for i in range(len(prop_votes)):
            vote_entropy.append(
                -np.sum(
                    [
                        1 - prop_votes[i] * np.log2(1 - prop_votes[i]),
                        prop_votes[i] * np.log2(prop_votes[i]),
                    ]
                )
            )
        idx_to_select = list(
            np.argsort(vote_entropy)[::-1][: self.num_to_select]
        )  # get the n highest entropies
        return idx_to_select


class greedyActiveSelector:
    def __init__(
        self,
        df: pd.DataFrame,
        activity_label: str,
        num_to_select: int = 20,
        seed: int = None,
    ) -> None:
        self.df = df
        self.num_to_select = num_to_select
        self.seed = seed
        self.activity_label = activity_label

    def get_selection(self) -> List[str]:
        """
        The function `get_selection` returns a list of indices selected using the `LazyPick` method from the
        `picker` object.
        """
        selection = (
            self.df[self.df[self.activity_label] > 0.5]
            .sort_values(by=[self.activity_label], ascending=False)
            .head(self.num_to_select)
        )
        idx_to_select = selection.index.values.tolist()
        return idx_to_select


class marginActiveSelector:
    def __init__(
        self,
        df: pd.DataFrame,
        activity_label: str,
        num_to_select: int = 20,
        seed: int = None,
    ) -> None:
        self.df = df
        self.num_to_select = num_to_select
        self.seed = seed
        self.activity_label = activity_label

    def get_selection(self) -> List[str]:
        predicted_prob = np.column_stack(
            [
                1
                - self.df[self.df[self.activity_label] > 0.5][
                    self.activity_label
                ].values,
                self.df[self.df[self.activity_label] > 0.5][self.activity_label].values,
            ]
        )
        rev = np.sort(predicted_prob, axis=1)[:, ::-1]
        values = rev[:, 0] - rev[:, 1]
        idx_to_select = list(np.argsort(values)[::-1][: self.num_to_select])
        return idx_to_select


class UCBSelector:
    def __init__(
        self,
        df: pd.DataFrame,
        activity_label: str,
        num_to_select: int = 20,
        c: float = 1.0,
        seed: int = None,
    ) -> None:
        self.df = df
        self.c = c
        self.num_to_select = num_to_select
        self.seed = seed
        self.activity_label = activity_label

    def get_selection(self) -> List[str]:
        # initialize arrays to store the total reward and the number of times each option is selected
        total_reward = np.zeros(len(self.df[self.df[self.activity_label] > 0.5]))
        num_selections = np.zeros(len(self.df[self.df[self.activity_label] > 0.5]))

        for _ in range(self.num_to_select):
            ucb_values = total_reward / (num_selections + 1e-5) + self.c * np.sqrt(
                np.log(sum(num_selections) + 1) / (num_selections + 1e-5)
            )
            selected_index = np.argmax(ucb_values)

            # simulate selecting the option and updating the total reward and selection count
            reward = self.df[self.df[self.activity_label] > 0.5].iloc[selected_index][
                self.activity_label
            ]
            total_reward[selected_index] += reward
            num_selections[selected_index] += 1

        # Return the indices of the selected options
        idx_to_select = list(np.argsort(num_selections)[::-1][-self.num_to_select :])
        return idx_to_select


def clean_smarts(smarts: str) -> str:
    """
    Cleans a SMARTS string by converting it to canonical SMARTS representation.

    Args:
        smarts (str): The input SMARTS string.

    Returns:
        str: The cleaned SMARTS string.
    """
    mol = Chem.MolFromSmarts(smarts)

    if mol is not None:
        try:
            canonical_smarts = Chem.MolToSmarts(
                Chem.MolFromSmiles(Chem.MolToSmiles(mol))
            )
        except:
            canonical_smarts = smarts
        return canonical_smarts
    else:
        # Handle the case when the input SMARTS is invalid
        raise ValueError("Invalid SMARTS string: {}".format(smarts))


def extract_liabilities_dataframes(df):
    """
    Extracts DataFrames for each SMARTS General liability from the input DataFrame.

    Args:
        df (pandas.DataFrame): Input DataFrame.

    Returns:
        List[pandas.DataFrame]: List containing DataFrames for each SMARTS General liability.
    """
    extracted_dataframes_alert = []
    extracted_dataframes_desired = []
    # Iterate through columns in the DataFrame
    for column_name in df.columns:
        # Check if the column represents a SMARTS General liability
        if column_name.endswith("SMARTS_General"):
            # Extract the liability name from the column
            liability_name = column_name.split("SMARTS_General")[0]

            # Create a temporary DataFrame with relevant columns
            relevant_columns = [
                "SMILES",
                f"{liability_name}AtomSelection_General",
                column_name,
            ]
            temp_df = df.loc[~df[column_name].isna(), relevant_columns]

            # Check if the temporary DataFrame is not empty
            if not temp_df.empty:
                if liability_name != "like":
                    extracted_dataframes_alert.append(temp_df)
                else:
                    extracted_dataframes_desired.append(temp_df)
    return [extracted_dataframes_alert, extracted_dataframes_desired]


def extract_and_process_liabilities(data: pd.DataFrame, liability_dict):
    """
    Extracts liabilities from the given data and processes them, updating the liability dictionary.

    Args:
        data (pandas.DataFrame):
            The input DataFrame containing chemical information.
        liability_dict (dict):
            A dictionary to store identified liabilities.

    Returns:
        dict:
            Updated liability dictionary after processing the data.
    """
    found_liabilities_dfs = extract_liabilities_dataframes(data)
    for i, df_list in enumerate(found_liabilities_dfs):
        for df in df_list:
            if i < 1:
                process_dataframe(df, liability_dict)
            else:
                process_dataframe(df, liability_dict, alert=False)

    return liability_dict


def process_dataframe(df: pd.DataFrame, liability_dict: Dict, alert: bool = True):
    """
    Processes a pandas DataFrame containing molecular information to identify substructures based on SMARTS patterns.

    Args:
        df (pd.DataFrame): The input DataFrame containing molecular information.
        liability_dict (Dict): A dictionary to store substructure matches and corresponding SMARTS patterns.
        alert (bool, optional): If True, alerts will be stored in the "alerts" key of the liability_dict;
                              if False, matches will be stored in the "desired" key. Default is True.

    Returns:
        None

    This function iterates through each row of the DataFrame, extracts SMARTS patterns from the last column,
    and identifies substructure matches in the molecular structure specified in the first column.
    The matched substructures are then converted to SMARTS notation and added to the liability_dict.
    """

    for i in range(len(df)):
        smarts_patterns = df.iloc[i, -1].split(".")

        for pat in smarts_patterns:
            smiles_i = df.iloc[i, 0]
            mol_i = Chem.MolFromSmiles(smiles_i)
            pattern_mol = Chem.MolFromSmarts(pat)

            substructure_matches = mol_i.GetSubstructMatches(pattern_mol)
            highlighted_atoms = set(df.iloc[i, 1])
            matches = [highlighted_atoms.intersection(x) for x in substructure_matches]

            converted_matches = [
                helper.MolFragmentToSmartsEnlarged(
                    mol_i, match, fancySmarts=False, includeEnlargedEnv=False
                )
                for match in matches
            ]

            cleaned_smarts = [clean_smarts(x) for x in converted_matches]
            if cleaned_smarts:
                liability_dict["alerts" if alert else "desired"].setdefault(
                    cleaned_smarts[0], set()
                ).add(pat)


def sample_training_data(path_training_data, current_smiles_query, sample_size=3):
    initial_data = pd.read_csv(path_training_data)
    initial_actives = initial_data.loc[initial_data["target"] == 1, :].copy()
    similarities = []
    # calculate Tanimoto similarity with query smiles
    maccs_query = MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(current_smiles_query))
    for initial_smiles in initial_actives.smiles.tolist():
        maccs_initial = MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(initial_smiles))
        similarities.append(DataStructs.TanimotoSimilarity(maccs_query, maccs_initial))

    initial_actives["query_similarity"] = similarities

    initial_actives = initial_actives.sort_values(
        by=["query_similarity"], ascending=False
    )

    smiles = [current_smiles_query] + initial_actives.head(sample_size).smiles.tolist()
    values = {0: {"Similarity": 1.0}}

    for i, val in enumerate(initial_actives.query_similarity.values):
        values[i + 1] = {"Similarity": np.round(val, 3)}

    return smiles, values
