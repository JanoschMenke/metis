from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdqueries
from rdkit.Chem import rdDepictor, rdMolDescriptors, rdFingerprintGenerator
from rdkit.DataStructs import TanimotoSimilarity
from rdkit.Chem.Draw import rdMolDraw2D
from copy import deepcopy
import numpy as np
import random
import string
import numpy as np
import pandas as pd
from typing import List, Dict
from pathlib import Path
import os
import PySide6


def get_random_string(length: int):
    # choose from all lowercase letter
    letters = string.ascii_lowercase
    result_str = "".join(random.choice(letters) for i in range(length))
    return result_str


def clear_current_files(path: str):
    # Convert path to a Path object
    path_obj = Path(path)

    # Check if the specified directory exists, if not, create it
    if not path_obj.exists():
        path_obj.mkdir(parents=True)
    try:
        # Iterate over files in the directory and delete them
        [f.unlink() for f in path_obj.glob("*") if f.is_file()]
    except Exception as e:
        print(f"An error occurred while clearing files: {e}")


def read_and_replace_css(css_path, pkg_path):
    with open(css_path, "r") as file:
        css_content = file.read()
        # Replace placeholder with actual package directory
        css_content = css_content.replace("__PACKAGE_DIR__", pkg_path)
    return css_content


def MolFragmentToSmartsEnlarged(
    mol, atomids: List[int], fancySmarts: bool = True, includeEnlargedEnv: bool = True
):
    """An extended Version of RDKits MolFragmentToSmarts(), can include the outgoing bonds and ring information

    Args:
        mol (rdkit mol objects): An rdkit molecule
        atom ids (list): Indices of Atoms part of the substructure to convert
        includeRingInfo (bool, optional): Can include the information whether the Atom is part of a Ring. Defaults to True.
        includeEnlargedEnv (bool, optional): Includes information on the bonds connect to atoms not part of the substructure. Defaults to True.

    Returns:
        smarts (str): A SMARTS String
    """

    # this is done to ensure that two distinct substructures do not become a
    # single substructe when extending the env
    trueEnv, enlargedEnv = getSubPaths(mol, atomids)
    intermediate_smarts = path2Smarts(mol, list(trueEnv))
    if "." in intermediate_smarts:
        actual_atom_envs = []
        smarts_patterns = intermediate_smarts.split(".")
        for sma in smarts_patterns:
            pattern = Chem.MolFromSmarts(sma)
            matches = mol.GetSubstructMatches(pattern)
            out = [list(x) for x in matches if set(list(x)).issubset(atomids)]
            actual_atom_envs.append(out[0])
    else:
        actual_atom_envs = [atomids]

    if fancySmarts:
        mol = fancySMARTS(Chem.RWMol(mol))
    smarts_list = []
    for atomids in actual_atom_envs:
        trueEnv, enlargedEnv = getSubPaths(mol, atomids)

        if includeEnlargedEnv == True:
            enlargedEnv = list(enlargedEnv)
        else:
            enlargedEnv = list(trueEnv)

        smarts_list.append(path2Smarts(mol, enlargedEnv))

    return ".".join(smarts_list)


def getSubPaths(mol, atomids):
    atomids = set(atomids)
    trueEnv = set()
    enlargedEnv = set()
    for atom in atomids:
        a = mol.GetAtomWithIdx(atom)
        for b in a.GetBonds():
            if b.GetIdx() in enlargedEnv:
                trueEnv.add(b.GetIdx())
            enlargedEnv.add(b.GetIdx())
    return trueEnv, enlargedEnv


def path2Smarts(mol, env: List[int]):
    """
    Converts a specified path in a molecule to a SMARTS pattern.

    Args:
        mol (Chem.Mol): The input molecule.
        env (List[int]): List of atom indices specifying the path in the molecule.

    Returns:
        str: SMARTS pattern corresponding to the specified path in the molecule.

    """
    amap = {}
    submol = Chem.PathToSubmol(mol, env, atomMap=amap)
    smarts = Chem.MolToSmarts(submol)
    return smarts


def hex2RGB(hex: str):
    "Convert Hex to RGB, we devide by 255 to make it readable by rdkit"
    hex = hex.lstrip("#")
    rgb = tuple(int(hex[i : i + 2], 16) for i in (0, 2, 4))
    return tuple(i / 255 for i in rgb)


def genScaffold(smiles: str, selection):
    mol = Chem.MolFromSmiles(smiles)
    atomsToKeep = set(range(mol.GetNumAtoms())) - set(selection)
    bondsConnected = []
    exitBonds = []
    for atomid in atomsToKeep:
        atom = mol.GetAtomWithIdx(atomid)
        for bond in atom.GetBonds():
            if getTrueEndAtom(bond, atom) in atomsToKeep:
                bondsConnected.append(bond.GetIdx())
            else:
                exitBonds.append(bond.GetIdx())
                mol.GetAtomWithIdx(getTrueEndAtom(bond, atom)).SetAtomicNum(0)
    bondsConnected = list(set(bondsConnected))
    new_mol = Chem.PathToSubmol(mol, bondsConnected + exitBonds, atomMap={})
    new_smiles = prep_smiles_libevent(Chem.MolToSmiles(new_mol))
    return new_smiles


def prep_smiles_libevent(smiles):
    """
    The function `prep_smiles_libevent` replaces each occurrence of "*" in a given string with "[*:x]",
    where x is a counter that increments for each occurrence of "*". This enables the Smiles to be used
    for libinvent
    """
    new_smiles = ""
    counter = 0
    for charakter in smiles:
        if charakter == "*":
            new_smiles += f"[*:{counter}]"
            counter += 1
        else:
            new_smiles += charakter
    return new_smiles


def fancySMARTS(mol, exitAtomIdx=[None]):
    """
    The function `fancySMARTS` takes a molecule and replaces each atom with a SMARTS query that includes
    the atom's atomic number, number of hydrogens, and whether it is in a ring or not.
    """
    for i in range(mol.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(i)
        numHs = atom.GetTotalNumHs()
        atomicNum = atom.GetAtomicNum()
        inRing = atom.IsInRing()

        queries = rdqueries.AtomNumEqualsQueryAtom(atomicNum)
        if i not in exitAtomIdx:
            a = Chem.AtomFromSmarts(f"[H{numHs}]")
            queries.ExpandQuery(a)
        b = Chem.AtomFromSmarts(f'[{"R" if inRing else "R0"}]')
        queries.ExpandQuery(b)
        mol.ReplaceAtom(i, queries)
    return mol


def findExitAtomIdx(mol, atomids):
    "identifies the atom ids of atom connect to the part that is to be removed"
    atomsToKeep = list(set(list(range(mol.GetNumAtoms()))) - set(atomids))

    atomsConnected = []
    for atomid in atomsToKeep:
        atom = mol.GetAtomWithIdx(atomid)
        for bond in atom.GetBonds():
            if getTrueEndAtom(bond, atom) in atomids:
                atomsConnected.append(atom.GetIdx())
    return atomsConnected


def getTrueEndAtom(bond, atom):
    "Returns the idx of the atom connected to the atom given by bond given"
    if not (bond.GetIdx() in [b.GetIdx() for b in atom.GetBonds()]):
        raise Exception(
            f"Bond {bond.GetIdx()} is not connected to Atom {atom.GetIdx()}"
        )

    if atom.GetIdx() == bond.GetEndAtomIdx():
        return bond.GetBeginAtomIdx()
    return bond.GetEndAtomIdx()


def getBondConnectingAtoms(mol, idx1, idx2):
    connectingBondIndex = None
    atom = mol.GetAtomWithIdx(idx1)
    for bond in atom.GetBonds():
        endAtom = getTrueEndAtom(bond, atom)
        if endAtom == idx2:
            connectingBondIndex = bond.GetIdx()
    return connectingBondIndex


def increase_resolution(mol, substructure, size=(400, 200)):
    mol = deepcopy(mol)
    substructure = deepcopy(substructure)
    drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])

    # highlightAtoms expects only one tuple, not tuple of tuples. So it needs to be merged into a single tuple
    matches = sum(mol.GetSubstructMatches(substructure), ())
    drawer.DrawMolecule(mol, highlightAtoms=matches)

    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()

    return svg.replace("svg:", "")


def calcSimilarityToOriginal(smiles1, smiles2):
    fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2)
    m1 = fpgen.GetFingerprint(Chem.MolFromSmiles(smiles1))
    m2 = fpgen.GetFingerprint(Chem.MolFromSmiles(smiles2))
    return TanimotoSimilarity(m1, m2)


def is_faulty_pyside_version():
    return False
    """
    "shitty hack to check the verson of Pyside"
    version = PySide2.__version__
    version = version.split(".")
    version = int("".join(version[:3]))
    if (version >= 5151) & (version <= 5152):
        return True
    else:

        return False
    """


def create_color_dictionary(settings):
    """
    The function `createRGBColorDict` takes a dictionary of settings and returns a dictionary where the
    keys are names and the values are RGB color values.
    """

    rgbColorDict = {}
    for name in settings:
        rgbColorDict[name] = hex2RGB(settings[name].color)

    rgbColorDict["other"] = (0.753, 0.753, 0.753)
    return rgbColorDict
