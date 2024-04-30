from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdqueries
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDepictor, rdMolDescriptors, rdFingerprintGenerator
from rdkit.DataStructs import TanimotoSimilarity
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG
from copy import deepcopy
import numpy as np
import random
import string
import numpy as np
import pandas as pd
from typing import List, Dict


def get_random_string(length: int):
    # choose from all lowercase letter
    letters = string.ascii_lowercase
    result_str = "".join(random.choice(letters) for i in range(length))
    return result_str


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


def splitMolecule(mol, atomids, smarts=False):
    """_summary_
    Removes the Selected Atoms from  molecule
    and provides the left over parts as SMILES

    Args:
        mol (rdkit mol): molecule
        atomids (list[int]): ids of atoms to remove

    Returns:
        list[str]: list of SMILES
    """
    bondsOnlyBetweenSelectedAtoms = set()
    bondsSelection = set()

    for atom in atomids:
        a = mol.GetAtomWithIdx(atom)
        for b in a.GetBonds():
            if getTrueEndAtom(b, a) in atomids:
                bondsOnlyBetweenSelectedAtoms.add(b.GetIdx())
            bondsSelection.add(b.GetIdx())
    leftoverBonds = list(set(range(mol.GetNumBonds())) - set(bondsSelection))

    mol = fixAromaticity(mol, atomids, bondsSelection=bondsSelection)
    if smarts == False:
        leftoverSmiles = Chem.MolToSmiles(Chem.PathToSubmol(mol, leftoverBonds))
        to_avoid_smiles = Chem.MolToSmiles(
            Chem.PathToSubmol(mol, list(bondsOnlyBetweenSelectedAtoms))
        )
    else:
        leftoverSmiles = Chem.MolToSmarts(Chem.PathToSubmol(mol, leftoverBonds))
        to_avoid_smiles = Chem.MolToSmarts(
            Chem.PathToSubmol(mol, list(bondsOnlyBetweenSelectedAtoms))
        )
    return leftoverSmiles.split("."), to_avoid_smiles


def fixAromaticity(mol, atomids, bondsSelection):
    exitAtoms = findExitAtomIdx(mol, atomids)
    if len(exitAtoms) == 2:
        bondBetween = getBondConnectingAtoms(mol, exitAtoms[0], exitAtoms[1])
        if bondBetween is not None:
            if mol.GetBondWithIdx(bondBetween).GetIsAromatic():
                ringInfo = mol.GetRingInfo()
                ringBonds = ringInfo.BondRings()
                ringsWithoutSelection = [
                    ring
                    for ring in ringBonds
                    if len(set(bondsSelection).intersection(ring)) == 0
                ]
                ringsWithSelection = [
                    ring for ring in ringsWithoutSelection if bondBetween in ring
                ]
                aromaticRing = all(
                    [
                        mol.GetBondWithIdx(ridx).GetIsAromatic()
                        for ridx in ringsWithSelection[0]
                    ]
                )
                if not aromaticRing:
                    bond = mol.GetBondWithIdx(bondBetween)
                    bond.SetBondType(Chem.BondType.SINGLE)
                    bond.SetIsAromatic(False)
                    mol.GetAtomWithIdx(bond.GetEndAtomIdx()).SetIsAromatic(False)
                    mol.GetAtomWithIdx(bond.GetBeginAtomIdx()).SetIsAromatic(False)
    return mol


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
        raise Exception(f"Bond {b.GetIdx()} is not connected to Atom {atom.GetIdx()}")

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


def fancySubstruct(smiles, atomids, fancy=True, smarts=True):
    """
    The function `fancySubstruct` generates the molecules does not contain the selected atoms
    which are given by atomids.

    Args:
        smiles: str, smiles to transform
        atomids: list[int], list of the selected atoms that should be removed
        fancy: bool, if True it creates a SMARTS/SMILES that includes more information
                     on connectivity and ring membership than a regular SMARTS
        smarts: bool, if True it outputs a SMARTS string else a SMILES


    """

    keepSmarts, _ = splitMolecule(Chem.MolFromSmiles(smiles), atomids, smarts=True)
    mol = Chem.RWMol(Chem.MolFromSmiles(smiles))
    if fancy == True:
        if len(keepSmarts) == 1:
            exitAtomIdx = findExitAtomIdx(mol, atomids)
            mol = fancySMARTS(mol, exitAtomIdx)
        else:
            mol = fancySMARTS(mol)

    return splitMolecule(mol, atomids, smarts=smarts)


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


def getDescriptorDataFrame(newSubstruct, toAvoid):
    descriptors = [
        calcSubstructDescriptors(smiles, toAvoid)
        for i, smiles in enumerate(newSubstruct)
    ]

    descriptors = pd.DataFrame(descriptors)
    descriptors.columns = [
        "hetero_atom_ratio",
        "rot_bonds_ratio",
        "num_atoms",
        "num_rings",
        "sim_to_original",
    ]

    original_values = calcSubstructDescriptors(toAvoid, toAvoid)
    return descriptors, original_values


def calcSubstructDescriptors(smiles, toAvoid):
    mol = Chem.MolFromSmiles(smiles)
    heteroRatio = getHeteroAtomRatio(mol)
    rotatableBondRatio = getRotatableBondRatio(mol)
    numAtoms = mol.GetNumAtoms()
    numRings = rdMolDescriptors.CalcNumRings(mol)
    sim2Original = calcSimilarityToOriginal(smiles, toAvoid)
    return [heteroRatio, rotatableBondRatio, numAtoms, numRings, sim2Original]


def getMolecularRest(smiles, toKeepSmarts):
    "Gets the SMARTS for the newly generated part of the molecule"
    mol = Chem.MolFromSmiles(smiles)
    matched_atoms = [list(mol.GetSubstructMatch(x)) for x in toKeepSmarts]
    matched_atoms = sum(matched_atoms, [])
    out = splitMolecule(mol, matched_atoms)
    return out[0][0]


def getHeteroAtomRatio(mol):
    atomTypes = np.array([atom.GetAtomicNum() for atom in mol.GetAtoms()])
    return np.mean(atomTypes != 6)


def getRotatableBondRatio(mol):
    numRotatableBonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if mol.GetNumBonds() == 0:
        return 0
    return numRotatableBonds / mol.GetNumBonds()


def calcSimilarityToOriginal(smiles1, smiles2):
    fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2)
    m1 = fpgen.GetFingerprint(Chem.MolFromSmiles(smiles1))
    m2 = fpgen.GetFingerprint(Chem.MolFromSmiles(smiles2))
    return TanimotoSimilarity(m1, m2)


def select_molecules_to_show(descriptor_df, matching, original_values):
    molecules_selected = {}

    index_to_select = descriptor_df[
        (descriptor_df.sim_to_original > 0) & (descriptor_df.sim_to_original < 1)
    ].sim_to_original.idxmax()
    molecules_selected["max_sim"] = matching.SMILES.iloc[index_to_select]
    descriptor_df.drop(index_to_select, inplace=True)
    index_to_select = descriptor_df[
        (descriptor_df.sim_to_original > 0)
        & (descriptor_df.sim_to_original < 1)
        & (descriptor_df.num_atoms < int(2.5 * original_values[2]))
    ].sim_to_original.idxmin()
    molecules_selected["min_sim"] = matching.SMILES.iloc[index_to_select]
    descriptor_df.drop(index_to_select, inplace=True)

    # select molecules with more rings
    if original_values[3] == 0:
        # from generated molecules select molecules with rings
        molecules_with_rings = descriptor_df[descriptor_df.num_rings > 0]
        # if such molecules exists
        if molecules_with_rings.shape[0] != 0:
            # check if one only contains a single ring
            if sum(molecules_with_rings.num_rings == 1) > 0:
                # find molecule with one ring that has as many atoms as the original
                index_to_select = select_ring_molecule(
                    descriptor_df, original_values[2], num_rings=1
                )
                molecules_selected["ring"] = matching.SMILES.iloc[index_to_select]
            else:
                index_to_select = select_ring_molecule(
                    descriptor_df, original_values[2], num_rings=None
                )
                molecules_selected["ring"] = matching.SMILES.iloc[index_to_select]

    # if selection has a ring select one without a ring
    else:
        index_to_select = select_ring_molecule(
            descriptor_df, original_values[2], num_rings=original_values[3] - 1
        )
        molecules_selected["ring"] = matching.SMILES.iloc[index_to_select]
    descriptor_df.drop(index_to_select, inplace=True)
    # SELECT HETERO ALTERNATIVE
    if original_values[0] <= np.nanmedian(descriptor_df.hetero_atom_ratio):
        min_ratio = np.percentile(descriptor_df.hetero_atom_ratio.dropna(), 80)
        max_ratio = np.percentile(descriptor_df.hetero_atom_ratio.dropna(), 90)

    else:
        min_ratio = np.percentile(descriptor_df.hetero_atom_ratio.dropna(), 10)
        max_ratio = np.percentile(descriptor_df.hetero_atom_ratio.dropna(), 20)
    try:
        index_to_select = (
            descriptor_df[
                (descriptor_df.hetero_atom_ratio >= min_ratio)
                & (descriptor_df.hetero_atom_ratio <= max_ratio)
            ]
            .sample(1)
            .index.tolist()[0]
        )
        molecules_selected["hetero"] = matching.SMILES.iloc[index_to_select]
    except:
        None

    return molecules_selected


def select_ring_molecule(df, num_atoms_original, num_rings=None):
    if num_rings is None:
        comparison = df.num_rings > 0
    else:
        comparison = df.num_rings == num_rings

    min_diff_atom_num = min(abs(df[comparison].num_atoms - num_atoms_original))
    index_to_select = (
        df[
            (comparison)
            & ((abs(df.num_atoms - num_atoms_original) == min_diff_atom_num))
        ]
        .sample(1)
        .index.tolist()[0]
    )
    return index_to_select
