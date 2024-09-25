from rdkit.Chem import Draw
from rdkit.Chem.Draw import SimilarityMaps
from rdkit.Chem import AllChem as Chem
import numpy as np
import io
from PIL import Image
import pickle
from PySide6.QtCore import QBuffer, QIODevice, Qt
from PySide6.QtGui import QPixmap, QImage
from cairosvg import svg2png
import os
from typing import Dict, List
from functools import partial
from metis.core.data import sample_training_data

from PySide6.QtCore import QObject, Signal, Slot, QRunnable
from PySide6 import QtCore
from PySide6.QtGui import QPixmap
from metis import PKGDIR


class DrawWorkerSignals(QObject):
    """
    Defines the signals available from a running worker thread.

    Supported signals are:

    finished
        No data
    """

    finished = Signal()


def show_png(data):
    bio = io.BytesIO(data)
    img = Image.open(bio)
    return img


def get_pred(fp, pred_function):
    fp = np.array([list(fp)])
    return pred_function(fp)[0][1]


def plot_explanation_map(mol, model, ecfp_settings: Dict):
    fpType = "count" if ecfp_settings.useCounts else "bv"
    fpType = "bv"  # currently only wokr "bv"
    fpfunc = partial(
        SimilarityMaps.GetMorganFingerprint,
        nBits=ecfp_settings.bitSize,
        radius=ecfp_settings.radius,
        fpType=fpType,
    )

    d = Draw.MolDraw2DSVG(600, 600)
    SimilarityMaps.GetSimilarityMapForModel(
        mol,
        fpfunc,
        lambda x: get_pred(x, model.predict_proba),
        draw2d=d,
    )

    d.FinishDrawing()
    svg = d.GetDrawingText()
    return svg  # .replace("svg:", "")


def save_explanation_map(
    model_path,
    smiles,
    ecfp_settings: Dict = {"radius": 2, "bitSize": 2048, "useCounts": False},
):
    if type(model_path) == str:
        svg = plot_explanation_map(
            Chem.MolFromSmiles(smiles),
            pickle.load(open(model_path, "rb")),
            ecfp_settings=ecfp_settings,
        )
    else:
        svg = plot_explanation_map(
            Chem.MolFromSmiles(smiles), model_path, ecfp_settings=ecfp_settings
        )

    png_data = svg2png(
        bytestring=svg,
    )
    image = QImage()
    image.loadFromData(png_data, "PNG")
    return image


def save_most_similar_active_map(trainings_data_path: str, current_smiles: str):
    mol_list = []
    smiles, molInfo = sample_training_data(
        path_training_data=trainings_data_path,  # ,
        current_smiles_query=current_smiles,  # self.df.SMILES.iloc[self.currentMolIndex],
    )

    legends = [f"Similarity: {molInfo[i]['Similarity']}" for i in molInfo]
    legends[0] = "Generated Molecule"
    legends.insert(0, "")
    legends.insert(2, "")
    smiles.insert(0, "")
    smiles.insert(2, "")
    for smi in smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            mol_list.append(mol)

    # Create a grid image from the list of molecules
    grid_img = Draw.MolsToGridImage(
        mol_list,
        molsPerRow=3,
        subImgSize=(200, 200),
        legends=legends,
        returnPNG=False,
    )

    return convert_to_qimage(grid_img)


def convert_to_qimage(pil_image) -> QImage:
    img_buffer = io.BytesIO()
    pil_image.save(img_buffer, format="PNG")
    img_data = img_buffer.getvalue()

    q_buffer = QBuffer()
    q_buffer.open(QIODevice.ReadWrite)
    q_buffer.write(img_data)
    q_buffer.seek(0)

    qimage = QImage()
    qimage.loadFromData(q_buffer.data(), "PNG")
    return qimage


def create_mol_image(smiles: str):
    img = Draw.MolToImage(Chem.MolFromSmiles(smiles), size=(240, 240), canvas=None)
    qimg = convert_to_qimage(img)
    return QPixmap.fromImage(qimg)


def set_image(
    image_type: str,
    index: int,
    smiles,
    data_path,
    ecfp_settings=None,
):
    if image_type == "mostSimilarActives":

        image = save_most_similar_active_map(
            data_path,
            smiles,
        )

    elif image_type == "atomContribution":
        image = save_explanation_map(
            data_path,
            smiles,
            ecfp_settings=ecfp_settings,
        )
    else:
        raise ValueError(f"Invalid image_type: {image_type}")

    pixmap = QPixmap.fromImage(image)
    pixmap = pixmap.scaled(
        600,
        600,
        QtCore.Qt.KeepAspectRatio,
        mode=QtCore.Qt.SmoothTransformation,
    )
    return pixmap
