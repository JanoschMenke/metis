#!/usr/bin/env python
from __future__ import print_function
from typing import List, Dict, Optional, Any
import os
from sys import platform
from pathlib import Path

# Import required modules
import sys, time, os, shutil
from PySide6.QtGui import *
from PySide6.QtWidgets import *
from PySide6.QtWidgets import QApplication
from PySide6 import QtCore, QtWidgets
from PySide6.QtCore import Slot

import argparse
from metis.config.settings import load_settings, BaseConfig
from metis.core import Controller
from metis.utils import helper

# from . import PKGDIR

PKGDIR = os.path.dirname(os.path.abspath(__file__))


# The main window class
class Metis(QtWidgets.QMainWindow):
    # Constructor function
    def __init__(
        self,
        settings: BaseConfig,
        output_path: Path,
        fileName=None,
        loglevel="WARNING",
    ):
        super(Metis, self).__init__()
        self.settings = settings

        self.controller = Controller(settings, output_path=output_path)
        self.setSize()  # set size of main window
        self.loglevels = []

        self.initGUI(fileName=fileName)
        self.load_initial_molecule()

    def setSize(self):
        self.setMinimumSize(1500, 800)
        self.setMaximumSize(1920, 1080)

    def initGUI(self, fileName=None):
        self.setCentralLayout()
        self.setWindowTitle("Metis")
        self.show()
        self.setFocus()

    def setCentralLayout(self):
        self.center = QtWidgets.QWidget()
        self.center.setLayout(self.controller.ui.setGlobalLayout())
        self.setCentralWidget(self.center)

    def load_initial_molecule(self):
        self.controller.backend.load_initial_molecule()

    def closeEvent(self, event):
        self.controller.ui.editor.logger.debug("closeEvent triggered")
        self.controller.exit_file()
        exit(0)
        event.ignore()


def launch(loglevel: str = "WARNING"):
    """
    Launch the Tetis application.

    Args:
        loglevel (str): Logging level.
    """
    if sys.platform == "darwin" and helper.is_faulty_pyside_version():
        os.environ["QT_MAC_WANTS_LAYER"] = "1"

    parser = argparse.ArgumentParser("Metis")
    parser.add_argument("-f", "--file", default="setting.yml")
    parser.add_argument("-o", "--output", default="../results/")

    args = parser.parse_args()

    if not os.path.isfile(args.file):
        parser.print_help()
        sys.exit(1)

    app = QApplication(sys.argv)
    settings: BaseConfig = load_settings(args.file)
    main_window = Metis(
        settings,
        output_path=Path(args.output),
        fileName=args.file,
        loglevel=loglevel,
    )

    stylesheet = helper.read_and_replace_css(
        f"{PKGDIR}/resources/design/style.css", PKGDIR
    )
    app.setStyleSheet(stylesheet)
    app.exec()

    sys.exit(0)


if __name__ == "__main__":
    launch(loglevel="WARNING")
