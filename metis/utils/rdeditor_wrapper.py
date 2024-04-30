from rdeditor import MainWindow


class MetisEditor(MainWindow):
    def __init__(self, **kwargs):
        super(MetisEditor, self).__init__(**kwargs)
        self.hide()
        self.saveAction.triggered.disconnect(self.saveFile)
        self.openAction.setVisible(False)
        self.saveAsAction.setVisible(False)
        self.original_mol = None

    def clearCanvas(self):
        self.editor.clearAtomSelection()
        self.editor.mol = self.original_mol
        self.fileName = None
        self.statusBar().showMessage("Original Molecule reloaded")

    def exitFile(self):
        response = self.msgApp(
            "Confirmation",
            "This will exit the editor and return to Metis. Do you want to Continue?",
        )
        if response == "Y":
            self.ptable.close()
            self.hide()
        else:
            self.editor.logger.debug("Abort closing")

    def show_and_raise(self):
        self.show()
        self.raise_()
        self.showNormal()
