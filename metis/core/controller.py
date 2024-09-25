from metis.core.backend import Backend
from metis.core.storing import AppState
from metis.ui import gui
from metis.models.denovo_generator import DeNovoRunner
from metis.models.random_forest import trainRF
from PySide6 import QtCore
from typing import Any, Dict
from pathlib import Path


class Controller(QtCore.QObject):
    """
    A class representing the main viewer for the Metis application.

    This class handles the connection between the UI and the backend,
    managing user interactions and state updates.

    Attributes:
        atom_saver (QtCore.Signal): Signal for saving atom selections.
        backend (Backend): The backend instance for data processing.
        settings (Any): Application settings.
        ui (gui.MetisGUI): The main GUI instance.
        current_text (str): The current text in the 'other' field.
    """

    atom_saver = QtCore.Signal(list, str)

    def __init__(self, settings: Any, output_path: Path):
        """
        Initialize the Viewer instance.

        Args:
            settings (Any): Application settings.
            output_path (str): Path for output files.
        """
        super().__init__()
        self.backend: Backend = Backend(settings, output_path)
        self.settings: Any = settings
        self.ui: gui.MetisGUI = gui.MetisGUI(self.settings)
        self.current_text: str = ""
        self.connect_signals()
        self.init_models()

    def connect_signals(self) -> None:
        """Connect all signals between UI elements and backend methods."""
        self._connect_navigation_signals()
        self._connect_backend_signals()
        self._connect_evaluation_signals()
        self._connect_atom_selection_signals()

    def _connect_navigation_signals(self) -> None:
        self.ui.navigation.directionSignal.connect(self.backend.handle_navigation)
        self.ui.navigation.openWindowSignal.connect(self._nav_bar_actions)

    def _connect_backend_signals(self) -> None:
        self.backend.signals.state_changed.connect(self.update_ui_state)
        self.backend.signals.molecule_updated.connect(self.ui.update_molecule)
        self.backend.signals.evaluation_updated.connect(self.update_evaluation)
        self.backend.signals.get_atoms.connect(self.get_selected_atom)
        self.backend.signals.atom_selection_updated.connect(
            self.ui.update_atom_selection
        )
        self.backend.signals.text_in_other.connect(self.set_text_in_other)
        self.backend.signals.properties.connect(
            self.ui.evaluationWindow.TPPDispay.updateMolProperties
        )
        self.backend.signals.explanation_image.connect(self.ui.load_explanation)
        self.backend.signals.similarity_image.connect(self.ui.load_similarity_map)

    def _connect_evaluation_signals(self) -> None:
        eval_widget = self.ui.evaluationWindow.evaluationWidget
        global_interaction = eval_widget.globalInteractionWindow

        global_interaction.sliderModule.valueChanged.connect(
            self.backend.handle_general_rating_change
        )
        global_interaction.listwidget.valueChanged.connect(
            self.backend.handle_multi_select_change
        )
        eval_widget.local_liability_changed.connect(
            self.backend.handle_liability_change
        )
        eval_widget.text_in_other_changed.connect(self.update_current_text_in_other)

        if self.ui.evaluationWindow.uses_tab:
            eval_widget.tab_changed.connect(self.backend.handle_tab_change)

    def _connect_atom_selection_signals(self) -> None:
        self.atom_saver.connect(self.backend.store_atom_selection)

    def init_models(self):
        """Initialize RF and De Novo models."""
        self.rf_model = (
            trainRF(self.settings.reward_model)
            if self.settings.reward_model is not None
            else None
        )
        self.denovo_model = (
            DeNovoRunner(self.settings.de_novo_model)
            if self.settings.de_novo_model is not None
            else None
        )

    def get_selected_atom(self, current_liability: str) -> None:
        """
        Get the currently selected atoms and emit the atom_saver signal.

        Args:
            current_liability (str): The current liability type.
        """
        text_in_other = self.current_text if current_liability == "other" else None
        self.atom_saver.emit(self.ui.editor.selectedAtoms, text_in_other)

    def update_current_text_in_other(self, text: str) -> None:
        """
        Update the current text in the 'other' field.

        Args:
            text (str): The new text to set.
        """
        self.current_text = text

    def set_text_in_other(self, current_tab: str, text: str) -> None:
        """
        Set the text in the 'other' field and update the UI.

        Args:
            current_tab (str): The current tab name.
            text (str): The text to set in the 'other' field.
        """
        self.current_text = text
        self.ui.update_text_in_other(current_tab, text)

    def exit_file(self):

        response = self.ui.exit_msg(self.backend.num_molecules_not_evaluated)
        if response == "Y":
            self.ui.navigation.directionSignal.emit(0)
            df = self.backend.molecule_handler.get_all_data()
            print(df)
            self.backend.save_molecular_data()
            exit(0)

    @QtCore.Slot(str)
    def _nav_bar_actions(self, action: str):
        self.ui.navigation.directionSignal.emit(0)
        if action == "view":
            self.open_mol_viewer()

        elif action == "unrated":
            #!self.findNextUnratedMol()
            pass
        elif action == "edit":
            pass  #!self.openMolEditor()
        elif action == "send":
            self.train_and_update_model()
        elif action == "finish":
            self.exit_file()

    @QtCore.Slot(AppState)
    def update_ui_state(self, state: AppState) -> None:
        """
        Update the UI state based on the current AppState.

        Args:
            state (AppState): The current application state.
        """
        self.ui.molCounter.updateCounter(state.current_mol_index)
        self.ui.navigation.backButton.setEnabled(True)
        self.ui.navigation.nextButton.setEnabled(True)
        # Update other UI elements based on the new state

    @QtCore.Slot(object)
    def update_evaluation(self, evaluation: Dict[str, Any]) -> None:
        """
        Update the evaluation widgets with the stored evaluation data.

        Args:
            evaluation (Dict[str, Any]): The evaluation data for the current molecule.
        """
        general_rating = evaluation.get("general_rating")
        multi_select = evaluation.get("multi_select")

        global_interaction = (
            self.ui.evaluationWindow.evaluationWidget.globalInteractionWindow
        )
        global_interaction.sliderModule.eval_score = general_rating
        global_interaction.listwidget.eval_score = multi_select

    def open_mol_viewer(self):

        self.ui.molWall.initGUI(
            self.backend.state.iteration, self.backend.mol_wall_molecules
        )

    def train_and_update_model(self):
        response = self.ui.train_confirmation()
        if response == "Y":
            self.ui.disable_inputs()
            self.backend.prepare_training_data()
            smiles, ratings = self.backend.prepare_training_data()
            if (self.settings.reward_model is not None) & (len(smiles) > 0):
                complete_mcc, new_mcc = self.RFTrainer.update_model(smiles, ratings)
