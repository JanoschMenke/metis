# %%
import subprocess as sb
import copy
import json
import os
from metis.utils.helper import get_random_string
from metis.core.data import extract_and_process_liabilities
import time
from PySide6.QtCore import QObject, Signal, Slot, QRunnable
import yaml
import copy
from metis import PKGDIR


class WorkerSignals(QObject):
    """
    Defines the signals available from a running worker thread.

    Supported signals are:

    finished
        No data
    """

    finished = Signal()


class Worker(QRunnable):
    def __init__(
        self,
        de_novo_runner=None,
        activity_label=None,
        output_path=None,
        human_component=None,
    ):
        super(Worker, self).__init__()
        # self.smiles = smiles
        # self.selectedAtoms = selectedAtoms
        self.de_novo_runner = de_novo_runner
        self.activity_label = activity_label
        self.dict = human_component
        self.output_path = output_path
        self.signals = WorkerSignals()
        print("Worker Initialization")

    @Slot()
    def run(self):
        print("Run on Core")
        self.de_novo_runner.run(self.activity_label, self.output_path, self.dict)
        self.signals.finished.emit()


class DeNovoRunner:
    def __init__(self, de_novo_config):
        self.cwd = f"{PKGDIR}/resources"
        self.settings = de_novo_config
        self.ssh_settings = yaml.safe_load(open(self.settings.ssh_settings))
        self.slurm_path = self.ssh_settings["path_remote_folder"]
        self.run_name = get_random_string(8)
        with open(self.ssh_settings["de_novo_json"]) as og_file:
            self.og_json_contents = json.loads(og_file.read())

        self.agent_path = None
        self.model_path = None

    def gen_denovo_file(self, activity_label: str, human_component=None):
        file_contents = copy.deepcopy(self.og_json_contents)
        if os.path.isfile(f"{self.cwd}/input_files/current_run/Agent.ckpt"):
            file_contents["parameters"]["reinforcement_learning"][
                "agent"
            ] = f"{self.slurm_path}/Agent.ckpt"
            self.agent_path = f"{self.cwd}/input_files/current_run/Agent.ckpt"
        else:
            file_contents["parameters"]["reinforcement_learning"][
                "agent"
            ] = f"{self.slurm_path}/Agent_Initial.ckpt"

        component_list = list()
        for component in file_contents["parameters"]["scoring_function"]["parameters"]:
            if self.settings.use_reward_model:
                if component["name"] == activity_label:
                    if os.path.isfile(f"{self.cwd}/input_files/current_run/Model.pkl"):
                        component["specific_parameters"][
                            "model_path"
                        ] = f"{self.slurm_path}/Model.pkl"
                        self.model_path = (
                            f"{self.cwd}/input_files/current_run/Model.pkl"
                        )
                    else:
                        component["specific_parameters"][
                            "model_path"
                        ] = f"{self.slurm_path}/Model_Initial.pkl"
            component_list.append(component)

        if human_component is not None:
            num_specified_components = len(component_list)
            num_human_components = len(human_component)
            scale_ratio = num_specified_components / num_human_components
            for comp in human_component:
                if comp is not None:
                    comp["weight"] *= scale_ratio
                    component_list.append(comp)

        file_contents["parameters"]["scoring_function"]["parameters"] = component_list
        with open(f"{self.cwd}/input_files/current_run/new_run.json", "w") as outfile:
            json.dump(file_contents, outfile, sort_keys=True, indent=4)

    def run(self, activity_label: str, output_path: str, human_component=None):
        self.gen_denovo_file(activity_label, human_component=human_component)

        generate_slurm_file(
            f"{self.cwd}/input_files/current_run/new_run.slurm",
            self.ssh_settings["default_slurm"],
            self.slurm_path,
            run_name=self.run_name,
        )

        start_remote_run(
            self.ssh_settings["ssh_login"],
            f"{self.cwd}/input_files/current_run/new_run.slurm",
            f"{self.cwd}/input_files/current_run/new_run.json",
            self.slurm_path,
            self.model_path,
            self.agent_path,
        )

        print("Slurmjob Submitted")
        check_reinvent_status(
            self.ssh_settings["ssh_login"], self.slurm_path, start=True
        )
        print("Training has Started")
        check_reinvent_status(
            self.ssh_settings["ssh_login"], self.slurm_path, start=False
        )
        print("Reinvent Finished")

        # checks if the file "scaffold_memory.csv" exists in the specified directory
        # `slurm_path/results`. If the file exists, it uses the `scp` command to copy the file from the remote
        # server to the local `../data/` directory. It also copies the file "Agent.ckpt" from the remote
        # server to the local `reinvent_connect/resources/current_run/` directory.
        if check_if_results_exist(
            self.ssh_settings["ssh_login"], self.slurm_path, "scaffold_memory.csv"
        ):
            os.system(
                f"scp {self.ssh_settings['ssh_login']}:{self.slurm_path}/results/scaffold_memory.csv {output_path}"
            )
            os.system(
                f"scp {self.ssh_settings['ssh_login']}:{self.slurm_path}/results/Agent.ckpt {self.cwd}/input_files/current_run/"
            )


def generate_substruct_flag(substruct_dict):
    component_custom_alerts = {
        "component_type": "custom_alerts",
        "name": "Custom_alerts",
        "weight": 1,
        "specific_parameters": {"smiles": []},
    }

    component_list = []
    component_custom_alerts["name"] = "alters"
    component_custom_alerts["specific_parameters"]["smiles"] = sum(
        [list(substruct_dict[name]) for name in substruct_dict], []
    )
    component_list.append(copy.deepcopy(component_custom_alerts))
    #    for name in substruct_dict:
    #        component_custom_alerts["name"] = f"alerts"
    #        component_custom_alerts["specific_parameters"]["smiles"] = substruct_dict[name]
    #        component_list.append(copy.deepcopy(component_custom_alerts))

    return component_list


def generate_slurm_file(
    slurm_file_name: str, from_file: str, slurm_path: str, run_name: str
):
    """
    The function `generate_slurm_file` generates a SLURM file by reading from a template file and
    replacing placeholders with provided values.
    """
    with open(from_file) as f:
        text = f.read()
    with open(slurm_file_name, "w") as rsh:
        rsh.write(text.format(**locals()))


def start_remote_run(
    ssh_login: str,
    slurm_file: str,
    json_file: str,
    slurm_path: str,
    model: str = None,
    agent: str = None,
):
    """
    The function `start_remote_run` submits a job to a remote server using SLURM and transfers necessary
    files to the server.

    Args:
        ssh_login: (str):  ssh login you need to have a password free login
        slurm_file (str):  path to the SLURM script file that will be used to submit the job to the
                           remote cluster
        json_file (str):   path to a JSON file that contains Reinvent configuration settings
        slurm_path (str): path to the directory where the SLURM files and other necessary files will
                          be copied to on the remote machine
        model (str, optional): path to a file containing the model that will be used in the remote run.
                               Defaults to None.
        agent (str, optional): path to agent used in a remote run. Defaults to None.


    """
    print("Submitting Job")
    os.system(f"scp {slurm_file} {ssh_login}:{slurm_path}")
    os.system(f"scp {json_file} {ssh_login}:{slurm_path}")

    if model is not None:
        os.system(f"scp {model} {ssh_login}:{slurm_path}")
    if agent is not None:
        os.system(f"scp {agent} {ssh_login}:{slurm_path}")

    command2 = f"sbatch {slurm_path}/new_run.slurm"
    command = f"cd {slurm_path};{command2}"
    result = sb.run(
        ["ssh", ssh_login, command],
        shell=False,
        stdout=sb.PIPE,
        stderr=sb.PIPE,
        check=True,
    )
    return True


def check_reinvent_status(ssh_login, path, start=False, sleepTime=10):
    # use start to check whether the job has started
    # checks whether the scaffold_memory.csv has been removed
    finished = start
    while finished == start:
        finished = check_if_results_exist(ssh_login, path, file="memory.csv")
        print("Wait ...")
        time.sleep(sleepTime)
    return True


def check_if_results_exist(ssh_login, path, file):
    available = sb.call(["ssh", ssh_login, f"test -f {path}/results/{file}"])
    if available == 0:
        return True
    if available == 1:
        return False


def generate_inception(smiles):
    inception = {"memory_size": 20, "sample_size": 5, "smiles": smiles}
    return inception
