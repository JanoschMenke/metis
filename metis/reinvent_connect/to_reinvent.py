# %%
import subprocess as sb
import pandas as pd
from rdkit.Chem import AllChem as Chem
import json
import os
from utils.helper import fancySubstruct, get_random_string
from utils.data import extract_and_process_liabilities
import time
from PySide2.QtCore import QObject, Signal, Slot, QRunnable
import yaml
import copy


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
        self, activity_label, human_component=None
    ):  # , smiles, selectedAtoms):
        super(Worker, self).__init__()
        # self.smiles = smiles
        # self.selectedAtoms = selectedAtoms
        self.activity_label = activity_label
        self.dict = human_component
        self.signals = WorkerSignals()
        print("Worker Initialization")

    @Slot()
    def run(self):
        print("Run on Core")
        startNextRun(self.activity_label, self.dict)
        self.signals.finished.emit()


def startNextRun(activity_label, human_component=None):
    """
    The function `startNextRun()` updates certain file paths in a JSON file based on the existence of
    specific files in the current directory.
    """
    cwd = f"{os.getcwd()}/reinvent_connect"
    sshSettings = yaml.safe_load(open(f"{cwd}/input_files/ssh_settings.yml"))

    run_name = get_random_string(8)
    slurm_path = sshSettings["path_remote_folder"]
    agent_path = None
    model_path = None

    with open(f"{cwd}/input_files/test_qsar.json") as og_file:
        file_contents = json.loads(og_file.read())

    # check if new agent file exists, otherwise use inital agent on server
    if os.path.isfile(f"{cwd}/input_files/current_run/Agent.ckpt"):
        file_contents["parameters"]["reinforcement_learning"][
            "agent"
        ] = f"{slurm_path}/Agent.ckpt"
        agent_path = f"{cwd}/input_files/current_run/Agent.ckpt"
    else:
        file_contents["parameters"]["reinforcement_learning"][
            "agent"
        ] = f"{slurm_path}/Agent_Initial.ckpt"

    # check if model file exists else use intial model on server
    component_list = list()
    for component in file_contents["parameters"]["scoring_function"]["parameters"]:
        if component["name"] == activity_label:
            if os.path.isfile(f"{cwd}/input_files/current_run/Model.pkl"):
                component["specific_parameters"][
                    "model_path"
                ] = f"{slurm_path}/Model.pkl"
                model_path = f"{cwd}/input_files/current_run/Model.pkl"
            else:
                component["specific_parameters"][
                    "model_path"
                ] = f"{slurm_path}/Model_Initial.pkl"
        component_list.append(component)

    if human_component is not None:
        num_specified_components = len(component_list)
        num_human_components = len(human_component)
        scale_ratio = num_specified_components / num_human_components
        print(human_component)
        for comp in human_component:
            if comp is not None:
                comp["weight"] *= scale_ratio
                component_list.append(comp)

    file_contents["parameters"]["scoring_function"]["parameters"] = component_list
    with open(f"{cwd}/input_files/current_run/new_run.json", "w") as outfile:
        json.dump(file_contents, outfile, sort_keys=True, indent=4)

    generate_slurm_file(
        f"{cwd}/input_files/current_run/new_run.slurm",
        f"{cwd}/input_files/standard_slurm.slurm",
        slurm_path,
        run_name=run_name,
    )
    start_remote_run(
        sshSettings["ssh_login"],
        f"{cwd}/input_files/current_run/new_run.slurm",
        f"{cwd}/input_files/current_run/new_run.json",
        slurm_path,
        model_path,
        agent_path,
    )

    print("Slurmjob Submitted")
    check_reinvent_status(sshSettings["ssh_login"], slurm_path, start=True)
    print("Training has Started")
    check_reinvent_status(sshSettings["ssh_login"], slurm_path, start=False)
    print("Reinvent Finished")

    # checks if the file "scaffold_memory.csv" exists in the specified directory
    # `slurm_path/results`. If the file exists, it uses the `scp` command to copy the file from the remote
    # server to the local `../data/` directory. It also copies the file "Agent.ckpt" from the remote
    # server to the local `reinvent_connect/input_files/current_run/` directory.
    if check_if_results_exist(
        sshSettings["ssh_login"], slurm_path, "scaffold_memory.csv"
    ):
        os.system(
            f"scp {sshSettings['ssh_login']}:{slurm_path}/results/scaffold_memory.csv ../data/"
        )
        os.system(
            f"scp {sshSettings['ssh_login']}:{slurm_path}/results/Agent.ckpt reinvent_connect/input_files/current_run/"
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
