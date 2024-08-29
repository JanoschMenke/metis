import copy
import json
import logging
import os
import time

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import paramiko
import yaml
from PySide6.QtCore import QObject, QRunnable, Signal, Slot

from metis import PKGDIR
from metis.utils.helper import get_random_string


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


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
        logger.debug("Worker Started De Novo Run")
        self.de_novo_runner.run(self.activity_label, self.output_path, self.dict)
        self.signals.finished.emit()
        logger.debug("Worker Finished De Novo Run")


class DeNovoRunner:
    def __init__(self, de_novo_config):
        self.cwd = f"{PKGDIR}/resources"
        self.settings = de_novo_config
        self.ssh_settings = yaml.safe_load(open(self.settings.ssh_settings))
        self.remote_executor = RemoteExecutor(
            hostname=self.ssh_settings["hostname"],
            username=self.ssh_settings["username"],
            password=self.ssh_settings.get("password"),
            key_filename=self.ssh_settings.get("key_filename"),
        )
        self.slurm_path = self.ssh_settings["path_remote_folder"]
        self.run_name = get_random_string(8)
        with open(self.ssh_settings["de_novo_json"]) as og_file:
            self.og_json_contents = json.loads(og_file.read())

        self.agent_path = None
        self.model_path = None
        logger.info("DeNovoRunner initialized")

    def gen_denovo_file(self, activity_label: str, human_component=None):
        logger.debug(f"Generating de novo file for activity label: {activity_label}")
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
        logger.debug("Reinvent .json file generated successfully")

    def run(self, activity_label: str, output_path: str, human_component=None):
        logger.info(f"Starting de novo run for activity label: {activity_label}")
        self.gen_denovo_file(activity_label, human_component=human_component)

        self.generate_slurm_file(
            f"{self.cwd}/input_files/current_run/new_run.slurm",
        )

        try:
            self.remote_executor.connect()
            local_files = [
                f"{self.cwd}/input_files/current_run/new_run.slurm",
                f"{self.cwd}/input_files/current_run/new_run.json",
                self.model_path,
                self.agent_path,
            ]
            local_files = [f for f in local_files if f is not None]
            self.remote_executor.transfer_files_to_remote(local_files, self.slurm_path)

            command = f"cd {self.slurm_path}; sbatch new_run.slurm"
            stdout, stderr = self.remote_executor.execute_remote_command(command)
            if stderr:
                logger.error(f"Error submitting job: {stderr}")
                return

            logger.info("Slurm job submitted")
            self.remote_executor.check_reinvent_status(start=True)
            logger.info("Training has started")
            self.remote_executor.check_reinvent_status(start=False)
            logger.info("Reinvent finished")

            if self.remote_executor.check_if_results_exist("scaffold_memory.csv"):
                self.remote_executor.transfer_files_from_remote(
                    [f"{self.slurm_path}/results/scaffold_memory.csv"], output_path
                )
                self.remote_executor.transfer_files_from_remote(
                    [f"{self.slurm_path}/results/Agent.ckpt"],
                    f"{self.cwd}/input_files/current_run/",
                )
                logger.info("Results transferred successfully")
            else:
                logger.warning("Results file not found on remote server.")

        finally:
            self.remote_executor.disconnect()
            logger.info("Remote connection closed")

    def generate_slurm_file(self, slurm_file_name):
        logger.debug(f"Generating SLURM file: {slurm_file_name}")
        with open(self.ssh_settings["default_slurm"]) as f:
            text = f.read()
        with open(slurm_file_name, "w") as rsh:
            rsh.write(text.format(slurm_path=self.slurm_path, run_name=self.run_name))
        logger.info("SLURM file generated successfully")


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
    return component_list


def generate_inception(smiles):
    inception = {"memory_size": 20, "sample_size": 5, "smiles": smiles}
    return inception


class RemoteExecutor:
    def __init__(
        self,
        hostname: str,
        username: str,
        password: Optional[str] = None,
        key_filename: Optional[str] = None,
        slurm_path: Optional[str] = None,
    ):
        self.hostname: str = hostname
        self.username: str = username
        self.password: Optional[str] = password
        self.key_filename: Optional[str] = key_filename
        self.ssh_client: Optional[paramiko.SSHClient] = None
        self.sftp_client: Optional[paramiko.SFTPClient] = None
        self.slurm_path: Optional[str] = slurm_path
        logger.info(f"RemoteExecutor initialized for {username}@{hostname}")

    def connect(self) -> None:
        logger.debug("Connecting to remote server")
        self.ssh_client = paramiko.SSHClient()
        self.ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())

        connect_kwargs: Dict[str, str] = {
            "hostname": self.hostname,
            "username": self.username,
        }
        if self.password:
            connect_kwargs["password"] = self.password
        elif self.key_filename:
            connect_kwargs["key_filename"] = self.key_filename

        self.ssh_client.connect(**connect_kwargs)
        self.sftp_client = self.ssh_client.open_sftp()
        logger.info("Connected to remote server")

    def disconnect(self) -> None:
        logger.debug("Disconnecting from remote server")
        if self.sftp_client:
            self.sftp_client.close()
        if self.ssh_client:
            self.ssh_client.close()
        logger.info("Disconnected from remote server")

    def transfer_files_to_remote(self, local_files: List[str], remote_dir: str) -> None:
        logger.debug(f"Transferring files to remote: {local_files}")
        for local_file in local_files:
            remote_path: str = f"{remote_dir}/{Path(local_file).name}"
            self.sftp_client.put(local_file, remote_path)
        logger.info(f"Transferred {len(local_files)} files to remote server")

    def transfer_files_from_remote(
        self, remote_files: List[str], local_dir: str
    ) -> None:
        logger.debug(f"Transferring files from remote: {remote_files}")
        for remote_file in remote_files:
            local_path: str = os.path.join(local_dir, Path(remote_file).name)
            self.sftp_client.get(remote_file, local_path)
        logger.info(f"Transferred {len(remote_files)} files from remote server")

    def execute_remote_command(self, command: str) -> Tuple[str, str]:
        logger.debug(f"Executing remote command: {command}")
        stdin, stdout, stderr = self.ssh_client.exec_command(command)
        return stdout.read().decode(), stderr.read().decode()

    def check_reinvent_status(self, start: bool = False, sleep_time: int = 10) -> bool:
        logger.debug(f"Checking Reinvent status, start={start}")
        finished: bool = start
        while finished == start:
            finished = self.check_if_results_exist("memory.csv")
            logger.info("Waiting for Reinvent to finish...")
            time.sleep(sleep_time)
        return True

    def check_if_results_exist(self, file: str) -> bool:
        logger.debug(f"Checking if file exists on remote: {file}")
        stdout, stderr = self.execute_remote_command(
            f"test -f {self.slurm_path}/results/{file}"
        )
        return stderr == ""
