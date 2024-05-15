# Metis - A Python-Based User Interface to Collect Expert Feedback for Generative Chemistry Models
---

Metis is GUI to enable the collection of accurate and detailed feedback on small molecules.
At its core, it is built around [Esben Bjerrums](https://github.com/EBjerrum) [rdEditor](https://github.com/EBjerrum/rdeditor) using PySide2.

You can find the preprint at [ChemRxiv](https://chemrxiv.org/engage/chemrxiv/article-details/66421031418a5379b0255d8a)

---

**Table of Contents**

- [Metis](#metis)
  * [Set Up](#set-up)
    + [Installation](#installation)
      - [Dependencies](#dependencies)
    + [SSH](#ssh)
  * [Usage](#usage)
    + [Examples](#examples)
      - [UI Only](#ui-only)
      - [Reward Model](#reward-model)
      - [De Novo Design](#de-novo-design)
  * [Settings](#settings)

---

## Set up

### Installation

Download the repository and navigate to the download location. You can install `metis` with `pip install .`. Make sure the environment you want to install into is activated.

#### Dependencies
Some notes on the dependencies. 

**scikit-learn**

The version `scikit-learn` constraints are only set to make sure that the examples given here work.
In theory, you could use any `scikit-learn` version. If you want to use Reinvent in the backend, you need to make sure that the version of `scikit-learn` Reinvent is using on the remote machine
should be updated to the version that matches your local installation used by metis.

**cairosvg**

Depending on the OS you are running installing `cairosvg` through pip can cause issues, as cairo is not found. On MacOS you can solve this by installing cairo using homebrew, or you can install `cairosvg` using conda-


### SSH
It is assumed you have a working version of Reinvent on a Server instance that is running Slurm and ssh.

1. Change the ssh settings in the `example_project/de_novo_files/ssh_settings.yml` file.
      - `ssh_login`: your login to SSH e.g. `username@remote_server` *you should be able to access your remote server without a password, for example, using an RSA Key*
      - `path_remote_folder`: path on the remote machine, from where Reinvent files will be loaded and stored.
      - `de_novo_json`: specify which default `reinvent.json` file to use
      - `default_slurm`: specify which default Slurm job to use

2. Copy and unzip the `metis_reinvent.zip` to the remote machine. Make sure that the `path_remote_folder` in the `ssh_settings.yml` file matches with the folder location and also in the `initial_reinvent.json`.


## Usage

After installation simply run:
```
metis -f path/to/settings.yml --output /path/where/to/save/
```
This will start the GUI. Examples can be found below.

### Examples

#### UI Only

In the most simple example, only the GUI will be started to collect feedback. No models are trained and no de novo run started.

```diff
- If you want to show the atom contributions to the predictions/model explanation
- (show_atom_contributions: render: true)
- you will experience heavy slowdowns when switching to a new molecule.
- The only solution at the moment is not to show them.
- You can set show_atom_contributions: render: False.
- This will yield a much smoother experience.    
```

```
cd example_project
metis -f settings_ui.yml --output results/
```

#### Reward Model

Here, next to collecting feedback, a reward model is also trained on the feedback. For this, we provided a QSAR model and Oracle model for JNK3 activity.
The setting `use_oracle_score: False`, will use the feedback of humans as the target variable that is to be predicted. If the setting is set to `True`, the molecules liked by the chemist will be scored by the oracle, and these scores will then be used as the target varible for the reward model. This can be thought of as an active learning setting, where the chemists decides which molecules are being "biologically validated".


```
cd example_project
metis -f settings_reward_model.yml --output results/
```

#### De Novo Design

With these settings, a REINVENT de novo run can be started directly using `Metis` on a remote machine.
The remote machine needs:
- a working installation of REINVENT 3. 
- update the REINVENTS scikit-learn to >1.0.0
- Wlurm
- access through SSH wih a key
- the unzipped `metis_reinvent.zip` folder

Once copied and unzipped, the paths and settings in the `de_novo_files` folder need to be adapted to fit to your paths on the remote machine.
```
cd example_project
metis -f settings_denovo.yml --output results/
```




## Settings

Here is a brief overview of all settings


| Name                     | Type                                                 | Required   | Default   |
|--------------------------|------------------------------------------------------|------------|-----------|
| <strong>seed</strong>                 | Union[int, None]                                     | False      |           |
| <strong>tutorial</strong>             | bool                                                 | False      | False     |
| <strong>debug</strong>                | bool                                                 | False      | False     |
| <strong>max_iterations</strong>       | int                                                  | True       | ...       |
| <strong>innerloop_iterations</strong> | Union[int, None]                                     | False      | None      |
| <strong>activity_label</strong>       | str                                                  | True       | ...       |
| <strong>introText</strong>            | str                                                  | True       | ...       |
| <strong>propertyLabels</strong>       | Dict                                                 | True       | ...       |
| <strong>data</strong>                 | <a href="#dataconfig">DataConfig</a>                            | True       | ...       |
| <strong>ui</strong>                   | <a href="#uiconfig">UIConfig</a>                                | True       | ...       |
| <strong>de_novo_model</strong>        | Union[<a href="#denovoconfig">DeNovoConfig</a>, None]           | False      | None      |
| <strong>reward_model</strong>         | Union[<a href="#rewardmodelconfig">RewardModelConfig</a>, None] | False      | None      |

- `debug`: if `True` will overwrite existing results folders
- `max_iterations` defines how often molecules are sampled, feedback collected and the model updated
- `innerloop_iteration` how often molecules are resampled from the same scaffold memory before the model is sent to the remote machine


<h3>DataConfig</h3>

| Name                   | Type   | Required   | Default   |
|------------------------|--------|------------|-----------|
| <strong>initial_path</strong>       | str    | True       | ...       |
| <strong>path</strong>               | str    | True       | ...       |
| <strong>selection_strategy</strong> | str    | True       | ...       |
| <strong>num_molecules</strong>      | int    | True       | ...       |
| <strong>run_name</strong>           | str    | True       | ...       |

- `initial_path`: path to inital dataset, the molecules that shall be evaluated first
- `path`:   path to subsequent datasets, these come from the server, and are generated by Reinvent, should end in `scaffold_memory.csv`
- `selection_strategy`: how to pick which molecules to show
- `num_molecules`: how many molecules to show
- `run_name`: what is the name of the run, under this name the results will be stored 


<h3>UIConfig</h3>

| Name                         | Type                                                | Required   | Default                                       |
|------------------------------|-----------------------------------------------------|------------|-----------------------------------------------|
| <strong>show_atom_contributions</strong>  | <a href="#additionalwindowsconfig">AdditionalWindowsConfig</a> | False      | {'render': False, 'path': None, 'ECFP': None} |
| <strong>show_reference_molecules</strong> | <a href="#additionalwindowsconfig">AdditionalWindowsConfig</a> | False      | {'render': False, 'path': None, 'ECFP': None} |
| <strong>tab</strong>                      | <a href="#tabconfig">TabConfig</a>                             | True       | ...                                           |
| <strong>navigationbar</strong>            | <a href="#navigationbarconfig">NavigationbarConfig</a>         | True       | ...                                           |
| <strong>general</strong>                  | <a href="#generalconfig">GeneralConfig</a>                     | True       | ...                                           |
| <strong>substructures</strong>            | <a href="#substructureconfig">SubstructureConfig</a>           | True       | ...                                           |
| <strong>global_properties</strong>        | <a href="#globalpropertiesconfig">GlobalPropertiesConfig</a>   | True       | ...                                           |

<h3>AdditionalWindowsConfig</h3>

| Name       | Type                                   | Required   | Default   |
|------------|----------------------------------------|------------|-----------|
| <strong>render</strong> | bool                                   | False      | False     |
| <strong>path</strong>   | Union[str, None]                       | False      |           |
| <strong>ECFP</strong>   | Union[<a href="#ecfpconfig">ECFPConfig</a>, None] | False      |           |</p>

<h3>ECFPConfig

| Name          | Type   | Required   | Default   |
|---------------|--------|------------|-----------|
| <strong>bitSize</strong>   | int    | True       | ...       |
| <strong>radius</strong>    | int    | True       | ...       |
| <strong>useCounts</strong> | bool   | False      | False     |

<h3>TabConfig</h3>

| Name          | Type   | Required   | Default   |
|---------------|--------|------------|-----------|
| <strong>render</strong>    | bool   | True       | ...       |
| <strong>tab_names</strong> | List   | True       | ...       |

- `render` if `False` it will not render the additional tabs

<h3>NavigationbarConfig</h3>

| Name           | Type                                | Required   | Default   |
|----------------|-------------------------------------|------------|-----------|
| <strong>sendButton</strong> | <a href="#navbuttonconfig">NavButtonConfig</a> | True       | ...       |
| <strong>editButton</strong> | <a href="#navbuttonconfig">NavButtonConfig</a> | True       | ...       |

<h3>NavButtonConfig</h3>

| Name       | Type   | Required   | Default   |
|------------|--------|------------|-----------|
| <strong>render</strong> | bool   | False      | False     |


<h3>GeneralConfig</h3>

| Name       | Type   | Required   | Default   |
|------------|--------|------------|-----------|
| <strong>render</strong> | bool   | False      | True      |
| <strong>slider</strong> | bool   | False      | False     |

<h3>SubstructureConfig</h3>

| Name            | Type   | Required   | Default   |
|-----------------|--------|------------|-----------|
| <strong>render</strong>      | bool   | False      | False     |
| <strong>liabilities</strong> | Dict   | True       | ...       |

Liablities control which properties you can select substructures for:
Keys such as `ugly` or `tox` are simply used within the script.
`name` will define how the button is called 
`color` will define the color of the button as well as the color of the atom highlight

```
liabilities:
      ugly:
        name: "Mutagenicity"
        color: "#ff7f7f"
      tox:
        name: "Toxicity" 
        color: "#51d67e"
      stability:
        name: "Stability"
        color: "#eed358"
      like:
        name: "Good"
        color: "#9542f5"
```


<h3>GlobalPropertiesConfig</h3>

| Name            | Type   | Required   | Default   |
|-----------------|--------|------------|-----------|
| <strong>render</strong>      | bool   | False      | False     |
| <strong>liabilities</strong> | List   | True       | ...       |

<h3>DeNovoConfig</h3>

| Name                       | Type   | Required   | Default   |
|----------------------------|--------|------------|-----------|
| <strong>ssh_settings</strong>           | str    | True       | ...       |
| <strong>use_human_scoring_func</strong> | bool   | False      | False     |
| <strong>use_reward_model</strong>       | bool   | False      | False     |


<h3>RewardModelConfig</h3>

| Name                   | Type                      | Required   | Default   |
|------------------------|---------------------------|------------|-----------|
| <strong>use_oracle_score</strong>   | bool                      | False      | True      |
| <strong>weight</strong>             | Union[str, None]          | False      | None      |
| <strong>oracle_path</strong>        | Union[str, None]          | False      | None      |
| <strong>qsar_model_path</strong>    | str                       | True       | ...       |
| <strong>training_data_path</strong> | str                       | True       | ...       |
| <strong>ECFP</strong>               | <a href="#ecfpconfig">ECFPConfig</a> | True       | ...       |

- `use_oracle_score` instead of using the feedback directly to train the reward model, one can use the oracle model to score molecules liked by the user. The reward model is then trained on the predictions of the oracle rather than on the direct feedback. This mimics an active learning scenario where the chemist can choose which molecules he wants to biologically validate



