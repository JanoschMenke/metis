from pydantic import BaseModel, Field, field_validator, validator
from typing import Dict, List, Optional
import re
import warnings
import yaml


class DataConfig(BaseModel):
    initial_path: str
    path: str
    selection_strategy: str
    num_molecules: int = Field(gt=0)
    run_name: str

    @field_validator("initial_path")
    def initial_path_validation(cls, v: str):
        if v.split("/")[-1] == "scaffold_memory.csv":
            raise ValueError('The initial data cannot be called "scaffold_memory.csv"')
        else:
            return v

    @field_validator("path")
    def path_validation(cls, v: str):
        if v.split("/")[-1] != "scaffold_memory.csv":
            raise ValueError(
                'The path specified in "path" must end in "scaffold_memory.csv"'
            )
        else:
            return v


class TabConfig(BaseModel):
    render: bool
    tab_names: List[str]


class NavButtonConfig(BaseModel):
    render: bool = Field(default=False)


class NavigationbarConfig(BaseModel):
    sendButton: NavButtonConfig
    editButton: NavButtonConfig


class GeneralConfig(BaseModel):
    render: bool = Field(default=True)
    slider: bool = Field(default=False)


class LiabilityConfig(BaseModel):
    name: str
    color: str

    @field_validator("color")
    def path_validation(cls, v: str):
        if re.search(r"^#(?:[0-9a-fA-F]{3}){1,2}$", v) is None:
            raise ValueError(
                f"The color {v} specified does not appear to be a valid hex color"
            )
        else:
            return v


class ECFPConfig(BaseModel):
    bitSize: int = Field(gt=0)
    radius: int = Field(gt=-1)
    useCounts: bool = False


class AdditionalWindowsConfig(BaseModel):
    render: bool = Field(default=False)
    path: Optional[str] = Field(default=None)
    ECFP: Optional[ECFPConfig] = Field(default=None)

    @validator("path", always=True)
    def path_validator(cls, v, values):
        if (v is not None) and (values["render"] == False):
            warnings.warn(
                f"Path {v} is specified for an additional window, but rendering for the Window is set to False"
            )
        elif (v is None) and (values["render"] == True):
            raise ValueError(
                f"A path must be specified if Window is set to be rendered."
            )
        return v


class SubstructureConfig(BaseModel):
    render: bool = Field(default=False)
    liabilities: Dict[str, LiabilityConfig]


class GlobalPropertiesConfig(BaseModel):
    render: bool = Field(default=False)
    liabilities: List[str]


class UIConfig(BaseModel):
    show_atom_contributions: AdditionalWindowsConfig = Field(
        default=AdditionalWindowsConfig()
    )
    show_reference_molecules: AdditionalWindowsConfig = Field(
        default=AdditionalWindowsConfig()
    )
    tab: TabConfig
    navigationbar: NavigationbarConfig
    general: GeneralConfig
    substructures: SubstructureConfig
    global_properties: GlobalPropertiesConfig


class RewardModelConfig(BaseModel):
    use_oracle_score: bool = True
    weight: Optional[str] = None
    oracle_path: Optional[str] = None
    qsar_model_path: str
    training_data_path: str
    ECFP: ECFPConfig


class DeNovoConfig(BaseModel):
    ssh_settings: str
    use_human_scoring_func: bool = False
    use_reward_model: bool = False


class BaseConfig(BaseModel):
    seed: Optional[int] = None
    tutorial: bool = False
    debug: bool = False
    max_iterations: int = Field(gt=0)
    innerloop_iterations: Optional[int] = Field(default=None, gt=0)
    activity_label: str
    introText: str
    propertyLabels: Dict[str, str]
    data: DataConfig
    ui: UIConfig
    de_novo_model: Optional[DeNovoConfig] = None
    reward_model: Optional[RewardModelConfig] = Field(default=None)

    @field_validator("innerloop_iterations")
    def innerloop_validator(cls, v=int):
        if (v == 0) or (v == 1):
            v = None
        return v

    @validator("de_novo_model")
    def de_novo_model_validator(cls, v):
        if (v.use_human_scoring_func == False) & (v.use_reward_model == False):
            warnings.warn(
                "No De Novo Model specified -> Feedback will not be included when Training Reinvent."
            )
        return v


def load_settings(settings_file: str) -> BaseConfig:
    """
    Load settings from a YAML file and validate them using the BaseConfig model.

    Args:
        settings_file (str): Path to the YAML settings file.

    Returns:
        BaseConfig: Validated settings object.

    Raises:
        ValueError: If the settings file is invalid or missing required fields.
    """
    try:
        with open(settings_file, "r") as f:
            settings_dict = yaml.safe_load(f)
        return BaseConfig(**settings_dict)
    except Exception as e:
        raise ValueError(f"Error loading settings: {str(e)}")


def save_settings(settings: BaseConfig, settings_file: str) -> None:
    """
    Save the current settings to a YAML file.

    Args:
        settings (BaseConfig): The current settings object.
        settings_file (str): Path to save the YAML settings file.
    """
    with open(settings_file, "w") as f:
        yaml.dump(settings.dict(), f)
