from pydantic import BaseModel, Field, field_validator
from typing import Dict, List, Optional
import re


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
    compareButton: NavButtonConfig


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


class SubstructureConfig(BaseModel):
    render: bool = Field(default=False)
    liabilities: Dict[str, LiabilityConfig]


class GlobalPropertiesConfig(BaseModel):
    render: bool = Field(default=False)
    liabilities: List[str]


class UIConfig(BaseModel):
    show_atom_contributions: bool
    show_reference_molecules: bool
    tab: TabConfig
    navigationbar: NavigationbarConfig
    general: GeneralConfig
    substructures: SubstructureConfig
    global_properties: GlobalPropertiesConfig


class ECFPConfig(BaseModel):
    bitSize: int = Field(gt=0)
    radius: int = Field(gt=-1)
    useCounts: bool


class InteractiveModelConfig(BaseModel):
    oracle_score: bool = True
    weight: str = "pseudo_confidence"
    oracle_path: str
    qsar_model_path: str
    training_data_path: str
    ECFP: ECFPConfig


class BaseConfig(BaseModel):
    tutorial: bool
    debug: bool
    max_iterations: int = Field(gt=0)
    innerloop_iterations: Optional[int] = Field(default=None, gt=0)
    use_human_component: bool
    activity_label: str
    introText: str
    propertyLabels: Dict[str, str]
    data: DataConfig
    ui: UIConfig
    interactive_model: InteractiveModelConfig

    @field_validator("innerloop_iterations")
    def innerloop_validator(cls, v=int):
        if (v == 0) or (v == 1):
            v = None
        return v
