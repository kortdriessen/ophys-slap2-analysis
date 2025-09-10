
from scbc.slap2.experiment_summary import ExperimentSummary
import pathlib as pl
import argparse
import h5py
import numpy as np

from datetime import datetime, timezone
import subprocess


from aind_data_schema.components.identifiers import Code
from aind_data_schema.core.processing import (
    DataProcess,
    Processing,
    ProcessName,
    ProcessStage,
    ResourceTimestamped,
    ResourceUsage,
)
from aind_data_schema_models.units import MemoryUnit
from aind_data_schema_models.system_architecture import OperatingSystem, CPUArchitecture

def to_serializable(val):
    if isinstance(val, (np.integer, np.int32, np.int64, np.uint8)):
        return int(val)
    elif isinstance(val, (np.floating, np.float32, np.float64)):
        return float(val)
    elif isinstance(val, np.ndarray):
        return val.tolist()
    return val


parser = argparse.ArgumentParser(description="Convert experimentsummary.mat to experiment_summary.h5.")
parser.add_argument('--mat_path', type=str, help='Path to the input .mat file')
parser.add_argument('--output_dir', type=str, help='Directory to save the outputted processing.json')
args = parser.parse_args()

mat_path = pl.Path(args.mat_path)
if not mat_path.exists():
    raise FileNotFoundError(f"The specified .mat file does not exist: {mat_path}")
expsum = ExperimentSummary(mat_path)
expsum._get_metadata()
extraction_metadata = {key: to_serializable(val) for key, val in expsum.metadata.items()}
print('extraction metadata:', extraction_metadata)

with h5py.File(mat_path, 'r') as expsum_h5py:
    align_params = {key: to_serializable(val[0][0]) for key, val in expsum_h5py['exptSummary/trialTable/alignParams'].items()}
    print('align params:', align_params)

try:
    version = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], text=True
    ).strip()
except:
    version = "Unknown"

example_code = Code(
    url="https://github.com/abcd",
    version=version,
    parameters={"size": 7},
)

p = Processing.create_with_sequential_process_graph(
    pipelines=[
        Code(
            name="Matlab Slap2 processing pipeline",
            url="https://github.com/AllenNeuralDynamics/ophys-slap2-analysis/tree/main/matlab/preprocessing",
            version="0.1.0",
        ),
    ],
    data_processes=[
        DataProcess(
            process_type=ProcessName.VIDEO_ROI_TIMESERIES_EXTRACTION,
            experimenters=[extraction_metadata.get('experimenter_name', 'Unknown')],
            stage=ProcessStage.PROCESSING,
            start_date_time=extraction_metadata.get('start_time', datetime.now(timezone.utc)),
            end_date_time=extraction_metadata.get('end_time', datetime.now(timezone.utc)),
            output_path="",
            pipeline_name="Matlab Slap2 processing pipeline",
            code=example_code.model_copy(
                update=dict(
                    parameters=extraction_metadata,
                )
            ),
        ),
        DataProcess(
            process_type=ProcessName.VIDEO_MOTION_CORRECTION,
            pipeline_name="Matlab Slap2 processing pipeline",
            experimenters=[align_params.get('experimenter_name', 'Unknown')],
            stage=ProcessStage.PROCESSING,
            start_date_time=align_params.get('start_time', datetime.now(timezone.utc)),
            end_date_time=align_params.get('end_time', datetime.now(timezone.utc)),
            output_path="",
            code=example_code.model_copy(
                update=dict(
                    parameters=align_params
                )
            ),
        ),
    ],
)

serialized = p.model_dump_json()
deserialized = Processing.model_validate_json(serialized)
p.write_standard_file(pl.Path(args.output_dir))