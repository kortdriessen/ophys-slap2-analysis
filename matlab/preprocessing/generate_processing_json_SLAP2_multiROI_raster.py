
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
extraction_metadata = {key: to_serializable(val) for key, val in expsum.metadata.items()}
print('extraction metadata:', extraction_metadata)

align_params = {key: to_serializable(val) for key, val in expsum.align_params.items()}
print('align params:', align_params)

try:
    version = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], text=True
    ).strip()
except:
    version = "Unknown"

example_code = Code(
    url="https://github.com/AllenNeuralDynamics/ophys-slap2-analysis",
    version=version,
    parameters={"size": 7},
)

p = Processing.create_with_sequential_process_graph(
    pipelines=[
        Code(
            name="SLAP2 multi-ROI raster processing pipeline (Matlab)",
            url="https://github.com/AllenNeuralDynamics/ophys-slap2-analysis/blob/main/matlab/preprocessing/processSLAP2.m",
            version=version,
        ),
    ],
    data_processes=[
        DataProcess(
            process_type=ProcessName.VIDEO_MOTION_CORRECTION,
            pipeline_name="SLAP2 multi-ROI raster processing pipeline (Matlab)",
            experimenters=[align_params.get('operator', 'Unknown')],
            stage=ProcessStage.PROCESSING,
            start_date_time=datetime.fromisoformat(align_params.get('startTime', datetime.now(timezone.utc).isoformat())),
            end_date_time=datetime.fromisoformat(align_params.get('endTime', datetime.now(timezone.utc).isoformat())),
            output_path="",
            code=example_code.model_copy(
                update=dict(
                    url="https://github.com/AllenNeuralDynamics/ophys-slap2-analysis/blob/main/matlab/preprocessing/multiRoiRegSLAP2.m",
                    parameters=align_params
                )
            ),
        ),
        DataProcess(
            process_type=ProcessName.VIDEO_ROI_TIMESERIES_EXTRACTION,
            experimenters=[extraction_metadata.get('operator', 'Unknown')],
            stage=ProcessStage.PROCESSING,
            start_date_time=datetime.fromisoformat(extraction_metadata.get('startTime', datetime.now(timezone.utc).isoformat())),
            end_date_time=datetime.fromisoformat(extraction_metadata.get('endTime', datetime.now(timezone.utc).isoformat())),
            output_path="",
            pipeline_name="SLAP2 multi-ROI raster processing pipeline (Matlab)",
            code=example_code.model_copy(
                update=dict(
                    url="https://github.com/AllenNeuralDynamics/ophys-slap2-analysis/blob/main/matlab/preprocessing/summarize_NoLoCo.m",
                    parameters=extraction_metadata,
                )
            ),
        ),
    ],
)

serialized = p.model_dump_json()
deserialized = Processing.model_validate_json(serialized)
p.write_standard_file(pl.Path(args.output_dir))