import pandas as pd, os, numpy as np
from scipy.io import loadmat
import analysis_an.process as an_proc  # fluorescence debleaching
from oasis.functions import deconvolve as oasis_deconvolve # oasis event detection
from collections import Iterable
import bottleneck as bn # speeds up certain operations on numpy arrays with nan


def read_data(dsets: str, dstore: str):
    """
    Loads extracted fluorescence traces from .mat files given a .xls spreadsheet containing a selection of datasets.

    The .mat files should contain these fields:

    sData.trace0 - raw ROI without decorrelating motion variables
    sData.n0 - noise for trace0
    sData.trace1 - raw ROI with motion variables decorrelated
    sData.n1 - noise for trace1
    sData.trace2 - NMF-based traces
    sData.n2 - noise for trace2
    sData.bleach - bleaching curve for this FOV downsampled time
    sData.frametime - sampling time in [s]
    sData.dsFac - downsampling factor
    sData.names - ROI names

    The .xls file should have at least these these columns:
    1. file - path to .mat file with extracted traces relative to dstore
    2. ROI field - number of ROI field (int)
    3. include - whether to include data from this file
    4. exclude ROI - exclude a single patch ROI from the scanfield (int)

    Parameters
    ----------
    dsets: str
        Path to .xls sheet
    dstore: str
        Folder containing recordings named as slap2_<mouse num>_<rec date>_<rec time>

    Returns
    =======
    pd.DataFrame
    """
    dsets_df = pd.read_excel(dsets)
    
    data = pd.concat([pd.DataFrame([
        {
            "file number": ds_row["file number"], 
            "file name": ds_row["file name"],
            "ROI field": ds_row["ROI field"],
            "ROI": roi,
            "frametime": m["sData"]["frametime"][0,0][0,0],
            "dsFac": m["sData"]["dsFac"][0,0][0,0],
            "trace0": m["sData"]["trace0"][0,0][:,roi_idx],
            "n0": m["sData"]["n0"][0,0][0,roi_idx],
            "trace1": m["sData"]["trace1"][0,0][:,roi_idx],
            "n1": m["sData"]["n1"][0,0][0,roi_idx],
            "trace2": m["sData"]["trace2"][0,0][:,roi_idx],
            "n2": m["sData"]["n2"][0,0][0,roi_idx],
            "bleach": m["sData"]["bleach"][0,0][0,:],
            "time": np.arange(len(m["sData"]["trace0"][0,0][:,roi_idx]))*m["sData"]["frametime"][0,0][0,0],
            "dstime": np.arange(len(m["sData"]["bleach"][0,0][0,:]))*m["sData"]["frametime"][0,0][0,0]*m["sData"]["dsFac"][0,0][0,0]
        } for roi_idx,roi in enumerate(m["sData"]["names"][0,0].astype(int).flatten())]) for m, ds_row in 
                [(loadmat(os.path.join(dstore, ds_row["file name"])), ds_row) for _, ds_row in dsets_df.iterrows() if ds_row["include"]]]
    )
    # exclude ROIs from recordings
    drop_rois = [(row["file name"], int(row["exclude ROI"])) for row_idx, row in dsets_df[~dsets_df["exclude ROI"].isnull()].iterrows()]
    for x in drop_rois:
        drop_idxs = data[(data["file name"] == x[0]) & (data["ROI"] == x[1])].index
        data.drop(drop_idxs, inplace = True)
    
    return data

def df_to_feather(df: pd.DataFrame, fpath: str):
    """
    Save a DataFrame to feather format.

    Parameters
    ----------
    df : pd.DataFrame
        pandas DataFrame.
    fpath : str
        File path and name where to store the data.
    """
    df.reset_index(drop = True).to_feather(fpath)
