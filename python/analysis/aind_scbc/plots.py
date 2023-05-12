import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from matplotlib import pyplot as plt


def plt_matlab_traces(df: pd.DataFrame, file_num_ROIs: list[tuple[int,int]] = []):
    """
    Plots all or a selection of traces.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with columns 'time', 'trace1', 'n1', 'file number', 'ROI'

    file_num_ROIs : list[tuple[int,int]]
        List of tuples of (<file number>, <ROI>)
    """
    fig = go.Figure()

    # select only a portion of the data
    if len(file_num_ROIs):
        df = df.loc[df.apply(lambda x: (x["file number"], x["ROI"]) in file_num_ROIs, axis = 1)]

    for idx, row in df.iterrows():
        fig.add_trace(go.Scatter(x = row["time"], y = row["trace1"]/row["n1"]**0.5, mode = "lines",
                                 hovertemplate = "File: %s<br>"%row["file number"] +\
                                "ROI: %d"%row["ROI"], name = "%d/%d"%(row["file number"], row["ROI"])))
    fig.show()