import os
import pandas as pd
from analysis_utils.io import load_tree_data
import constants as C

FILTERED_CACHE_DIR = "df_cache_filtered"


def filter_df(df):
    mask = df["zum"] >= C.FILTER_DEPTH_UM
    if "nInteractions" in df.columns:
        mask &= df["nInteractions"] == 1
    x = df["xum"]
    y = df["yum"]
    for xmin, xmax, ymin, ymax in C.FILTER_REGIONS_EXCLUDE_XY_UM:
        mask &= ~((x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax))
    return df[mask].reset_index(drop=True)


def load_filtered(filename):
    filepath = f"{C.ROOT_FILES_DIR}/{filename}.root"
    cache_path = os.path.join(FILTERED_CACHE_DIR, f"{filename}.pkl")
    os.makedirs(FILTERED_CACHE_DIR, exist_ok=True)
    if os.path.exists(cache_path):
        src_mtime = os.path.getmtime(filepath)
        if src_mtime <= os.path.getmtime(cache_path):
            print(f"Loading cached filtered DataFrame: {cache_path}")
            return pd.read_pickle(cache_path)
    df = load_tree_data(filepath, tree_name="bef_tree", cache_dir=None)
    df = filter_df(df)
    df.to_pickle(cache_path)
    print(f"Cached filtered DataFrame to {cache_path}")
    return df
