import os
from pathlib import Path

import pandas as pd
import seqdata as sd

from ._utils import get_download_path, try_download_urls

HERE = Path(__file__).parent
pkg_resources = None

def get_dataset_info():
    """
    Return a pandas DataFrame with info about EUGENe's built in datasets.

    Returns
    -------
    df : pd.DataFrame
        Info about builtin datasets indexed by dataset name.
    """
    global pkg_resources
    if pkg_resources is None:
        import pkg_resources
    stream = pkg_resources.resource_stream(__name__, "datasets.csv")
    return pd.read_csv(stream, index_col=0)

def random1000():
    """
    Reads in the built-in random1000 dataset into a SeqData object.

    This dataset is designed for testing and development purposes.
    It contains 1000 random sequences of length 100.

    Returns
    -------
    sdata : SeqData
        SeqData object for the random1000 dataset.
    """
    infile = f"{HERE}/random1000/random1000_seqs.tsv"
    outzarr = f"{HERE}/random1000/random1000_seqs.zarr"
    sdata = sd.read_table(
        name="seq",
        tables=infile,
        out=outzarr,
        seq_col="seq",
        batch_size=1000,
        overwrite=True
    )
    return sdata

def jores21(
    dataset="leaf", 
    download_dir: str = None,
    batch_size: int = 1000,
    fixed_length=False,
    return_sdata: bool = True
):
    """
    Reads in the jores21 dataset into a SeqData object.

    Parameters
    ----------
    dataset : str, optional
        Dataset to read. Either "leaf" or "proto". The default is "leaf".
        The "leaf" dataset is from promoters tested in tobacco leaves
        and the "proto" dataset is from promoters tested in maize protoplasts.
    add_metadata : bool, optional
        If True, add metadata to the SeqData object. The default is False.
    return_sdata : bool, optional
        If True, return SeqData object for the jores21 dataset. The default is True.
        If False, return a list of paths to the downloaded files.
    **kwargs : kwargs, dict
        Keyword arguments to pass to read_csv.

    Returns
    -------
    sdata : SeqData
        SeqData object for the Jores21 dataset.
    """
    urls_list = [
        "https://raw.githubusercontent.com/tobjores/Synthetic-Promoter-Designs-Enabled-by-a-Comprehensive-Analysis-of-Plant-Core-Promoters/main/CNN/CNN_test_leaf.tsv",
        "https://raw.githubusercontent.com/tobjores/Synthetic-Promoter-Designs-Enabled-by-a-Comprehensive-Analysis-of-Plant-Core-Promoters/main/CNN/CNN_train_leaf.tsv",
        "https://raw.githubusercontent.com/tobjores/Synthetic-Promoter-Designs-Enabled-by-a-Comprehensive-Analysis-of-Plant-Core-Promoters/main/CNN/CNN_train_proto.tsv",
        "https://raw.githubusercontent.com/tobjores/Synthetic-Promoter-Designs-Enabled-by-a-Comprehensive-Analysis-of-Plant-Core-Promoters/main/CNN/CNN_test_proto.tsv",
        "https://static-content.springer.com/esm/art%3A10.1038%2Fs41477-021-00932-y/MediaObjects/41477_2021_932_MOESM3_ESM.xlsx",
    ]
    if dataset == "leaf":
        data_idxs = [0, 1]
    elif dataset == "proto":
        data_idxs = [2, 3]
    else:
        raise ValueError("dataset must be either 'leaf' or 'proto'.")
    download_dir = get_download_path() if download_dir is None else download_dir
    paths = try_download_urls(data_idxs, urls_list, download_dir, "jores21")
    if return_sdata:
        outzarr = os.path.join(download_dir, "jores21", f"jores21_{dataset}.zarr")
        sdata = sd.read_table(
            name="seq",
            tables=paths,
            out=outzarr,
            seq_col="sequence",
            fixed_length=fixed_length,
            batch_size=batch_size,
            overwrite=True
        )
        return sdata
    else:
        return paths
    