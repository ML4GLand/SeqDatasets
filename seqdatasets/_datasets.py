import os
from pathlib import Path

import xarray as xr
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
    It contains 1000 random sequences of length 100 as well as 10
    random continuouts value "activities" between 0 and 1 and 10
    binary "labels" of 0 or 1.

    Returns
    -------
    sdata : SeqData
        SeqData object for the random1000 dataset.
    """
    infile = f"{HERE}/random1000/random1000_seqs.tsv"
    outzarr = f"{HERE}/random1000/random1000_seqs.zarr"
    if not os.path.exists(outzarr):
        sdata = sd.read_table(
            name="seq",
            tables=infile,
            out=outzarr,
            seq_col="seq",
            fixed_length=True,
            batch_size=1000,
            overwrite=True
        )
    else:
        sdata = sd.open_zarr(outzarr)
    return sdata

def jores21(
    dataset="leaf", 
    download_dir: str = None,
    batch_size: int = 1000,
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
    download_dir : str, optional
        Directory to download the dataset to. The default is None.
        If EUGENe is installed, the default is the EUGENe download directory.
        If EUGENe is not installed, the default is the current working directory.
    batch_size : int, optional
        Batch size to use when reading in the dataset. The default is 1000.
        This controls the number of sequences read in at a time and should be 
        adjusted based on the amount of memory available.
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
        if not os.path.exists(outzarr):
            print("Zarr file not found. Creating new zarr file.")
            sdata = sd.read_table(
                name="seq",
                tables=paths,
                out=outzarr,
                seq_col="sequence",
                fixed_length=True,
                batch_size=batch_size,
                overwrite=True
            )
        else:
            print("Zarr file found. Opening zarr file.")
            sdata = sd.open_zarr(outzarr)
        return sdata
    else:
        return paths
    
def ray13(
    dataset="norm", 
    download_dir: str = None,
    batch_size: int = 1000,
    return_sdata=True
):
    """
    Reads in the ray13 dataset into a SeqData object.

    Parameters
    ----------
    dataset : str
        Dataset to read, can either be "norm" or "raw". The default is "norm".
        The "norm" dataset is the normalized probe intensities from Ray et al 2013
        and the "raw" dataset is the raw probe intensities.
    download_dir : str, optional
        Directory to download the dataset to. The default is None.
        If EUGENe is installed, the default is the EUGENe download directory.
        If EUGENe is not installed, the default is the current working directory.
    return_sdata : bool, optional
        If True, return SeqData object for the ray13 dataset. 
        If False, return a list of paths to the downloaded files. The default is True.
    **kwargs : kwargs, dict
        Keyword arguments to pass to read_csv.

    Returns
    -------
    sdata : SeqData
        SeqData object for the ray13 dataset.
    """
    urls_list = [
        "http://hugheslab.ccbr.utoronto.ca/supplementary-data/RNAcompete_eukarya/norm_data.txt.gz",
        "http://hugheslab.ccbr.utoronto.ca/supplementary-data/RNAcompete_eukarya/raw_data.txt.gz",
    ]

    if dataset == "norm":
        data_idxs = [0]
    elif dataset == "raw":
        data_idxs = [1]
    else:
        raise ValueError("dataset must be either 'norm' or 'raw'.")
    
    download_dir = get_download_path() if download_dir is None else download_dir
    paths = try_download_urls(data_idxs, urls_list, download_dir, "ray13")

    # Gunzip the downloaded file if has .gz
    if return_sdata:
        outzarr = os.path.join(download_dir, "ray13", f"ray13_{dataset}.zarr")
        if not os.path.exists(outzarr):
            print("Zarr file not found. Creating new zarr file.")
            sdata = sd.read_table(
                name="seq",
                tables=paths,
                out=outzarr,
                seq_col="RNA_Seq",
                fixed_length=False,
                batch_size=batch_size,
                overwrite=True,
                delim_whitespace=True
            )
        else:
            print("Zarr file found. Opening existing zarr file.")
            sdata = sd.open_zarr(outzarr)
        return sdata
    else:
        return paths

def deBoer20(
    datasets: list, 
    download_dir: str = None,
    batch_size: int = 1000,
    return_sdata=True
):
    """
    Reads in the deBoer20 dataset.

    Parameters
    ----------
    datasets : list of ints
        List of datasets indices to read. There are 8 in total (0-7), each from:
        https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104878 
    return_sdata : bool, optional
        If True, return SeqData object for the deBoer20 dataset. The default is True.
        If False, return a list of paths to the downloaded files.
    **kwargs : kwargs, dict
        Keyword arguments to pass to read_csv.

    Returns
    -------
    sdata : SeqData
        SeqData object for the deBoer20 dataset.
    """
    urls_list = [
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104878/suppl/GSE104878_20160503_average_promoter_ELs_per_seq_atLeast100Counts.txt.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104878/suppl/GSE104878_20160609_average_promoter_ELs_per_seq_Abf1TATA_ALL.shuffled.txt.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104878/suppl/GSE104878_20160609_average_promoter_ELs_per_seq_pTpA_ALL.shuffled.txt.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104878/suppl/GSE104878_20161024_average_promoter_ELs_per_seq_3p1E7_Gal_ALL.shuffled.txt.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104878/suppl/GSE104878_20161024_average_promoter_ELs_per_seq_3p1E7_Gly_ALL.shuffled.txt.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104878/suppl/GSE104878_20170811_average_promoter_ELs_per_seq_OLS_Glu_goodCores_ALL.txt.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104878/suppl/GSE104878_20180808_processed_Native80_and_N80_spikein.txt.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104878/suppl/GSE104878_pTpA_random_design_tiling_etc_YPD_expression.txt.gz",
    ]

    if type(datasets) is int:
        data_idxs = [datasets]

    download_dir = get_download_path() if download_dir is None else download_dir
    paths = try_download_urls(data_idxs, urls_list, download_dir, "deBoer20")

    if return_sdata:
        dataset_str = "-".join([str(i) for i in data_idxs])
        outzarr = os.path.join(download_dir, "deBoer20", f"deBoer20_datasets{dataset_str}.zarr")
        if not os.path.exists(outzarr):
            print("Zarr file not found. Creating new zarr file.")
            sdata = sd.read_table(
                name="seq",
                tables=paths,
                out=outzarr,
                seq_col="seq",
                fixed_length=False,
                batch_size=batch_size,
                overwrite=True,
                header=None,
                names=["seq", "target"],
            )
        else:
            print("Zarr file found. Opening existing zarr file.")
            sdata = sd.open_zarr(outzarr)
        return sdata
    else:
        return paths
    
def deAlmeida22(
    dataset="train", 
    download_dir: str = None,
    batch_size: int = 1000,
    return_sdata=True, 
):
    """
    Reads the deAlmeida22 dataset into a SeqData object.

    Parameters
    ----------
    dataset : str, optional
        Dataset to read. Either "train", "val", "test". The default is "train".
        These are the different sets used in the deAlmeida et al 2022 paper.
    return_sdata : bool, optional
        If True, return SeqData object for the deAlmeida22 dataset. The default is True.
        If False, return a list of paths to the downloaded files.
    **kwargs : kwargs, dict
        Keyword arguments to pass to read_csv.

    Returns
    -------
    sdata : SeqData
        SeqData object for the deAlmeida22 dataset.
    """
    urls_list = [
        "https://zenodo.org/record/5502060/files/Sequences_Train.fa?download=1",
        "https://zenodo.org/record/5502060/files/Sequences_Val.fa?download=1",
        "https://zenodo.org/record/5502060/files/Sequences_Test.fa?download=1",
        "https://zenodo.org/record/5502060/files/Sequences_activity_Train.txt?download=1",
        "https://zenodo.org/record/5502060/files/Sequences_activity_Val.txt?download=1",
        "https://zenodo.org/record/5502060/files/Sequences_activity_Test.txt?download=1",
    ]
    if dataset == "train":
        data_idxs = [0, 3]
    elif dataset == "val":
        data_idxs = [1, 4]
    elif dataset == "test":
        data_idxs = [2, 5]

    download_dir = get_download_path() if download_dir is None else download_dir
    paths = try_download_urls(data_idxs, urls_list, download_dir, "deAlmeida22")

    if return_sdata:
        outzarr = os.path.join(download_dir, "deAlmeida22", f"deAlmeida22_{dataset}.zarr")
        if not os.path.exists(outzarr):
            print("Zarr file not found. Creating new zarr file.")
            sdata = sd.read_flat_fasta(
                name="seq",
                out=outzarr,
                fasta=paths[0],
                batch_size=batch_size,
                fixed_length=True,
                overwrite=True
            )
            metadata = pd.read_csv(paths[1], sep="\t")
            meta_xr = xr.Dataset.from_dataframe(metadata)
            meta_xr = meta_xr.drop_vars("index").rename_dims({"index": "_sequence"}).assign(_sequence=sdata["_sequence"])
            sdata = sdata.merge(meta_xr)
            sd.to_zarr(sdata, outzarr, mode="w")
        else:
            print("Zarr file found. Opening existing zarr file.")
            sdata = sd.open_zarr(outzarr)
        return sdata
    else:
        return paths