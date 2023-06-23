import gzip
import io
import os
from pathlib import Path

import pandas as pd
import wget


def get_download_path():
    """
    Returns a path to download files to
    
    If EUGENe is installed, returns the path specified in settings.dataset_dir
    Otherwise, creates a directory called 'seqdatasets' in the current working directory
    """
    try:
        from eugene import settings
        return settings.dataset_dir
    except:
        return os.path.join(os.getcwd(), 'seqdatasets')
    
def try_download_urls(
    data_idxs: list, 
    url_list: list, 
    dataset_dir: str,
    ds_name: str, 
    processing: str = None
) -> list:
    """Download the data from the given urls.

    Parameters
    ----------
    data_idxs : list
        The indices of the data to be downloaded.
    url_list : list
        The urls of the data to be downloaded.
    ds_name : str
        The name of the dataset to be downloaded.
    processing : str, optional
        The processing of the data to be downloaded. The default is None.
    
    Returns
    -------
    list
        The downloaded data.
    """
    ds_path = os.path.join(dataset_dir, ds_name)
    paths = []
    if processing is not None:
        processing = "." + processing
    for i in data_idxs:
        base_name = os.path.basename(url_list[i]).split("?")[0]
        search_path = os.path.join(dataset_dir, ds_name, base_name)
        if not os.path.exists(search_path):
            if not os.path.isdir(ds_path):
                print(f"Path {ds_path} does not exist, creating new folder.")
                os.makedirs(ds_path)
            print(f"Downloading {ds_name} {os.path.basename(url_list[i])} to {ds_path}...")
            path = wget.download(url_list[i], os.path.relpath(ds_path))
            print(f"Finished downloading {os.path.basename(url_list[i])}")
            if processing == ".gz":
                print("Processing gzip file...")
                with gzip.open(path) as gz:
                    with io.TextIOWrapper(gz, encoding="utf-8") as file:
                        file = pd.read_csv(file, delimiter=r"\t", engine="python", header=None)
                        save_path = os.path.join(ds_path, base_name)
                        print(f"Saving file to {save_path}...")
                        file.to_csv(save_path, index = False, compression="gzip")
                        print(f"Saved file to {save_path}")
                        paths.append(save_path)
            else:
                paths.append(path)
        else:
            print(f"Dataset {ds_name} {base_name} has already been downloaded.")
            paths.append(os.path.join(ds_path, base_name))

    return paths
