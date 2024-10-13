from __future__ import absolute_import, division, print_function
from typing import List, Optional, Tuple, Union, Literal
import os
import subprocess
import logging
import numpy as np
import xarray as xr
from os import PathLike
from gkmSVM import gkmSvmModule


# Verbosity dictionary defines logging level based on integer verbosity
verbosity_dict = {
    0: logging.ERROR,
    1: logging.INFO,
    2: logging.DEBUG,
}


def fit_clf(
    model: gkmSvmModule,
    pos_fasta_path: PathLike,
    neg_fasta_path: PathLike,
    folds: Optional[int] = None,
    fold: Optional[int] = None,
    log_dir: PathLike = "./",
    prefix: str = "seqs",
    overwrite: bool = False,
    seed: Optional[int] = 1,
):
    """Fits a gkmSVM classifier to fasta file inputs

    Returns
    -------
    None
    """
    # Check if model file exists
    if os.path.exists(os.path.join(log_dir, f"{model.model_name}.model.txt")):
        if overwrite:
            os.remove(os.path.join(log_dir, f"{model.model_name}.model.txt"))
        else:
            raise Exception(
                f"Model file {os.path.join(log_dir, f'{model.model_name}')} already exists, set overwrite=True to overwrite"
            )

    # Update gkmSVM model options to include cross-validation
    if folds is not None:
        assert fold is None, "Cannot specify both folds and fold"
        model.options = f"{model.options} -x {int(folds)}"
        model.options = f"{model.options} -r 1"
    
    if fold is not None:
        assert folds is None, "Cannot specify both folds and fold"
        model.options = f"{model.options} -i {int(fold)}"
        model.options = f"{model.options} -r 1"
    
    # Make log directory to house model and log files
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    # Set-up a logger
    logging.basicConfig(
        format="%(asctime)s %(levelname)-8s %(message)s",
        level=verbosity_dict[model.verbosity],
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    logger = logging.getLogger(__name__)

    # Set up log and model file
    log_file_path = os.path.join(log_dir, f"{prefix}_fit.log")
    log_file = open(log_file_path, "w")
    model_file = os.path.join(log_dir, model.model_name)

    # Open the process
    cmd = " ".join(
        [
            "gkmtrain",
            pos_fasta_path,  # positive training file
            neg_fasta_path,  # negative training file
            model_file,  # model file
            model.options,  # options
        ]
    )
    logger.info(f"Running gkmtrain with command: {cmd}")
    process = subprocess.Popen(cmd, stdout=log_file, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = process.communicate()
    if stderr.decode("utf-8") != "":
        err_file = open(os.path.join(log_dir, f"{prefix}_fit.err"), "w")
        err_file.write(stderr.decode("utf-8"))
        raise Exception(
            "Error in gkmtrain, check error fil"
            + os.path.join(log_dir, f"{prefix}_fit.err")
        )
    logger.info(f"Model fit, log file saved to {log_file_path}")
    log_file.close()
    model_file = os.path.abspath(f"{model_file}.model.txt")
    
    # Update model with fit parameters
    if folds is not None:
        model.update_after_fit(model_file)
