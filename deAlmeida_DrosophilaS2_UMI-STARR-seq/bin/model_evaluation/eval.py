from os import PathLike
from typing import List, Optional
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from scipy.stats import gaussian_kde
from sklearn.metrics import r2_score, roc_curve, precision_recall_curve
from scipy.stats import pearsonr, spearmanr


def training_summary(
    log_path: PathLike,
    metrics: List[str] = ["r2score"],
    logger: Optional[str] = "tensorboard",
    figsize: Optional[tuple] = (12, 6),
    save: Optional[str] = None,
    return_axes: Optional[bool] = False,
    **kwargs,
) -> None:
    if logger == "csv":
        metrics_df = pd.read_csv(os.path.join(log_path, "metrics.csv"))
        if "step" in metrics_df.columns:
            metrics_df = metrics_df.set_index("step")
            
        _, ax = plt.subplots(1, 2, figsize=figsize)

        if "train_loss" in metrics_df.columns:
            train_loss_step = metrics_df["train_loss"].dropna()
            ax[0].plot(train_loss_step, label="train_loss_step")

        if "train_loss_epoch" in metrics_df.columns:
            train_loss_epoch = metrics_df["train_loss_epoch"].dropna()
            ax[0].plot(train_loss_epoch, label="train_loss_epoch")

        if "val_loss_epoch" in metrics_df.columns:
            val_loss_epoch = metrics_df["val_loss_epoch"].dropna()
            ax[0].plot(val_loss_epoch, label="val_loss_epoch")

        for metric in metrics:
            if f"val_{metric}_epoch" in metrics_df.columns:
                train_metric_epoch = metrics_df[f"train_{metric}_epoch"].dropna().reset_index(drop=True)
                ax[1].plot(train_metric_epoch, label=f"train_{metric}_epoch")

            if f"val_{metric}_epoch" in metrics_df.columns:
                val_metric_epoch = metrics_df[f"val_{metric}_epoch"].dropna().reset_index(drop=True)
                ax[1].plot(val_metric_epoch, label=f"val_{metric}_epoch")

        # Add labels
        ax[0].legend()
        ax[0].set_xlabel("Step")
        ax[0].set_ylabel("Loss")
        ax[1].legend()
        ax[1].set_xlabel("Epoch")
        ax[1].set_ylabel(metric)

    elif logger == "tensorboard":
        raise NotImplementedError("Tensorboard not implemented yet")
    
    else:
        raise ValueError(f"Invalid logger: {logger}")
    
    if save:
        plt.savefig(save)

    if return_axes:
        return ax


def scatter(
    x,
    y,
    ax=None,
    density=False,
    c="b",
    alpha=1,
    s=10,
    xlabel="Observed",
    ylabel="Predicted",
    figsize=(4, 4),
    save=None,
    add_reference_line=True,
    rasterized=False,
):
    # Set up the axes
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Drop NA values if any
    x_nas = np.isnan(x)
    y_nas = np.isnan(y)
    x = x[~x_nas & ~y_nas]
    y = y[~x_nas & ~y_nas]

    if density:
        # Get point densities
        xy = np.vstack([x,y])
        z = gaussian_kde(xy)(xy)

        # Sort the points by density, so that the densest points are plotted last
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]
        c=z

    # Plot the points
    ax.scatter(x, y, c=c, s=s, rasterized=rasterized, alpha=alpha)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    # Add scores
    r2 = r2_score(x, y)
    pearson_r = pearsonr(x, y)
    spearman_r = spearmanr(x, y)
    ax.annotate(f"R2: {r2:.3f}", (0.05, 0.95), xycoords="axes fraction")
    ax.annotate(f"Pearson: {pearson_r[0]:.3f}", (0.05, 0.90), xycoords="axes fraction")
    ax.annotate(f"Spearman: {spearman_r[0]:.3f}", (0.05, 0.85), xycoords="axes fraction")
    
    # Add y=x line for reference but make the mins and maxes extend past the data
    if add_reference_line:
        min_val = min(min(x), min(y))
        max_val = max(max(x), max(y))
        ax.plot([min_val, max_val], [min_val, max_val], c="k", ls="--", lw=1)
    
    # Plt
    plt.tight_layout()

    # Save
    if save:
        plt.savefig(save, dpi=300)
        plt.close()
    else:
        plt.show()


def auroc(
    y_true,
    y_pred,
    ax=None,
    c="b",
    lw=2,
    figsize=(4, 4),
    save=None,
    rasterized=False,
):
    # Set up the axes
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Plot the points
    fpr, tpr, _ = roc_curve(y_true, y_pred)
    ax.plot(fpr, tpr, c=c, lw=lw, rasterized=rasterized)
    ax.plot([0, 1], [0, 1], c="k", ls="--", lw=1)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")

    # Add scores
    auc = np.trapz(tpr, fpr)
    ax.annotate(f"AUC: {auc:.3f}", (0.05, 0.95), xycoords="axes fraction")

    # Plt
    plt.tight_layout()

    # Save
    if save:
        plt.savefig(save, dpi=300)
        plt.close()
    else:
        plt.show()


def auprc(
    y_true,
    y_pred,
    ax=None,
    c="b",
    lw=2,
    figsize=(4, 4),
    save=None,
    rasterized=False,
):
    # Set up the axes
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Plot the points
    precision, recall, _ = precision_recall_curve(y_true, y_pred)
    ax.plot(recall, precision, c=c, lw=lw, rasterized=rasterized)
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")

    # Add scores
    auprc = np.trapz(precision, recall)
    ax.annotate(f"AUPRC: {auprc:.3f}", (0.05, 0.95), xycoords="axes fraction")

    # Plt
    plt.tight_layout()

    # Save
    if save:
        plt.savefig(save, dpi=300)
        plt.close()
    else:
        plt.show()
