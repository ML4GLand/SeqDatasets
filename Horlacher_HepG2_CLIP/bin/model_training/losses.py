from typing import Mapping


import torch
import torch.nn as nn
import torch.nn.functional as F


def MNLLLoss(logps, true_counts):
    """A loss function based on the multinomial negative log-likelihood.

    This loss function takes in a tensor of normalized log probabilities such
    that the sum of each row is equal to 1 (e.g. from a log softmax) and
    an equal sized tensor of true counts and returns the probability of
    observing the true counts given the predicted probabilities under a
    multinomial distribution. Can accept tensors with 2 or more dimensions
    and averages over all except for the last axis, which is the number
    of categories.

    Adapted from Alex Tseng.
    
    Parameters
    ----------
    logps: torch.tensor, shape=(n, ..., L)
        A tensor with `n` examples and `L` possible categories. 

    true_counts: torch.tensor, shape=(n, ..., L)
        A tensor with `n` examples and `L` possible categories.

    Returns
    -------
    loss: float
        The multinomial log likelihood loss of the true counts given the
        predicted probabilities, averaged over all examples and all other
        dimensions.
    """

    log_fact_sum = torch.lgamma(torch.sum(true_counts, dim=-1) + 1)
    log_prod_fact = torch.sum(torch.lgamma(true_counts + 1), dim=-1)
    log_prod_exp = torch.sum(true_counts * logps, dim=-1)
    return -log_fact_sum + log_prod_fact - log_prod_exp


def rbpnet_loss(
    outputs: Mapping[str, torch.Tensor], 
    targets: Mapping[str, torch.Tensor],
    l: int =1.0
):
    y_clip, y_ctl = F.log_softmax(outputs["eCLIP_profile"], dim=-1), F.log_softmax(outputs["control_profile"], dim=-1)
    y_obs, y_obs_ctl = targets["signal"], targets["control"]
    clip_loss = MNLLLoss(y_clip, y_obs)
    ctl_loss = MNLLLoss(y_ctl, y_obs_ctl)
    loss = clip_loss + l*ctl_loss
    return {
        "loss": loss, 
        "eCLIP_loss": clip_loss, 
        "control_loss": ctl_loss
    }

