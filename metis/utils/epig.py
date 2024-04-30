import numpy as np
import logging
import math


def conditional_epig_from_probs(probs_pool, probs_targ):
    """
    See conditional_epig_from_logprobs.

    Arguments:
        probs_pool: ndarray[float], shape [N_p, K, Cl]
        probs_targ: ndarray[float], shape [N_t, K, Cl]

    Returns:
        ndarray[float], shape [N_p, N_t]
    """
    # Estimate the joint predictive distribution.
    probs_pool = np.transpose(probs_pool, (1, 0, 2))  # [K, N_p, Cl]
    probs_targ = np.transpose(probs_targ, (1, 0, 2))  # [K, N_t, Cl]
    probs_pool = np.expand_dims(probs_pool, axis=(2, 4))  # [K, N_p, 1, Cl, 1]
    probs_targ = np.expand_dims(probs_targ, axis=(1, 3))  # [K, 1, N_t, 1, Cl]
    probs_pool_targ_joint = probs_pool * probs_targ
    probs_pool_targ_joint = np.mean(probs_pool_targ_joint, axis=0)

    # Estimate the marginal predictive distributions.
    probs_pool = np.mean(probs_pool, axis=0)
    probs_targ = np.mean(probs_targ, axis=0)

    # Estimate the product of the marginal predictive distributions.
    probs_pool_targ_indep = probs_pool * probs_targ

    # Estimate the conditional expected predictive information gain for each pair of examples.
    # This is the KL divergence between probs_pool_targ_joint and probs_pool_targ_joint_indep.
    nonzero_joint = probs_pool_targ_joint > 0
    log_term = np.copy(probs_pool_targ_joint)
    log_term[nonzero_joint] = np.log(probs_pool_targ_joint[nonzero_joint])
    log_term[nonzero_joint] -= np.log(probs_pool_targ_indep[nonzero_joint])
    scores = np.sum(probs_pool_targ_joint * log_term, axis=(-2, -1))
    return scores  # [N_p, N_t]


def conditional_mse_from_predictions(predictions_pool, predictions_targ):
    """
    Calculate the mean squared error (MSE) between the predicted values for pairs of examples.
    Suitable for regression models.

    Arguments:
        predictions_pool: ndarray[float], shape [N_p]
        predictions_targ: ndarray[float], shape [N_t]

    Returns:
        ndarray[float], shape [N_p, N_t]
    """
    predictions_pool = np.expand_dims(predictions_pool, axis=1)  # [N_p, 1]
    predictions_targ = np.expand_dims(predictions_targ, axis=0)  # [1, N_t]
    mse = np.mean((predictions_pool - predictions_targ) ** 2, axis=-1)  # [N_p, N_t]
    return mse


def check(scores, max_value=math.inf, epsilon=1e-6, score_type=""):
    """
    Warn if any element of scores is negative, a nan or exceeds max_value.

    We set epsilon = 1e-6 based on the fact that np.finfo(np.float).eps ~= 1e-7.
    """
    if not np.all((scores + epsilon >= 0) & (scores - epsilon <= max_value)):
        min_score = np.min(scores)
        max_score = np.max(scores)

        logging.warning(
            f"Invalid {score_type} score (min = {min_score}, max = {max_score})"
        )

    return scores


def epig_from_conditional_scores(scores):
    """
    Arguments:
        scores: ndarray[float], shape [N_p, N_t]

    Returns:
        ndarray[float], shape [N_p,]
    """
    scores = np.mean(scores, axis=-1)  # [N_p,]
    scores = check(scores, score_type="EPIG")  # [N_p,]
    return scores  # [N_p,]


def epig_from_probs(probs_pool, probs_targ, classification=True):
    """
    See epig_from_logprobs.

    Arguments:
        probs_pool: ndarray[float], shape [N_p, K, Cl]
        probs_targ: ndarray[float], shape [N_t, K, Cl]

    Returns:
        ndarray[float], shape [N_p,]
    """
    if classification:
        scores = conditional_epig_from_probs(probs_pool, probs_targ)  # [N_p, N_t]
    else:
        scores = conditional_mse_from_predictions(probs_pool, probs_targ)
    return epig_from_conditional_scores(scores)  # [N_p,]
