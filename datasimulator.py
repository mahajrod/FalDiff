#! /usr/bin/env python3

"""
Simulates a quantitative data set and outputs its compositional picture
"""

import random
import numpy as np
import scipy.spatial.distance as distance
import scipy.stats as stats
from collections import namedtuple


class ConnectedComponent:
    """
    Basic components of an AbsoluteSample graph.
    """
    def __init__(self, means, sds, corr_sd=0.2):
        if len(means) != len(sds):
            raise ValueError("Unmatched number of means and standard deviations")

        self._means = np.array(means)
        self._sds = np.array(sds)
        self._size = len(means)

        self._corr_matrix = self._gen_correlation_matrix(corr_sd)
        self._corr_matrix_root = np.linalg.cholesky(self._corr_matrix)

    def __len__(self):
        return self._size

    def _gen_correlation_matrix(self, corr_sd):
        corr_min, corr_max = -1, 1
        corr_mean = 0
        corr_matrix = distance.squareform(stats.truncnorm.rvs(a=corr_min/corr_sd,
                                                              b=corr_max/corr_sd,
                                                              loc=corr_mean, scale=corr_sd,
                                                              size=len(self)))
        corr_matrix += np.identity(len(self))

        return corr_matrix

    def take_sample(self, n, size=500):
        indep_samples = np.vstack([np.ceil(stats.truncnorm.rvs(a=0, b=(5 + mean/sd),
                                                               loc=mean, scale=sd, size=size))
                                   for sd, mean in zip(self._means, self._sds)])

        dep_samples = self._corr_matrix_root.dot(indep_samples)
        return dep_samples[:, random.sample(range(size), n)], np.corrcoef(dep_samples)


def draw_means_from_normal(mean, sd, size):
    return stats.truncnorm.rvs(a=mean/sd, b=(5 + mean/sd), loc=mean, scale=sd, size=size)


def draw_cv_from_normal(mean, sd, size):
    return stats.truncnorm.rvs(a=mean/sd, b=(5 + mean/sd), loc=mean, scale=sd, size=size)


def draw_cv_from_uniform(min, max, size):
    return np.random.uniform(min, max, size)


def cv_to_sd(means, sds):
    return means * sds





class AbsoluteSample(object):
    def __init__(self, size=100, poisson_lam=1):
        pass