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
    def __init__(self, means, sds, component_n, corr_sd=0.2):
        if len(means) != len(sds):
            raise ValueError("Unmatched numbers of means and standard deviations")
        self._component_n = component_n

        self._means = np.array(means)
        self._sds = np.array(sds)
        self._size = len(means)

        self._names = ["{}{}".format(component_n, i) for i in range(len(self))]

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

    def _sample_from_distributions(self, means, sds, n, size):
        indep_samples = np.vstack([np.ceil(stats.truncnorm.rvs(a=0, b=(5 + mean/sd),
                                                               loc=mean, scale=sd, size=size))
                                   for sd, mean in zip(means, sds)])

        dep_samples = self._corr_matrix_root.dot(indep_samples)
        return dep_samples[:, random.sample(range(size), n)], np.corrcoef(dep_samples)

    def sample(self, n, size=500):
        return self._sample_from_distributions(self._means, self._sds, n, size)

    def sample_mutated(self, factors, n, size=500):
        if len(factors) > len(self):
            raise ValueError("More factors than components to mutate")
        if not all([0 <= x <= 1 for x in factors]):
            raise ValueError("All factors must be real numbers in [0, 1]")
        iter_factors = iter(factors)

        elements_to_mutate = random.sample(range(len(self)), len(factors))
        mutator_vector = np.fromiter((next(iter_factors) if i in elements_to_mutate else 1
                                      for i in range(len(self))), dtype=np.float64)
        mut_means = np.copy(self._means) * mutator_vector
        mut_sds = np.copy(self._sds) * mutator_vector

        return self._sample_from_distributions(mut_means, mut_sds, n, size), elements_to_mutate

    @property
    def names(self):
        return self._names


class AbsoluteSample:
    def __init__(self, size=100, poisson_lam=1):
        pass


def draw_means_from_normal(mean, sd, size):
    return stats.truncnorm.rvs(a=mean/sd, b=(5 + mean/sd), loc=mean, scale=sd, size=size)


def draw_cv_from_normal(mean, sd, size):
    return stats.truncnorm.rvs(a=mean/sd, b=(5 + mean/sd), loc=mean, scale=sd, size=size)


def draw_cv_from_uniform(min, max, size):
    return np.random.uniform(min, max, size)


def cv_to_sd(means, sds):
    return means * sds


