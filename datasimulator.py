#! /usr/bin/env python3

"""
Simulates a quantitative data set and outputs its compositional picture
"""

import random
import numpy as np
from scipy.spatial.distance import squareform
import scipy.stats as stats
import scipy.misc as misc
from collections import namedtuple


class ConnectedComponent:
    """
    Basic components of an AbsoluteSample graph.
    """
    def __init__(self, means, sds, component_id, corr_sd=0.2):
        if len(means) != len(sds):
            raise ValueError("Unmatched numbers of means and standard deviations")
        # Number of elements in the component
        self._size = len(means)

        # Component id
        self._component_id = component_id
        self._names = ["{}{}".format(component_id, i) for i in range(len(self))]

        # Distribution data
        self._means = np.array(means)
        self._sds = np.array(sds)

        # Correlation matrix and its Cholesky decomposition
        self._corr_matrix = self._gen_correlation_matrix(corr_sd)
        self._corr_matrix_root = np.linalg.cholesky(self._corr_matrix)

    def __len__(self):
        return self._size

    def _gen_correlation_matrix(self, corr_sd):
        """
        Draws correlation coefficients from a truncated normal distribution with mean = 0 and given
        standard deviation. The numbers drawn belong to the range [-1, 1].
        :param corr_sd: float; standard deviation
        :return: 2D numpy array
        """
        corr_min, corr_max = -1, 1
        corr_mean = 0
        return (squareform(stats.truncnorm.rvs(a=corr_min/corr_sd, b=corr_max/corr_sd,
                                               loc=corr_mean, scale=corr_sd,
                                               size=misc.comb(len(self), 2, exact=True))) +
                np.identity(len(self)))

    def _sample(self, means, sds, n, size):
        """
        Draws n observations from a set of truncated normal distributions representing the
        abundances of elements in a ConnectedComponent. Each row corresponds to an element, hence
        each column represents a sample. The values can't be negative. It uses Cholesky
        decomposition to simulate linear dependency as given by the correlation matrix stored in a
        ConnectedComponent.
        :param means: a vector of mean values
        :param sds: a vector of standard deviations
        :param n: the number of samples to return
        :param size: the size of sampling pool used to simulate correlation; the source of samples
        to be returned; bigger values lead to better correlation simulation and estimation. Can be
        regarded as a pool of biological replicates.
        :return: numpy array with as many columns as samples demanded and a post-hoc estimated
        correlation matrix.
        """
        indep_samples = np.vstack([np.ceil(stats.truncnorm.rvs(a=0, b=(5 + mean/sd),
                                                               loc=mean, scale=sd, size=size))
                                   for sd, mean in zip(means, sds)])

        dep_samples = self._corr_matrix_root.dot(indep_samples)
        return dep_samples[:, random.sample(range(size), n)], np.corrcoef(dep_samples)

    def sample(self, n, size=500):
        """
        Draws n observations from a set of truncated normal distributions representing the
        abundances of elements in a ConnectedComponent. Each row corresponds to an element, hence
        each column represents a sample. The values can't be negative. It uses Cholesky
        decomposition to simulate linear dependency as given by the correlation matrix stored in a
        ConnectedComponent.
        :param n: the number of samples to return
        :param size: the size of sampling pool used to simulate correlation; the source of samples
        to be returned; bigger values lead to better correlation simulation and estimation. Can be
        regarded as a pool of biological replicates.
        :return: numpy array with as many columns as samples demanded and a post-hoc estimated
        correlation matrix.
        """

        return self._sample(self._means, self._sds, n, size)

    def sample_mod(self, factors, n, size=500):
        """
        Copies distribution data of all elements and modifies distribution parameters of as many
        elements as numbers given in a vector of multiplication factors. The elements to be
        modified are picked randomly from a uniform distribution. Draws n observations from
        the resulting set of distributions. Each row corresponds to an element, hence each column
        represents a sample. The values can't be negative. It uses Cholesky
        decomposition to simulate linear dependency as given by the correlation matrix stored in a
        ConnectedComponent.
        :param factors: a vector of multiplication factors (numbers in range [0, +inf))
        :param n: the number of samples to return
        :param size: the size of sampling pool used to simulate correlation; the source of samples
        to be returned; bigger values lead to better correlation simulation and estimation. Can be
        regarded as a pool of biological replicates.
        :return: numpy array with as many columns as samples demanded, a post-hoc estimated
        correlation matrix, indices of elements with modified distributions with corresponding
        multiplication factors
        """
        if len(factors) > len(self):
            raise ValueError("More factors than components to mutate")
        if not all([x >= 0 for x in factors]):
            raise ValueError("All factors must be real numbers in [0, +inf)")
        iter_factors = iter(factors)

        elements_to_mutate = random.sample(range(len(self)), len(factors))
        mutator_vector = np.fromiter((next(iter_factors) if i in elements_to_mutate else 1
                                      for i in range(len(self))), dtype=np.float64)
        mut_means = np.copy(self._means) * mutator_vector
        mut_sds = np.copy(self._sds) * mutator_vector

        return self._sample(mut_means, mut_sds, n, size) + (tuple(zip(elements_to_mutate, factors)), )

    @property
    def names(self):
        """
        Return a list element ids
        :return:
        """
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


def draw_connected_component_size_from_poisson(lam):
    return np.random.poisson(lam) or 1


def cv_to_sd(means, sds):
    return means * sds


