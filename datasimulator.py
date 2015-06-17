#! /usr/bin/env python

"""
Simulates a quantitative data set and outputs its compositional picture
"""


from __future__ import division, print_function
import numpy as np
import random


class Component(object):
    """
    Basic component of an AbsoluteSample set.
    """
    def __init__(self, mean_abund, stdev):
        """
        :param mean_abund: mean absolute abundance of the component
        :param stdev: standard deviation of absolute abundance
        :return:
        """
        self._mean = mean_abund
        self._stdev = stdev

    @property
    def mean(self):
        return self._mean

    @property
    def stdev(self):
        return self.stdev()


class AbsoluteSample(object):
    def __init__(self, size=100, poisson_lam=1):
        pass